"""snRNA-VEP v2: mechanism-aware, rank-normalised evidence model.

Rationale (from scripts/10 and the mechanism analysis): dominant snRNA variants (e.g. ReNU, altered function) are
invisible to conservation but visible to cryo-EM contact profiles; recessive loss-of-function variants sit at
conserved positions. Evidence is therefore combined per disease mechanism, as ACMG/ClinGen interpretation requires.

Features (unsupervised, computed over ALL possible variants of each gene, so no labels are used): within-gene
percentile ranks of
  contact   : fraction of resolved splicing states with protein contact; max protein contacts; snRNA-snRNA contacts
  graph     : the same contact fraction averaged over 3D neighbours in the multi-state contact graph (1 hop)
  conserv.  : phyloP447
Context    : the disease mechanism the variant is being assessed for (dominant / recessive). Training rows: patient
             variants carry their own mode; controls of genes with both modes appear once per context.
Model      : L2 logistic regression on [features, features x recessive-context] (11 coefficients).
Evaluation : held-out unit (gene; RNU6 RP genes together), separately for dominant and recessive positives against
             the unit's controls; mean AUROC across units; rank-pooled paired bootstrap vs baselines.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import bootstrap_auroc, metrics, paired_bootstrap_delta

src = open(Path(__file__).resolve().parent / "09_train_net.py").read().split("ARGS = dict(")[0]
G = {"__file__": str(Path(__file__).resolve().parent / "09_train_net.py"), "__name__": "v2"}
exec(compile(src, "09_train_net", "exec"), G)
d, idx, X, A, nodes_f = G["d"].copy(), G["idx"], G["X"], G["A"], G["nodes_f"]

L2 = pd.read_csv(config.CURATION / config.LABELS, sep="\t")
modes = L2[L2.label == 1].groupby("key")["mode"].agg(lambda s: "|".join(sorted(set(s))))
d["mode"] = d.key.map(modes).fillna("")
d["unit"] = np.where(d.gene_name.str.startswith("RNU6-"), "RNU6", d.gene_name)

# graph-smoothed contact fraction (1 hop over the contact graph, row-normalised adjacency incl. self loop)
fp = torch.tensor(nodes_f.set_index(nodes_f.family + ":" + nodes_f.ref_pos.astype(int).astype(str))
                  .reindex(list(idx)).frac_prot.fillna(0).to_numpy(), dtype=torch.float32).unsqueeze(1)
nb = torch.sparse.mm(A, fp).squeeze(1).numpy()
d["nb_frac_prot"] = nb[[idx[n] for n in d.node_name]]

FEATS = ["frac_states_protein_contact", "max_protein_res", "frac_states_snrna_contact", "nb_frac_prot", "phylop447"]
for c in FEATS:                                                    # unsupervised within-gene ranks (all variants)
    d[c + "_r"] = d.groupby("gene_name")[c].rank(pct=True)
R = [c + "_r" for c in FEATS]

gene_modes = d[d.y_path == 1].groupby("gene_name")["mode"].agg(lambda s: set("|".join(s).split("|")))


def rows_for_training(df):
    """Expand labelled rows into (features, context) training rows."""
    out = []
    for r in df[df.y_path.notna()].itertuples():
        if r.y_path == 1:
            ctxs = {("recessive" if m.startswith("AR") else "dominant") for m in r.mode.split("|") if m}
        else:
            gm = gene_modes.get(r.gene_name, set())
            ctxs = {("recessive" if m.startswith("AR") else "dominant") for m in gm} or {"dominant"}
        for c in ctxs:
            out.append((r.Index, c, r.y_path))
    return pd.DataFrame(out, columns=["i", "ctx", "y"])


def design(df, idxs, ctx):
    F = df.loc[idxs, R].to_numpy()
    rec = (np.asarray(ctx) == "recessive").astype(float)[:, None]
    return np.concatenate([F, F * rec, rec], 1)


def model():
    return make_pipeline(StandardScaler(), LogisticRegression(C=0.3, class_weight="balanced", max_iter=5000))


results, pooled = [], []
for unit in ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6", "RNU4ATAC", "RNU12", "RNU6ATAC"]:
    train = d[d.unit != unit]
    tr = rows_for_training(train)
    m = model().fit(design(d, tr.i, tr.ctx), tr.y)
    for ctx, pref in [("dominant", "AD-"), ("recessive", "AR-")]:
        t = d[(d.unit == unit) & ((d.y_path == 0) | ((d.y_path == 1) & d["mode"].str.contains(pref)))].copy()
        if (t.y_path == 1).sum() < 3 or (t.y_path == 0).sum() < 3:
            continue
        t["snrnavep_v2"] = m.predict_proba(design(d, t.index, [ctx] * len(t)))[:, 1]
        t["combo_untrained"] = t[["frac_states_protein_contact_r", "phylop447_r"]].mean(1)
        for s in ["snrnavep_v2", "combo_untrained", "frac_states_protein_contact", "phylop447", "cadd_phred"]:
            mm = metrics(t.y_path, t[s])
            results.append(dict(unit=unit, ctx=ctx, score=s, auroc=mm["auroc"], n_pos=mm["n_pos"], n=mm["n"]))
            t[s + "_rk"] = t[s].rank(pct=True)
        pooled.append(t.assign(ctx=ctx))

res = pd.DataFrame(results)
res.to_csv(config.result("v2_heldout.tsv"), sep="\t", index=False)
print(res.pivot_table(index=["ctx", "score"], columns="unit", values="auroc").round(3).to_string())
print("\nmean AUROC across held-out units:")
print(res.groupby(["ctx", "score"]).auroc.mean().unstack(0).round(3).to_string())
P = pd.concat(pooled)
print("\nrank-pooled paired bootstrap (all held-out units, both contexts):")
for b in ["combo_untrained", "frac_states_protein_contact", "phylop447"]:
    dl, lo, hi, p = paired_bootstrap_delta(P.y_path.values, P.snrnavep_v2_rk.values, P[b + "_rk"].values, n=2000)
    print(f"  v2 - {b:28s}: {dl:+.3f} [{lo:+.3f}, {hi:+.3f}] p {p:.3f}")
s = P[P.cadd_phred.notna()]
dl, lo, hi, p = paired_bootstrap_delta(s.y_path.values, s.snrnavep_v2_rk.values, s.cadd_phred_rk.values, n=2000)
print(f"  v2 - CADD (SNVs, n {len(s)}, pathogenic {int(s.y_path.sum())}): {dl:+.3f} [{lo:+.3f}, {hi:+.3f}] p {p:.3f}")
lo, hi = bootstrap_auroc(P.y_path, P.snrnavep_v2_rk, n=1000)
print(f"  v2 rank-pooled AUROC {metrics(P.y_path, P.snrnavep_v2_rk)['auroc']:.3f} [{lo:.3f}, {hi:.3f}]")
P.to_csv(config.result("v2_heldout_predictions.tsv.gz"), sep="\t", index=False)

# final model on all data: coefficients for interpretation
tr = rows_for_training(d)
mf = model().fit(design(d, tr.i, tr.ctx), tr.y)
coef = mf[-1].coef_[0]
names = R + [c + " x recessive" for c in R] + ["recessive context"]
print("\nfinal-model coefficients (standardised):")
for n, c in zip(names, coef):
    print(f"  {n:42s} {c:+.3f}")


# ---------------------------------------------------------------------------------------------------------------
# Nested selection (snRNA-VEP v2-nested): for each held-out unit and context, choose among candidate scorers using
# an inner leave-one-unit-out loop over the remaining units only, then apply the winner to the held-out unit.
CANDIDATES = ["snrnavep_v2", "combo_untrained", "frac_states_protein_contact", "phylop447"]
UNITS = ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6", "RNU4ATAC", "RNU12", "RNU6ATAC"]


def scores_for(train_df, test_df, ctx):
    tr = rows_for_training(train_df)
    m = model().fit(design(d, tr.i, tr.ctx), tr.y)
    out = test_df.copy()
    out["snrnavep_v2"] = m.predict_proba(design(d, test_df.index, [ctx] * len(test_df)))[:, 1]
    out["combo_untrained"] = test_df[["frac_states_protein_contact_r", "phylop447_r"]].mean(1)
    return out


def ctx_rows(df, unit, ctx):
    pref = "AD-" if ctx == "dominant" else "AR-"
    return df[(df.unit == unit) & ((df.y_path == 0) | ((df.y_path == 1) & df["mode"].str.contains(pref)))]


nested = []
for unit in UNITS:
    rest = [u for u in UNITS if u != unit]
    for ctx in ["dominant", "recessive"]:
        test = ctx_rows(d, unit, ctx)
        if (test.y_path == 1).sum() < 3 or (test.y_path == 0).sum() < 3:
            continue
        inner = {c: [] for c in CANDIDATES}
        for v in rest:
            vt = ctx_rows(d, v, ctx)
            if (vt.y_path == 1).sum() < 3 or (vt.y_path == 0).sum() < 3:
                continue
            sc = scores_for(d[(d.unit != unit) & (d.unit != v)], vt, ctx)
            for c in CANDIDATES:
                inner[c].append(metrics(sc.y_path, sc[c])["auroc"])
        choice = max(CANDIDATES, key=lambda c: np.nanmean(inner[c]) if inner[c] else -1)
        sc = scores_for(d[d.unit != unit], test, ctx)
        nested.append(dict(unit=unit, ctx=ctx, chosen=choice, auroc=metrics(sc.y_path, sc[choice])["auroc"],
                           n_pos=int(sc.y_path.sum())))
N = pd.DataFrame(nested)
N.to_csv(config.result("v2_nested.tsv"), sep="\t", index=False)
print("\n=== snRNA-VEP v2-nested (scorer chosen per context by inner leave-one-unit-out on training units only) ===")
print(N.round(3).to_string(index=False))
print(N.groupby("ctx").auroc.mean().round(3).to_string())
