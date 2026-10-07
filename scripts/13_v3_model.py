"""snRNA-VEP v3 (specified before looking at its results): v2 + RNA language model + per-mechanism models.

Changes vs v2 (scripts/11_v2_model.py), fixed in advance:
  1. Separate L2 logistic models per mechanism context (dominant / recessive) instead of one model with interaction
     terms, so the recessive model is not dominated by the many dominant/RNU4-2 rows.
  2. Two extra unsupervised within-gene rank features: RNA-FM zero-shot score (masked marginal for SNVs, sequence
     log-likelihood change otherwise; data/external/rnafm.tsv) and folding destabilisation (ViennaRNA ddG_fold).
Same held-out-unit protocol, same baselines, plus nested selection among {v3, v2, combo, contacts, phyloP, RNA-FM}.
Added after seeing v2 on labels v3 (disclosed): v3-balanced, identical to v3 but with sample weights that give every
(gene unit, class) the same total weight, so the 322 RNU4-2 controls do not dominate the recessive model. It enters
the nested selection as one more candidate.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import metrics, paired_bootstrap_delta

src = open(Path(__file__).resolve().parent / "11_v2_model.py").read().split("results, pooled = [], []")[0]
G = {"__file__": str(Path(__file__).resolve().parent / "11_v2_model.py"), "__name__": "v3"}
exec(compile(src, "11_v2_model", "exec"), G)
d, rows_for_training, design, model, R = G["d"], G["rows_for_training"], G["design"], G["model"], G["R"]

fm = pd.read_csv(config.ROOT / "data/external/rnafm.tsv", sep="\t").drop_duplicates("key").set_index("key")
d["rnafm"] = -d.key.map(fm.rnafm_mm).fillna(d.key.map(fm.rnafm_dll))     # higher = more damaging
d["rnafm_mm"] = -d.key.map(fm.rnafm_mm)
d["destab"] = d.ddG_fold
for c in ["rnafm", "destab"]:
    d[c + "_r"] = d.groupby("gene_name")[c].rank(pct=True)
R3 = R + ["rnafm_r", "destab_r"]
UNITS = ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6", "RNU4ATAC", "RNU12", "RNU6ATAC"]


def ctx_rows(df, unit, ctx):
    pref = "AD-" if ctx == "dominant" else "AR-"
    return df[(df.unit == unit) & ((df.y_path == 0) | ((df.y_path == 1) & df["mode"].str.contains(pref)))]


def v3_scores(train_df, test_df, ctx, balanced=False):
    tr = rows_for_training(train_df)
    tr = tr[tr.ctx == ctx]
    X = d.loc[tr.i, R3].fillna(0.5).to_numpy()
    w = None
    if balanced:
        grp = d.loc[tr.i, "unit"].to_numpy() + "|" + tr.y.astype(str).to_numpy()
        cnt = pd.Series(grp).value_counts()
        w = (1.0 / cnt[grp].to_numpy()) * len(grp) / cnt.size
    m = model()
    m.fit(X, tr.y, logisticregression__sample_weight=w)
    return m.predict_proba(test_df[R3].fillna(0.5).to_numpy())[:, 1], m


def all_scores(train_df, test_df, ctx):
    out = test_df.copy()
    out["snrnavep_v3"], _ = v3_scores(train_df, test_df, ctx)
    out["snrnavep_v3b"], _ = v3_scores(train_df, test_df, ctx, balanced=True)
    tr = rows_for_training(train_df)
    m2 = model().fit(design(d, tr.i, tr.ctx), tr.y)
    out["snrnavep_v2"] = m2.predict_proba(design(d, test_df.index, [ctx] * len(test_df)))[:, 1]
    out["combo_untrained"] = test_df[["frac_states_protein_contact_r", "phylop447_r"]].mean(1)
    return out


SC = ["snrnavep_v3b", "snrnavep_v3", "snrnavep_v2", "combo_untrained", "frac_states_protein_contact", "phylop447", "rnafm", "cadd_phred"]
rows, pooled, nested = [], [], []
for unit in UNITS:
    for ctx in ["dominant", "recessive"]:
        test = ctx_rows(d, unit, ctx)
        if (test.y_path == 1).sum() < 3 or (test.y_path == 0).sum() < 3:
            continue
        t = all_scores(d[d.unit != unit], test, ctx)
        for s in SC:
            rows.append(dict(unit=unit, ctx=ctx, score=s, auroc=metrics(t.y_path, t[s])["auroc"], n_pos=int(t.y_path.sum())))
            t[s + "_rk"] = t[s].rank(pct=True)
        pooled.append(t.assign(ctx=ctx))
        # nested choice among candidates (inner leave-one-unit-out over training units only)
        cand = ["snrnavep_v3b", "snrnavep_v3", "snrnavep_v2", "combo_untrained", "frac_states_protein_contact", "phylop447", "rnafm"]
        inner = {c: [] for c in cand}
        for v in [u for u in UNITS if u != unit]:
            vt = ctx_rows(d, v, ctx)
            if (vt.y_path == 1).sum() < 3 or (vt.y_path == 0).sum() < 3:
                continue
            sv = all_scores(d[(d.unit != unit) & (d.unit != v)], vt, ctx)
            for c in cand:
                inner[c].append(metrics(sv.y_path, sv[c])["auroc"])
        choice = max(cand, key=lambda c: np.nanmean(inner[c]) if inner[c] else -1)
        nested.append(dict(unit=unit, ctx=ctx, chosen=choice, auroc=metrics(t.y_path, t[choice])["auroc"], n_pos=int(t.y_path.sum())))

res = pd.DataFrame(rows)
res.to_csv(config.result("v3_heldout.tsv"), sep="\t", index=False)
print(res.pivot_table(index=["ctx", "score"], columns="unit", values="auroc").round(3).to_string())
print("\nmean AUROC across held-out units:")
print(res.groupby(["ctx", "score"]).auroc.mean().unstack(0).round(3).to_string())
P = pd.concat(pooled)
P.to_csv(config.result("v3_heldout_predictions.tsv.gz"), sep="\t", index=False)
print("\nrank-pooled paired bootstrap:")
for ctxs in [["dominant", "recessive"], ["dominant"], ["recessive"]]:
    Q = P[P.ctx.isin(ctxs)]
    for a, b in [("snrnavep_v3b", "snrnavep_v3"), ("snrnavep_v3b", "phylop447"), ("snrnavep_v3b", "combo_untrained"), ("snrnavep_v3", "snrnavep_v2"), ("snrnavep_v3", "phylop447"), ("snrnavep_v3", "rnafm"),
                 ("snrnavep_v3", "combo_untrained")]:
        dl, lo, hi, p = paired_bootstrap_delta(Q.y_path.values, Q[a + "_rk"].values, Q[b + "_rk"].values, n=2000)
        print(f"  [{'+'.join(ctxs)}] {a} - {b:28s}: {dl:+.3f} [{lo:+.3f}, {hi:+.3f}] p {p:.3f}")
    s = Q[Q.cadd_phred.notna()]
    dl, lo, hi, p = paired_bootstrap_delta(s.y_path.values, s.snrnavep_v3_rk.values, s.cadd_phred_rk.values, n=2000)
    print(f"  [{'+'.join(ctxs)}] snrnavep_v3 - CADD (SNVs n {len(s)}): {dl:+.3f} [{lo:+.3f}, {hi:+.3f}] p {p:.3f}")
N = pd.DataFrame(nested)
N.to_csv(config.result("v3_nested.tsv"), sep="\t", index=False)
print("\n=== nested ===")
print(N.round(3).to_string(index=False))
print(N.groupby("ctx").auroc.mean().round(3).to_string())

# SGE zero-shot of v3 dominant model and RNA-FM
test = d[(d.gene_name == "RNU4-2") & d.y_sge.notna()]
p, m = v3_scores(d[d.unit != "RNU4-2"], test, "dominant")
print(f"\nSGE zero-shot: v3 (no RNU4-2) rho {spearmanr(p, test.y_sge)[0]:.3f}; RNA-FM rho "
      f"{spearmanr(test.rnafm, test.y_sge, nan_policy='omit')[0]:.3f}; RNA-FM masked marginal (SNVs) rho "
      f"{spearmanr(test.rnafm_mm, test.y_sge, nan_policy='omit')[0]:.3f}")
for ctx in ["dominant", "recessive"]:
    for b in [False, True]:
        _, m = v3_scores(d, d.iloc[:5], ctx, balanced=b)
        print(f"final {ctx} {'balanced' if b else ''} coefficients:", dict(zip(R3, np.round(m[-1].coef_[0], 2))))
