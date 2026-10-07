"""Reviewer-driven robustness analyses (labels from config.LABELS).

(a) Clustered uncertainty for the final system (v2-nested) vs baselines: (i) position-block bootstrap (resample snRNA
    positions with replacement within each gene unit; all variants at a position move together), (ii) gene-level
    bootstrap of the mean per-unit ΔAUROC, (iii) per-unit sign consistency.
(b) Stronger baselines under the identical held-out protocol: secondary structure (ViennaRNA ΔΔG of folding,
    pairing probability, U4/U6 or partner duplex ΔΔG), local constraint (gnomAD observed/expected in 10-nt windows;
    circular where controls come from gnomAD, flagged), and a TRAINED secondary-structure model (same design and
    regularisation as v2, structural contact features replaced by the ViennaRNA features).
(c) Within-gene dominant vs recessive: in RNU4-2 and RNU2-2 (genes with both modes), compare dominant and recessive
    pathogenic variants on structural contact and conservation (logistic regression, position-clustered SEs), and
    their contacts with named proteins.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import metrics

src = open(Path(__file__).resolve().parent / "11_v2_model.py").read().split("results, pooled = [], []")[0]
G = {"__file__": str(Path(__file__).resolve().parent / "11_v2_model.py"), "__name__": "rob"}
exec(compile(src, "11_v2_model", "exec"), G)
d, rows_for_training, design, model, R = G["d"], G["rows_for_training"], G["design"], G["model"], G["R"]
rng = np.random.default_rng(0)
UNITS = ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6", "RNU4ATAC", "RNU12"]

# ---- (b) baselines: unsupervised ranks + trained secondary-structure model
d["vienna_destab"] = d.ddG_fold
d["vienna_paired"] = 1 - d.p_unpaired_wt
d["vienna_duplex"] = d.ddG_duplex
d["constraint_1_minus_oe"] = 1 - d.rel_oe_w10
SS = ["vienna_destab", "vienna_paired", "vienna_duplex", "phylop447"]
for c in SS + ["constraint_1_minus_oe"]:
    d[c + "_r"] = d.groupby("gene_name")[c].rank(pct=True).fillna(0.5)
saved = list(R)


def fit_scores(cols, train_df, test_df, ctx):
    R[:] = cols
    tr = rows_for_training(train_df)
    m = model().fit(design(d, tr.i, tr.ctx), tr.y)
    p = m.predict_proba(design(d, test_df.index, [ctx] * len(test_df)))[:, 1]
    R[:] = saved
    return p


Pv = pd.read_csv(config.result("v2_heldout_predictions.tsv.gz"), sep="\t", low_memory=False)
N = pd.read_csv(config.result("v2_nested.tsv"), sep="\t")
choice = {(u, c): s for u, c, s in zip(N.unit, N.ctx, N.chosen)}
rows, pooled = [], []
for unit in UNITS:
    for ctx, pref in [("dominant", "AD-"), ("recessive", "AR-")]:
        t = d[(d.unit == unit) & ((d.y_path == 0) | ((d.y_path == 1) & d["mode"].str.contains(pref)))].copy()
        if (t.y_path == 1).sum() < 3 or (t.y_path == 0).sum() < 3:
            continue
        pv = Pv[(Pv.unit == unit) & (Pv.ctx == ctx)].set_index("key")
        t["final"] = t.key.map(pv[choice[(unit, ctx)]])
        t["trained_secondary_structure"] = fit_scores([c + "_r" for c in SS], d[d.unit != unit], t, ctx)
        for s in ["final", "trained_secondary_structure", "vienna_destab", "vienna_paired", "vienna_duplex",
                  "constraint_1_minus_oe", "cadd_phred", "phylop447", "frac_states_protein_contact"]:
            rows.append(dict(unit=unit, ctx=ctx, score=s, auroc=metrics(t.y_path, t[s])["auroc"]))
        pooled.append(t.assign(ctx=ctx))
B = pd.DataFrame(rows)
B.to_csv(config.result("baselines_extended.tsv"), sep="\t", index=False)
print("=== (b) held-out AUROC, extended baselines ===")
print(B.pivot_table(index="score", columns=["ctx", "unit"], values="auroc").round(3).to_string())
print(B.groupby(["ctx", "score"]).auroc.mean().unstack(0).round(3).sort_values("dominant", ascending=False).to_string())

# ---- (a) clustered uncertainty
P = pd.concat(pooled)


def auc_units(df, s):
    return {u: metrics(g.y_path, g[s])["auroc"] for u, g in df.groupby("unit")}


print("\n=== (a) clustered uncertainty: final minus baseline ===")
out = []
for ctx in ["dominant", "recessive"]:
    Q = P[P.ctx == ctx]
    for b in ["cadd_phred", "phylop447", "frac_states_protein_contact", "trained_secondary_structure"]:
        Qb = Q[Q[b].notna()]
        obs = {u: auc_units(Qb, "final")[u] - auc_units(Qb, b)[u] for u in Qb.unit.unique()}
        point = np.mean(list(obs.values()))
        # position-block bootstrap within units
        bs = []
        for _ in range(1000):
            parts = []
            for u, g in Qb.groupby("unit"):
                pos = g.node_name.unique()
                pick = rng.choice(pos, len(pos))
                parts.append(pd.concat([g[g.node_name == p_] for p_ in pick]))
            S = pd.concat(parts)
            try:
                bs.append(np.mean([auc_units(S, "final")[u] - auc_units(S, b)[u] for u in S.unit.unique()]))
            except Exception:
                continue
        bs = np.array([x for x in bs if np.isfinite(x)])
        # gene-level bootstrap of per-unit deltas
        dl = np.array(list(obs.values()))
        gb = np.array([rng.choice(dl, len(dl)).mean() for _ in range(5000)])
        out.append(dict(ctx=ctx, baseline=b, mean_delta=point, block_ci=f"[{np.quantile(bs, .025):+.3f}, {np.quantile(bs, .975):+.3f}]",
                        block_p=2 * min((bs <= 0).mean(), (bs >= 0).mean()), gene_ci=f"[{np.quantile(gb, .025):+.3f}, {np.quantile(gb, .975):+.3f}]",
                        units_better=f"{(dl > 0).sum()}/{len(dl)}"))
A = pd.DataFrame(out)
A.to_csv(config.result("clustered_tests.tsv"), sep="\t", index=False)
print(A.round(3).to_string(index=False))

# ---- (c) within-gene dominant vs recessive
print("\n=== (c) within-gene: dominant vs recessive pathogenic variants ===")
named = pd.read_csv(config.result("mechanism_proteins.tsv"), sep="\t") if config.result("mechanism_proteins.tsv").exists() else None
rowsc = []
for gene in ["RNU4-2", "RNU2-2P"]:
    pth = d[(d.gene_name == gene) & (d.y_path == 1)].copy()
    pth["dom"] = pth["mode"].str.contains("AD-").astype(int)
    pth = pth[pth["mode"].str.contains("AD-") ^ pth["mode"].str.contains("AR-")]       # unambiguous mode only
    X = sm.add_constant(pth[["frac_states_protein_contact_r", "nb_frac_prot_r", "phylop447_r"]])
    fit = sm.GLM(pth.dom, X, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": pth.node_name.astype("category").cat.codes})
    for c in X.columns[1:]:
        rowsc.append(dict(gene=gene, n_dom=int(pth.dom.sum()), n_rec=int((1 - pth.dom).sum()), feature=c,
                          coef=fit.params[c], p_cluster=fit.pvalues[c],
                          median_dom=pth.loc[pth.dom == 1, c].median(), median_rec=pth.loc[pth.dom == 0, c].median()))
Cw = pd.DataFrame(rowsc)
Cw.to_csv(config.result("within_gene_dom_vs_rec.tsv"), sep="\t", index=False)
print(Cw.round(3).to_string(index=False))
