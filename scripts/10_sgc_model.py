"""snRNA-VEP v1 model: Simplified Graph Convolution (SGC; Wu et al. 2019) over the multi-state spliceosome contact
graph, then L2-regularized logistic regression. Low capacity by design (labels: ~90 pathogenic variants in 7 genes).

  node features H0 = family-agnostic structure summaries (snrna_vep.graph.node_feature_matrix_agnostic)
  H_k = A_norm^k H0 for k = 1, 2   (A = RNA-RNA contacts in any cryo-EM state + backbone + self loops)
  variant x = [H0, H1, H2 at the affected position, variant type]   -> logistic regression

Evaluation as in scripts/09_train_net.py (gene-out / family-out; RNU6 RP genes as one unit; identical test sets),
plus SGE zero-shot. Pooled comparisons use within-unit percentile ranks so separately trained models are comparable.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import bootstrap_auroc, metrics, paired_bootstrap_delta
from snrna_vep.net import build_graph

src = open(Path(__file__).resolve().parent / "09_train_net.py").read().split("ARGS = dict(")[0]
G = {"__file__": str(Path(__file__).resolve().parent / "09_train_net.py"), "__name__": "sgc"}
exec(compile(src, "09_train_net", "exec"), G)
d, idx, X, A, FAMILY = G["d"], G["idx"], G["X"], G["A"], G["FAMILY"]


def unit_mask(df, unit):
    return df.gene_name.str.startswith("RNU6-") if unit == "RNU6" else (df.gene_name == unit)

# SGC propagation
H = [X]
for _ in range(2):
    H.append(torch.sparse.mm(A, H[-1]))
Hcat = torch.cat(H, 1).numpy()
node_i = np.array([idx[n] for n in d.node_name])
type_cols = d[["is_snv", "is_ins", "is_del"]].astype(float).to_numpy()
XV = np.concatenate([Hcat[node_i], type_cols], 1)


def model(C=0.05):
    return make_pipeline(StandardScaler(), LogisticRegression(C=C, class_weight="balanced", max_iter=5000))


# SGE pseudo-labels from the SGE paper's own calibrated ranges (PS3: score <= -0.46; BS3: score >= -0.27)
sge_raw = d.sge_score if "sge_score" in d else pd.Series(np.nan, index=d.index)
y_sgelab = np.where(sge_raw <= -0.46, 1.0, np.where(sge_raw >= -0.27, 0.0, np.nan))


def fit_pred(train_mask, test_mask, cols=None, use_sge=False):
    y = d.y_path.values.copy()
    if use_sge:
        extra = train_mask & np.isnan(y) & ~np.isnan(y_sgelab)
        y = np.where(extra, y_sgelab, y)
    tr = train_mask & ~np.isnan(y)
    Xtr, Xte = XV[tr], XV[test_mask]
    if cols is not None:
        Xtr, Xte = Xtr[:, cols], Xte[:, cols]
    m = model().fit(Xtr, y[tr].astype(int))
    return m.predict_proba(Xte)[:, 1]


k0 = X.shape[1]
COLSETS = {"sgc_k2": None, "no_graph_k0": list(range(k0)) + list(range(3 * k0, 3 * k0 + 3))}
rows, pooled = [], []
for ev in ["gene-out", "family-out"]:
    for unit in ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6", "RNU4ATAC", "RNU12"]:
        test = (unit_mask(d, unit) & d.y_path.notna()).values
        if d.y_path.values[test].sum() < 3:
            continue
        drop = unit_mask(d, unit).values if ev == "gene-out" else (d.paralog_family == FAMILY[unit]).values
        t = d[test].copy()
        for name, cols in COLSETS.items():
            t[name] = fit_pred(~drop, test, cols)
        t["sgc_k2_sge"] = fit_pred(~drop, test, None, use_sge=(unit != "RNU4-2"))
        for s in list(COLSETS) + ["sgc_k2_sge", "frac_states_protein_contact", "max_protein_res", "cadd_phred", "phylop447"]:
            rows.append(dict(eval=ev, unit=unit, score=s, **metrics(t.y_path, t[s])))
        # within-unit percentile ranks for pooling
        for s in list(COLSETS) + ["sgc_k2_sge", "frac_states_protein_contact", "cadd_phred", "phylop447"]:
            t[s + "_r"] = t[s].rank(pct=True)
        pooled.append(t.assign(eval_mode=ev, unit=unit))

res = pd.DataFrame(rows)
res.to_csv(config.RESULTS / "sgc_heldout.tsv", sep="\t", index=False)
print("=== Held-out AUROC (identical test variants) ===")
print(res.pivot_table(index="score", columns=["eval", "unit"], values="auroc").round(3).to_string())
P = pd.concat(pooled)
for ev in ["gene-out", "family-out"]:
    t = P[P.eval_mode == ev]
    mean_auc = res[(res["eval"] == ev)].groupby("score").auroc.mean().round(3).to_dict()
    print(f"\n[{ev}] mean AUROC across units: {mean_auc}")
    for a, b in [("sgc_k2_sge", "sgc_k2"), ("sgc_k2_sge", "frac_states_protein_contact"), ("sgc_k2", "frac_states_protein_contact"), ("sgc_k2", "no_graph_k0"), ("sgc_k2", "phylop447"),
                 ("frac_states_protein_contact", "phylop447")]:
        dl, lo, hi, pv = paired_bootstrap_delta(t.y_path.values, t[a + "_r"].values, t[b + "_r"].values, n=2000)
        print(f"   {a} - {b} (rank-pooled): {dl:+.3f} [{lo:+.3f}, {hi:+.3f}] p {pv:.3f}")
    s = t[t.cadd_phred.notna()]
    for a in ["sgc_k2", "sgc_k2_sge", "frac_states_protein_contact"]:
        dl, lo, hi, pv = paired_bootstrap_delta(s.y_path.values, s[a + "_r"].values, s.cadd_phred_r.values, n=2000)
        print(f"   {a} - CADD (SNVs n {len(s)}, path {int(s.y_path.sum())}): {dl:+.3f} [{lo:+.3f}, {hi:+.3f}] p {pv:.3f}")

print("\n=== RNU4-2 SGE zero-shot ===")
test = ((d.gene_name == "RNU4-2") & d.y_sge.notna()).values
for ev, drop in [("gene-out", (d.gene_name == "RNU4-2").values), ("family-out", (d.paralog_family == "U4").values)]:
    for name, cols in COLSETS.items():
        p = fit_pred(~drop, test, cols)
        print(f"  {ev:10s} {name:12s} spearman {spearmanr(p, d.y_sge.values[test])[0]:+.3f}")
tt = d[test]
for s in ["frac_states_protein_contact", "cadd_phred", "phylop447"]:
    x = tt[[s, "y_sge"]].dropna()
    print(f"  baseline   {s:28s} spearman {spearmanr(x[s], x.y_sge)[0]:+.3f}")
