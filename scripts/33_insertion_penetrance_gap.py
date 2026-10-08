"""Does a saturation functional assay explain the insertion/SNV penetrance gap in the RNU4-2 critical region?

Chen et al. (Nature 2024) report, for the same 18-nt critical region, single-base insertions at OR 1,531 in NDD
(54/8,841 vs 2/490,132) but SNVs at OR 9.51 (6/8,841 vs 35/490,132) - a ~160-fold difference in effect size.
The saturation genome editing assay of De Jonghe et al. (Nature 2026) scores the two classes comparably.
This script quantifies that discrepancy per variant.

Data: data/curation/renu_variant_counts.tsv (Chen Extended Data Table 1, per-variant carrier counts in NDD and
population cohorts) joined to the SGE function scores and the multi-state structural features.

Analyses
  1. Class-level: SGE score distributions for insertions vs SNVs in the critical region (Mann-Whitney), against the
     published odds ratios - the discrepancy itself.
  2. Variant-level: logistic regression of "observed in NDD" on SGE score, with and without an insertion term,
     over critical-region variants that the SGE assay scored. If the insertion term stays large and significant
     after conditioning on the assay score, the assay does not explain the excess.
  3. Case-count model: negative-binomial / Poisson regression of NDD carrier count on SGE score + insertion.
  4. Matched comparison: SNVs whose SGE score is at least as damaging as n.64_65insT - how many are seen in patients?
     This is the cleanest statement of the gap and needs no model.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

CR = (120291825, 120291842)
NDD_N = 8841          # GEL NDD cohort
POP_N = 490132 + 76215 + 245400   # UK Biobank + gnomAD v4.0 + All of Us

C = pd.read_csv(config.CURATION / "renu_variant_counts.tsv", sep="\t", comment="#")
C["ndd"] = C.gel_ndd + C.nongel_ndd
C["pop"] = C.gnomad + C.ukb + C.aou

feat = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t", low_memory=False)
sge = pd.read_csv(config.DATA / "external" / "sge_rnu4-2.tsv", sep="\t")
F = feat[feat.gene_name == "RNU4-2"].drop(columns=["sge_score", "sge_type"], errors="ignore").merge(sge, on="key", how="left")
F = F[(F.pos >= CR[0]) & (F.pos <= CR[1])].copy()
D = F.merge(C[["key", "ndd", "pop", "gel_ndd"]], on="key", how="left").fillna({"ndd": 0, "pop": 0, "gel_ndd": 0})
D = D[D.sge_score.notna() & D.vtype.isin(["ins", "snv"])].copy()
D["is_ins"] = (D.vtype == "ins").astype(int)
D["seen_ndd"] = (D.ndd > 0).astype(int)
print(f"critical-region variants with an SGE score: {len(D)} ({int(D.is_ins.sum())} insertions, {int((1 - D.is_ins).sum())} SNVs)")
print(f"  observed in NDD: {int(D.seen_ndd.sum())} ({int(D[D.is_ins == 1].seen_ndd.sum())} ins, {int(D[D.is_ins == 0].seen_ndd.sum())} snv)")

# ---- 1. the discrepancy
print("\n=== 1. class-level: does the assay separate the two classes? ===")
a, b = D[D.is_ins == 1].sge_score, D[D.is_ins == 0].sge_score
u, p = mannwhitneyu(a, b)
print(f"  SGE score, insertions n={len(a)} median {a.median():+.3f} vs SNVs n={len(b)} median {b.median():+.3f}: Mann-Whitney p={p:.3f}")
print(f"  published human odds ratios for the same region: insertions 1,531 (p=3.3e-92); SNVs 9.51 (p=8.2e-5)")
ct = [[int(D[D.is_ins == 1].ndd.sum()), int(D[D.is_ins == 1]["pop"].sum())],
      [int(D[D.is_ins == 0].ndd.sum()), int(D[D.is_ins == 0]["pop"].sum())]]
orr, pf = fisher_exact(ct)
print(f"  carrier counts here: ins {ct[0]} vs snv {ct[1]} (NDD, population) -> OR {orr:.1f}, Fisher p={pf:.3g}")

# ---- 2. variant-level logistic: is the insertion effect explained by the assay?
print("\n=== 2. observed-in-NDD ~ SGE score (+ insertion) ===")
rows = []
for name, cols in [("SGE only", ["sge_score"]), ("insertion only", ["is_ins"]), ("SGE + insertion", ["sge_score", "is_ins"])]:
    X = sm.add_constant(D[cols].astype(float))
    m = sm.Logit(D.seen_ndd, X).fit(disp=0)
    r = dict(model=name, pseudo_r2=m.prsquared, llf=m.llf)
    for c in cols:
        r[f"beta_{c}"], r[f"p_{c}"] = m.params[c], m.pvalues[c]
    rows.append(r)
R = pd.DataFrame(rows)
print(R.round(4).to_string(index=False))
full = sm.Logit(D.seen_ndd, sm.add_constant(D[["sge_score", "is_ins"]].astype(float))).fit(disp=0)
print(f"  odds ratio for insertion, holding the assay score constant: {np.exp(full.params['is_ins']):.1f} "
      f"[{np.exp(full.conf_int().loc['is_ins', 0]):.1f}, {np.exp(full.conf_int().loc['is_ins', 1]):.1f}]")

# ---- 3. count model
print("\n=== 3. NDD carrier count ~ SGE score + insertion (negative binomial) ===")
try:
    nb = sm.GLM(D.ndd, sm.add_constant(D[["sge_score", "is_ins"]].astype(float)),
                family=sm.families.NegativeBinomial()).fit()
    print(pd.DataFrame({"coef": nb.params, "se": nb.bse, "p": nb.pvalues}).round(4).to_string())
    print(f"  rate ratio for insertion at equal assay score: {np.exp(nb.params['is_ins']):.1f}")
except Exception as e:
    print("  negative binomial failed:", e)

# ---- 4. the model-free statement
print("\n=== 4. matched on assay score: SNVs at least as damaging as n.64_65insT ===")
ins_t = D[D.hgvs == "n.64_65insT"].iloc[0]
worse = D[(D.is_ins == 0) & (D.sge_score <= ins_t.sge_score)]
print(f"  n.64_65insT SGE score {ins_t.sge_score:+.3f}; NDD carriers {int(ins_t.ndd)}")
print(f"  SNVs with an SGE score at least as damaging: {len(worse)}")
print(f"  of those, carriers in NDD: {int(worse.ndd.sum())} across {int((worse.ndd > 0).sum())} distinct variants")
print(f"  of those, carriers in population cohorts: {int(worse['pop'].sum())}")
print(worse[["hgvs", "sge_score", "ndd", "pop"]].sort_values("sge_score").to_string(index=False))

D[["key", "hgvs", "vtype", "n_pos", "sge_score", "ndd", "pop", "is_ins", "seen_ndd",
   "frac_states_protein_contact", "frac_states_snrna_contact", "max_protein_res"]].to_csv(
    config.RESULTS / "renu_insertion_gap.tsv", sep="\t", index=False)
R.to_csv(config.RESULTS / "renu_insertion_gap_models.tsv", sep="\t", index=False)
