"""Verify every computable number in MANUSCRIPT_RECURRENCE.md against its source data.

Written to be run before submission. Each check prints PASS or FAIL with the manuscript's claim and the value
recomputed from the committed data, so a mismatch is visible rather than inferred. Claims that are taken from the
literature (parental origin, cohort shares) are listed at the end as unverifiable here and must be checked
against the cited papers by hand.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, linregress, poisson, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

OK = []


def check(label, claimed, actual, tol=0.02, fmt="{:.3f}"):
    if isinstance(claimed, (int, np.integer)) and isinstance(actual, (int, np.integer)):
        good = claimed == actual
    else:
        good = abs(float(claimed) - float(actual)) <= tol * max(abs(float(claimed)), 1e-9) + 1e-12
    OK.append(good)
    tag = "PASS" if good else "**FAIL**"
    c = fmt.format(claimed) if not isinstance(claimed, (int, np.integer)) else str(claimed)
    a = fmt.format(actual) if not isinstance(actual, (int, np.integer)) else str(actual)
    print(f"  [{tag}] {label:<58s} manuscript {c:>12s}   recomputed {a:>12s}")


RO = pd.read_csv(config.PROCESSED / "roulette_rnu4_2.tsv", sep="\t")
C = pd.read_csv(config.CURATION / "renu_variant_counts.tsv", sep="\t", comment="#")
C["ndd"] = C.gel_ndd + C.nongel_ndd
C["pop"] = C.gnomad + C.ukb + C.aou
C["gpos"] = C.key.str.split(":").str[1].astype(int)
SR = pd.read_csv(config.RESULTS / "homopolymer_insertion_rate.tsv", sep="\t")
LR = pd.read_csv(config.RESULTS / "homopolymer_insertion_rate_longread.tsv", sep="\t")
M = pd.read_csv(config.RESULTS / "homopolymer_rate_lr_vs_sr.tsv", sep="\t")

print("\n=== 1. Roulette mutation rate ===")
check("423 Roulette records", 423, int(len(RO)))
check("mean rate inside critical region 0.063", 0.063, RO[RO.in_CR].MR.mean())
check("mean rate outside critical region 0.088", 0.088, RO[~RO.in_CR].MR.mean())

S = C.merge(RO.groupby("key").MR.mean().rename("rate"), on="key", how="left").dropna(subset=["rate"])
S = S[(S.ndd + S["pop"]) > 0]
r1, r2 = spearmanr(S.rate, S.ndd), spearmanr(S.rate, S["pop"])
check("n = 23 variants in the rate correlation", 23, int(len(S)))
check("rho rate vs patient carriers +0.18", 0.18, r1.statistic, tol=0.06)
check("p = 0.41", 0.41, r1.pvalue, tol=0.08)
check("rho rate vs population carriers -0.23", -0.23, r2.statistic, tol=0.08)
check("p = 0.29", 0.29, r2.pvalue, tol=0.10)
sep = ((S[S.ndd > 0]["pop"] == 0).all() and (S[S["pop"] > 0].ndd == 0).all())
OK.append(bool(sep))
print(f"  [{'PASS' if sep else '**FAIL**'}] patient-carrier and population-carrier variant sets are disjoint")

print("\n=== 2. Homopolymer scaling ===")
AT = SR[SR.base.isin(["A", "T"])]
lr_ = linregress(AT.length, np.log10(AT.rel_to_snv))
rs = spearmanr(AT.length, AT.rel_to_snv)
check("1.90-fold per added base", 1.90, 10 ** lr_.slope)
check("R2 = 0.90", 0.90, lr_.rvalue ** 2, tol=0.03)
check("Spearman rho = 0.99", 0.99, rs.statistic, tol=0.02)
check("p = 2.6e-18", 2.6e-18, rs.pvalue, tol=0.6, fmt="{:.2e}")
check("0.10 at tract length 2", 0.104, AT[AT.length == 2].rel_to_snv.mean(), tol=0.03)
check("24.5 at tract length 12", 24.5, AT[AT.length == 12].rel_to_snv.mean(), tol=0.12)

print("\n=== 3. The T4 tract ===")
t4 = SR[(SR.length == 4) & (SR.base == "T")].iloc[0]
a4 = SR[(SR.length == 4) & (SR.base == "A")].iloc[0]
check("short-read T4 relative frequency 0.257", 0.257, t4.rel_to_snv)
check("163 insertions at T4", 163, int(t4.n_insertions))
check("10,992 T4 tracts", 10992, int(t4.n_tracts))
check("short-read A4 relative frequency 0.208", 0.208, a4.rel_to_snv)
check("130 insertions at A4", 130, int(a4.n_insertions))
check("10,839 A4 tracts", 10839, int(a4.n_tracts))
check("21,831 A/T tracts at length 4 (Table 1)", 21831, int(t4.n_tracts + a4.n_tracts))

lt4 = LR[(LR.length == 4) & (LR.base == "T")].iloc[0]
la4 = LR[(LR.length == 4) & (LR.base == "A")].iloc[0]
check("long-read T4 0.543", 0.543, lt4.rel_to_snv)
check("long-read A4 0.262", 0.262, la4.rel_to_snv)
check("long-read mean at length 4 = 0.403", 0.403, (lt4.rel_to_snv + la4.rel_to_snv) / 2)
sr_mean = (t4.rel_to_snv + a4.rel_to_snv) / 2
check("long read is 1.73x the short-read estimate", 1.73, ((lt4.rel_to_snv + la4.rel_to_snv) / 2) / sr_mean)
check("short-read mean at length 4 = 0.233", 0.233, sr_mean)

MT = M[M.base.isin(["A", "T"])]
lo_rng = MT[MT.length.between(2, 4)].ratio_lr_over_sr
hi_rng = MT[MT.length.between(10, 12)].ratio_lr_over_sr
print(f"  [info] lr/sr ratio at lengths 2-4: {lo_rng.min():.1f}-{lo_rng.max():.1f} "
      f"(manuscript says 1.3-2.2)")
print(f"  [info] lr/sr ratio at lengths 10-12: {hi_rng.min():.1f}-{hi_rng.max():.1f} "
      f"(manuscript says 5.9-7.4)")

print("\n=== 4. Carrier counts and the class-level statistic ===")
d = C[C.ndd > 0]
ins, snv = d[d.vtype == "ins"], d[d.vtype == "snv"]
check("99 insertion carriers", 99, int(ins.ndd.sum()))
check("5 insertion alleles", 5, int(len(ins)))
check("15 substitution carriers", 15, int(snv.ndd.sum()))
check("6 substitution alleles", 6, int(len(snv)))
check("19.8 insertion carriers per allele", 19.8, ins.ndd.sum() / len(ins))
check("2.5 substitution carriers per allele", 2.5, snv.ndd.sum() / len(snv))
check("7.9-fold class-level excess", 7.9, (ins.ndd.sum() / len(ins)) / (snv.ndd.sum() / len(snv)))
loo = d[d.hgvs != "n.64_65insT"]
li, ls = loo[loo.vtype == "ins"], loo[loo.vtype == "snv"]
check("leave-one-out ratio = 1.00", 1.00, (li.ndd.sum() / len(li)) / (ls.ndd.sum() / len(ls)))

print("\n=== 5. The position-matched test (the primary analysis) ===")
f = C[C.gpos == 120291839]
check("4 variants observed at chr12:120,291,839", 4, int(len(f)))
check("89 carriers of n.64_65insT", 89, int(f[f.hgvs == "n.64_65insT"].ndd.iloc[0]))
check("2 carriers of n.64_65insG", 2, int(f[f.hgvs == "n.64_65insG"].ndd.iloc[0]))
check("0 carriers of n.64_65insC", 0, int(f[f.hgvs == "n.64_65insC"].ndd.iloc[0]))
check("2 carriers of n.65A>G", 2, int(f[f.vtype == "snv"].ndd.sum()))
check("93 patient carriers at the position", 93, int(f.ndd.sum()))

EXP = 0.403
nsnv = int(f[f.vtype == "snv"].ndd.sum())
bt = binomtest(89, 89 + nsnv, p=EXP / (1 + EXP))
lo, hi = bt.proportion_ci(0.95)
check("89 of 91 in the insertion-vs-substitution test", 91, int(89 + nsnv))
check("observed share 97.8%", 0.978, 89 / (89 + nsnv))
check("CI lower bound 92.3%", 0.923, lo, tol=0.01)
check("CI upper bound 99.7%", 0.997, hi, tol=0.01)
check("expected share 28.7%", 0.287, EXP / (1 + EXP))
check("exact binomial p = 1.3e-45", 1.3e-45, bt.pvalue, tol=0.25, fmt="{:.2e}")
expc = (89 + nsnv) * EXP / (1 + EXP)
check("Poisson expectation 26.1 carriers", 26.1, expc)
check("Poisson P(X>=89) = 5.2e-22", 5.2e-22, poisson.sf(88, expc), tol=0.3, fmt="{:.2e}")
check("position-matched residual 110x", 110, (89 / nsnv) / EXP, tol=0.02, fmt="{:.0f}")
check("44-fold vs the substitution at the same base", 44, 89 / nsnv, tol=0.03, fmt="{:.0f}")

print("\n=== 6. Gene conversion ===")
GC = pd.read_csv(config.RESULTS / "gene_conversion_test.tsv", sep="\t")
check("11 observed patient alleles tested", 11, int(len(GC)))
check("114 patient carriers covered", 114, int(GC.ndd_carriers.sum()))
check("0 alleles conversion-explicable", 0, int(GC.conversion_explicable.sum()))

COMP = str.maketrans("ACGTN", "TGCAN")
RF = pd.read_csv(config.PROCESSED / "reference_core_loci.tsv.gz", sep="\t").set_index("gene_name")
G = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t").set_index("gene_name")
mat = lambda g: RF.loc[g].seq.upper()[2000:len(RF.loc[g].seq) - 2000].translate(COMP)[::-1]
a2, a1 = mat("RNU4-2"), mat("RNU4-1")
diffs = [i + 1 for i, (x, y) in enumerate(zip(a2, a1)) if x != y]
check("137 of 141 identical", 137, int(sum(x == y for x, y in zip(a2, a1))))
check("97.2% identity", 0.972, sum(x == y for x, y in zip(a2, a1)) / 141, tol=0.005)
check("separation 1,194 bp", 1194, int(G.loc["RNU4-1"].start - G.loc["RNU4-2"].end))
print(f"  [info] difference positions: {diffs} (manuscript says n.37, n.88, n.99, n.113)")
in_cr = [p for p in diffs if 62 <= p <= 79]
OK.append(not in_cr)
print(f"  [{'PASS' if not in_cr else '**FAIL**'}] no paralogue difference inside the critical region")

print("\n=== 7. Cohort replication ===")
for arm, claim in [("gel_ndd", 11.2), ("nongel_ndd", 7.5)]:
    a = C[C[arm] > 0]
    i, s = a[a.vtype == "ins"], a[a.vtype == "snv"]
    check(f"{arm} per-allele ratio {claim}", claim, (i[arm].sum() / len(i)) / (s[arm].sum() / len(s)), tol=0.03)

print("\n=== 8. Paralogue control (RNU4-1) ===")
GN = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")
GN = GN[GN["filter"] == "PASS"].copy()
GN["isins"] = (GN.alt.str.len() == GN.ref.str.len() + 1)
GN["issnv"] = (GN.alt.str.len() == 1) & (GN.ref.str.len() == 1)
check("2,545 gnomAD PASS records at RNU4-1", 2545, int((GN.gene_name == "RNU4-1").sum()))
tr = GN[(GN.gene_name == "RNU4-1") & GN.pos.between(120293170, 120293173)]
check("1 insertion allele at the RNU4-1 tract", 1, int(tr.isins.sum()))
check("2 substitution alleles at the RNU4-1 tract", 2, int(tr.issnv.sum()))
check("insertion AC 5 at chr12:120,293,173", 5, int(tr[tr.isins].AC.iloc[0]))
check("tract ratio 0.50 insertions per substitution", 0.50, int(tr.isins.sum()) / int(tr.issnv.sum()))
for _g, _lo, _hi, _claim in [("RNU4-2", 120291763, 120291903, 0.091), ("RNU4-1", 120293097, 120293237, 0.064)]:
    b = GN[(GN.gene_name == _g) & GN.pos.between(_lo, _hi)]
    check(f"{_g} gene-body ratio {_claim}", _claim, int(b.isins.sum()) / int(b.issnv.sum()), tol=0.02)

print("\n=== 9. Selection requirement ===")
check("17% per division over 30 oogonial divisions", 17.0, 100 * (110 ** (1 / 30) - 1), tol=0.03, fmt="{:.1f}")
check("26% over 20 divisions", 26.5, 100 * (110 ** (1 / 20) - 1), tol=0.03, fmt="{:.1f}")
check("13% over 40 divisions", 12.5, 100 * (110 ** (1 / 40) - 1), tol=0.05, fmt="{:.1f}")
check("8% over 60 divisions", 8.1, 100 * (110 ** (1 / 60) - 1), tol=0.05, fmt="{:.1f}")
check("1.1% per division, selfish spermatogonial", 1.1, 100 * (1000 ** (1 / 610) - 1), tol=0.05, fmt="{:.1f}")
check("15-fold larger per division", 15, (110 ** (1 / 30) - 1) / (1000 ** (1 / 610) - 1), tol=0.05, fmt="{:.0f}")

print("\n" + "=" * 100)
print(f"{sum(OK)} of {len(OK)} automated checks PASS")
print("=" * 100)
print("""
NOT verifiable from these data - check by hand against the cited papers:
  - 77.4% share in the discovery cohort; 72.6% in the French cohort
  - parental origin 54/54 maternal (discovery); 47 maternal / 4 paternal of 51 (French)
  - the four paternal variants being n.62T>C, n.68A>C, n.76del, n.92C>G
  - saturation-editing claim that n.64_65insT is not unusually damaging
  - cohort sizes quoted in the header of renu_variant_counts.tsv
""")
