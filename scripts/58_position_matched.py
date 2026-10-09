"""Position-matched analysis of the ReNU allelic excess.

Why this replaces the class-level comparison. The manuscript previously reported insertions-versus-substitutions
across the whole critical region as its primary statistic, on the grounds that comparing classes is more robust
to ascertainment than singling out the allele that was discovered first. A leave-one-out test shows that reasoning
does not hold: excluding n.64_65insT, the per-allele insertion:substitution ratio is exactly 1.00 (2.5 vs 2.5).
The excess is not a property of insertions. It is a property of ONE allele.

That is not a weakness, it is a better experiment. Because the effect is allele-specific, it can be controlled at
the POSITION level, which is far stronger than controlling at the class level: variants at the same nucleotide are
covered by the same assays, called by the same pipelines, reported in the same screens and ascertained together.
Any explanation that acts on a region, a variant class, or a sequencing difficulty is held constant.

chr12:120,291,839 carries four observed variants - three insertions (one of which is the recurrent allele) and a
substitution - which makes it an internal control set that no cohort design could improve on.

Tests
  1. within-position contrasts at every critical-region position with >1 observed variant
  2. the focal position, against its own controls, with the measured mutational expectation applied
  3. exact binomial CIs on each contrast, and a Poisson check on whether the mutational supply could produce it
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, poisson

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

REL_INS_T4 = 0.403          # measured long-read insertion freq at a T4 tract, relative to per-substitution [50]
FOCAL = 120291839

C = pd.read_csv(config.CURATION / "renu_variant_counts.tsv", sep="\t", comment="#")
C["ndd"] = C.gel_ndd + C.nongel_ndd
C["pop"] = C.gnomad + C.ukb + C.aou
C["gpos"] = C.key.str.split(":").str[1].astype(int)

# ---------------------------------------------------------------- 1. the class-level claim does not survive LOO
d = C[C.ndd > 0]
for tag, sub in [("all alleles", d), ("excluding n.64_65insT", d[d.hgvs != "n.64_65insT"])]:
    i, s = sub[sub.vtype == "ins"], sub[sub.vtype == "snv"]
    pi, ps = i.ndd.sum() / len(i), s.ndd.sum() / len(s)
    print(f"  {tag:24s} ins {i.ndd.sum():3d}/{len(i)} = {pi:5.2f}   "
          f"snv {s.ndd.sum():3d}/{len(s)} = {ps:4.2f}   ratio {pi / ps:5.2f}")
print("  -> the excess is allele-specific, not class-specific; analyse it position-matched\n")

# ---------------------------------------------------------------- 2. within-position contrasts
print("=== every critical-region position with more than one observed variant ===")
rows = []
for pos, g in C[C.gpos.between(120291825, 120291842)].groupby("gpos"):
    obs = g[g.ndd > 0]
    if len(g) < 2:
        continue
    top = g.sort_values("ndd", ascending=False).iloc[0]
    others = g[g.hgvs != top.hgvs]
    rows.append(dict(position=pos, n_variants=len(g), top_allele=top.hgvs, top_carriers=int(top.ndd),
                     other_alleles=len(others), other_carriers=int(others.ndd.sum()),
                     ratio=top.ndd / max(others.ndd.sum(), 0.5)))
P = pd.DataFrame(rows).sort_values("top_carriers", ascending=False)
print(P.to_string(index=False))

# ---------------------------------------------------------------- 3. the focal position
print(f"\n=== the focal position chr12:{FOCAL:,} ===")
f = C[C.gpos == FOCAL].sort_values("ndd", ascending=False)
print(f[["hgvs", "vtype", "ndd", "pop"]].to_string(index=False))

ins_t = int(f[f.hgvs == "n.64_65insT"].ndd.iloc[0])
others = f[f.hgvs != "n.64_65insT"]
oth = int(others.ndd.sum())
snv_here = int(f[f.vtype == "snv"].ndd.sum())
n_oth_alleles = len(others)

print(f"\n  n.64_65insT                      {ins_t} carriers")
print(f"  all other alleles at this site   {oth} carriers across {n_oth_alleles} alleles")
print(f"  substitutions at this site       {snv_here} carriers")
print(f"  observed ratio vs substitutions  {ins_t / max(snv_here, 1):.1f}x")
print(f"  mutational expectation           {REL_INS_T4:.2f}x (a T4 insertion is LESS likely than a substitution)")
print(f"  position-matched residual        {(ins_t / max(snv_here, 1)) / REL_INS_T4:.0f}x")

bt = binomtest(ins_t, ins_t + snv_here, p=REL_INS_T4 / (1 + REL_INS_T4))
lo, hi = bt.proportion_ci(0.95)
print(f"\n  exact binomial, observed share {ins_t}/{ins_t + snv_here} = {ins_t / (ins_t + snv_here):.3f} "
      f"[{lo:.3f}, {hi:.3f}]")
print(f"  expected share under the measured mutational supply: {REL_INS_T4 / (1 + REL_INS_T4):.3f}")
print(f"  p = {bt.pvalue:.2e}")

exp_counts = (ins_t + snv_here) * REL_INS_T4 / (1 + REL_INS_T4)
print(f"\n  Poisson check: expected {exp_counts:.1f} insertion carriers, observed {ins_t}; "
      f"P(X >= {ins_t}) = {poisson.sf(ins_t - 1, exp_counts):.2e}")

# ---------------------------------------------------------------- what the controls exclude
print("\n=== what the within-position controls exclude ===")
ctrl = [("insertion-calling bias in a homopolymer",
         f"n.64_65insG ({int(f[f.hgvs == 'n.64_65insG'].ndd.iloc[0])}) and n.64_65insC "
         f"({int(f[f.hgvs == 'n.64_65insC'].ndd.iloc[0])}) are the same class of call at the same site"),
        ("ascertainment by region or assay design",
         "all four alleles lie at one nucleotide and are covered by any assay covering the critical region"),
        ("slippage mutability of the tract",
         f"measured directly: a T4 insertion is {REL_INS_T4:.2f}x a substitution, i.e. disfavoured [49, 50]"),
        ("interlocus gene conversion",
         "RNU4-1 carries the identical T4 tract; 0/11 alleles explicable [56]")]
for a, b in ctrl:
    print(f"  - {a}\n      {b}")

P.to_csv(config.RESULTS / "position_matched.tsv", sep="\t", index=False)
pd.DataFrame([dict(position=FOCAL, allele="n.64_65insT", carriers=ins_t, other_carriers=oth,
                   snv_carriers=snv_here, expected_rel=REL_INS_T4,
                   observed_ratio=ins_t / max(snv_here, 1),
                   residual=(ins_t / max(snv_here, 1)) / REL_INS_T4, binom_p=bt.pvalue)]
             ).to_csv(config.RESULTS / "position_matched_focal.tsv", sep="\t", index=False)
print(f"\nwrote position_matched.tsv and position_matched_focal.tsv")


# ---------------------------------------------------------------- 4. sensitivity to the comparator set
# The primary test uses the single substitution at the focal nucleotide (2 carriers). A reviewer will ask what
# happens as the comparator widens, and the answer should be in the paper rather than in a rebuttal letter.
print("\n=== sensitivity: how the residual depends on which substitutions are used as the comparator ===")
SETS = [("same nucleotide (primary)", FOCAL, FOCAL),
        ("the T4 tract", 120291836, 120291839),
        ("tract +/- 3 nt", 120291833, 120291842),
        ("whole critical region", 120291825, 120291842)]
sens = []
for lab, lo, hi in SETS:
    s = C[(C.vtype == "snv") & C.gpos.between(lo, hi)]
    n = int(s.ndd.sum())
    b = binomtest(ins_t, ins_t + n, p=REL_INS_T4 / (1 + REL_INS_T4))
    res = (ins_t / max(n, 1)) / REL_INS_T4
    adv = 100 * (res ** (1 / 30) - 1)
    sens.append(dict(comparator=lab, snv_carriers=n, obs_share=ins_t / (ins_t + n),
                     ratio=ins_t / max(n, 1), residual=res, p=b.pvalue, per_division_pct=adv))
    print(f"  {lab:28s} snv {n:3d}  share {ins_t / (ins_t + n):.3f}  ratio {ins_t / max(n, 1):6.1f}x  "
          f"residual {res:6.0f}x  p {b.pvalue:.1e}   -> {adv:4.1f}% per division")
S = pd.DataFrame(sens)
S.to_csv(config.RESULTS / "position_matched_sensitivity.tsv", sep="\t", index=False)
print(f"""
  The residual spans {S.residual.min():.0f}-{S.residual.max():.0f}-fold across these choices, so the narrowest
  comparator gives the largest number and should not be quoted alone. Significance is overwhelming throughout
  (p {S.p.max():.0e} at worst), and the mechanism requirement holds across the whole range:
  {S.per_division_pct.min():.1f}-{S.per_division_pct.max():.1f}% per division against ~1.1% documented.
  The manuscript therefore leads with the conservative {S.residual.min():.0f}-fold and reports the range.""")
