"""How many genes look constrained in gnomAD because they are UNCALLABLE rather than selected against?

Logic. Synonymous variants are close to neutral, so a gene's observed/expected synonymous count should sit near 1.
A large synonymous DEFICIT is therefore not evidence of selection - it is missing data: reads that could not be
mapped, or sites that failed QC. Crucially, LOEUF and pLI are also observed/expected ratios computed over the same
reads, so the same missing data pushes a gene toward looking MORE constrained, not less. A gene can therefore be
promoted into the "highly constrained, likely disease-relevant" tier purely by being hard to sequence.

This script quantifies that class genome-wide from the gnomAD v4.0 constraint metrics (MANE Select transcripts),
and asks how many of them carry a misleadingly low LOEUF.

Follow-up (scripts/41) validates a sample of these loci against the HPRC pangenome, where long-read haplotype
assemblies are not subject to short-read mappability.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("cm.tsv")
C = pd.read_csv(SRC, sep="\t", low_memory=False)
print(f"rows: {len(C)}")
C = C[C.mane_select.astype(str).str.lower().isin(["true", "1", "t"])].copy()
print(f"MANE Select transcripts: {len(C)}")

for c in ["syn.obs", "syn.exp", "syn.oe", "lof.obs", "lof.exp", "lof.oe", "lof.oe_ci.upper", "lof.pLI", "mis.oe"]:
    C[c] = pd.to_numeric(C[c], errors="coerce")
C = C[(C["syn.exp"] >= 5) & C["syn.oe"].notna() & C["lof.oe_ci.upper"].notna()].copy()
print(f"with usable synonymous expectation (exp >= 5): {len(C)}\n")

# Poisson lower tail: how unlikely is the observed synonymous count given expectation?
from scipy.stats import poisson
C["syn_p_depleted"] = poisson.cdf(C["syn.obs"], C["syn.exp"])
from statsmodels.stats.multitest import multipletests
C["syn_fdr"] = multipletests(C.syn_p_depleted.clip(1e-300, 1), method="fdr_bh")[1]

print("=== synonymous observed/expected across MANE genes ===")
print(C["syn.oe"].describe(percentiles=[.001, .01, .05, .25, .5, .75]).round(3).to_string())

LOEUF_CUT = 0.6          # the conventional "constrained" threshold
for cut, lab in [(0.5, "syn o/e < 0.50 (half the expected synonymous variation missing)"),
                 (0.7, "syn o/e < 0.70"),
                 (0.9, "syn o/e < 0.90")]:
    sub = C[C["syn.oe"] < cut]
    sig = sub[sub.syn_fdr < 0.05]
    con = sig[sig["lof.oe_ci.upper"] < LOEUF_CUT]
    print(f"\n{lab}: {len(sub)} genes; FDR<0.05 depleted: {len(sig)}; "
          f"of those, LOEUF < {LOEUF_CUT} (would be called constrained): {len(con)}")

SUS = C[(C["syn.oe"] < 0.5) & (C.syn_fdr < 0.05)].copy()
SUS["loeuf"] = SUS["lof.oe_ci.upper"]
print(f"\n=== the suspect class: {len(SUS)} genes with >=50% of synonymous variation missing (FDR<0.05) ===")
print(f"  median syn o/e {SUS['syn.oe'].median():.3f}; median LOEUF {SUS.loeuf.median():.3f} "
      f"(all MANE genes: {C['lof.oe_ci.upper'].median():.3f})")
print(f"  of these, {int((SUS.loeuf < 0.6).sum())} have LOEUF < 0.6 and {int((SUS['lof.pLI'] > 0.9).sum())} have pLI > 0.9")
print("\n  most extreme 25 by synonymous deficit:")
print(SUS.nsmallest(25, "syn.oe")[["gene", "syn.obs", "syn.exp", "syn.oe", "lof.obs", "lof.exp", "loeuf", "lof.pLI"]]
      .round(3).to_string(index=False))

# Is the suspect class enriched for apparent constraint relative to matched genes?
from scipy.stats import mannwhitneyu
rest = C[~C.gene.isin(SUS.gene)]
u, p = mannwhitneyu(SUS.loeuf.dropna(), rest["lof.oe_ci.upper"].dropna(), alternative="less")
print(f"\n  LOEUF of suspect genes vs the rest: Mann-Whitney one-sided p = {p:.3g} "
      f"({'suspect genes look MORE constrained' if p < 0.05 else 'no shift'})")

C[["gene", "transcript", "syn.obs", "syn.exp", "syn.oe", "syn_fdr", "lof.obs", "lof.exp",
   "lof.oe", "lof.oe_ci.upper", "lof.pLI", "mis.oe"]].to_csv(
    config.RESULTS / "constraint_callability_all.tsv", sep="\t", index=False)
SUS[["gene", "transcript", "syn.obs", "syn.exp", "syn.oe", "syn_fdr", "loeuf", "lof.pLI"]].to_csv(
    config.RESULTS / "constraint_callability_suspect.tsv", sep="\t", index=False)
print(f"\nwrote {len(SUS)} suspect genes to results/constraint_callability_suspect.tsv")
