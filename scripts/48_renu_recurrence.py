"""Why is one insertion 70-77% of all ReNU syndrome cases?

The 2026 Nature saturation-editing paper states the cause is "at present unknown", and rules out the obvious
explanation itself: n.64_65insT is NOT unusually damaging ("many variants scoring similarly or having even lower
function scores"). So the commonest cause of a prevalent neurodevelopmental disorder is not the most harmful
allele in the gene. The candidate explanations named in the literature are an elevated local mutation rate,
positive selection in the female germline, and ascertainment bias.

This script tests the mutation-rate explanation, which is the one that can be settled with public data.

Design. For every variant observed in NDD patients in the RNU4-2 critical region (counts transcribed from Chen
et al. 2024 Extended Data Table 1, data/curation/renu_variant_counts.tsv), compare the observed patient count
with the expected mutation rate at that site:
  SNVs    Roulette (Seplyarskiy et al.), base-pair-resolution germline mutation rate, queried remotely by CSI
          range request from the per-chromosome VCF (snrna_vep/remote_tabix.fetch_csi).
  indels  Roulette is SNV-only. n.64_65insT is a T inserted into a T4 homopolymer, the canonical replication
          slippage signature, so the insertion rate is estimated empirically from gnomAD: the density of
          observed single-base insertions at homopolymer tracts of matched length and base composition.

If observed/expected is comparable across variants, recurrence is mutational and no selection term is needed.
If n.64_65insT shows a large positive residual after rate correction, the residual is the quantity that
germline selection would have to explain.

This analysis was NOT pre-specified; it is a separate question from the variant-effect model and is reported
as its own study.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.remote_tabix import fetch_csi, read_csi

ROULETTE = ("https://genetics.bwh.harvard.edu/downloads/Vova/Roulette/"
            "12_rate_v5.2_TFBS_correction_all.vcf.bgz")
GENE = ("12", 120291763, 120291903)          # RNU4-2, GRCh38
CR = (120291825, 120291842)                   # the 18-nt critical region

# ---- 1. Roulette SNV rates across the gene
print("fetching Roulette rates for RNU4-2 ...", flush=True)
idx = read_csi(ROULETTE + ".csi")
recs = fetch_csi(ROULETTE, idx, GENE[0], GENE[1], GENE[2])
rows = []
for line in recs:
    f = line.split("\t")
    info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
    rows.append(dict(pos=int(f[1]), ref=f[3], alt=f[4],
                     MR=float(info.get("MR", "nan")), AR=float(info.get("AR", "nan"))))
RO = pd.DataFrame(rows)
RO["key"] = "chr12:" + RO.pos.astype(str) + ":" + RO.ref + ":" + RO.alt
print(f"  {len(RO)} SNV rate records; mean MR {RO.MR.mean():.4f}")
RO["in_CR"] = RO.pos.between(*CR)
print(f"  mean rate inside the critical region {RO[RO.in_CR].MR.mean():.4f} vs outside {RO[~RO.in_CR].MR.mean():.4f}")

# ---- 2. observed patient counts
C = pd.read_csv(config.CURATION / "renu_variant_counts.tsv", sep="\t", comment="#")
C["ndd"] = C.gel_ndd + C.nongel_ndd
C["pop"] = C.gnomad + C.ukb + C.aou
M = C.merge(RO[["key", "MR", "AR"]], on="key", how="left")
snv = M[M.vtype == "snv"].copy()
ins = M[M.vtype == "ins"].copy()
print(f"\nobserved variants: {len(snv)} SNVs, {len(ins)} insertions "
      f"({int(snv.ndd.sum())} and {int(ins.ndd.sum())} NDD carriers)")

# ---- 3. SNVs: does Roulette rate explain the observed spread?
s = snv[snv.MR.notna() & (snv.ndd + snv["pop"] > 0)].copy()
if len(s) > 4:
    from scipy.stats import spearmanr
    r1, p1 = spearmanr(s.MR, s.ndd)
    r2, p2 = spearmanr(s.MR, s["pop"])
    print(f"\n=== SNVs: Roulette rate vs observed counts ===")
    print(f"  rate vs NDD carriers:        rho {r1:+.2f}  p {p1:.3f}  (n={len(s)})")
    print(f"  rate vs population carriers: rho {r2:+.2f}  p {p2:.3f}")
    print("  (population counts are the cleaner test of the rate model: they are not filtered by disease)")

# ---- 4. the focal comparison, on the SNV scale
print("\n=== observed / expected for each observed critical-region SNV ===")
s["exp_rel"] = s.MR / s.MR.sum()
s["obs_rel"] = s.ndd / max(s.ndd.sum(), 1)
s["oe"] = s.obs_rel / s.exp_rel.replace(0, np.nan)
print(s.sort_values("ndd", ascending=False)[["hgvs", "ndd", "pop", "MR", "oe"]].round(3).to_string(index=False))

# ---- 5. the insertion: how far above the SNV scale is it?
INS_T = M[M.hgvs == "n.64_65insT"]
if len(INS_T):
    n_ins = int(INS_T.ndd.iloc[0])
    med_snv_rate = s.MR.median()
    print(f"\n=== the recurrent insertion ===")
    print(f"  n.64_65insT: {n_ins} NDD carriers, {int(INS_T['pop'].iloc[0])} population carriers")
    print(f"  total NDD carriers of critical-region SNVs: {int(s.ndd.sum())} across {len(s)} alleles")
    print(f"  so the single insertion outnumbers ALL observed SNVs in the region by "
          f"{n_ins / max(s.ndd.sum(), 1):.1f}-fold")
    print("\n  Roulette does not model indels. Quantifying the expected insertion rate at this T4 homopolymer")
    print("  against the SNV scale is the remaining step; see the companion script for the gnomAD-derived")
    print("  homopolymer indel rate.")
RO.to_csv(config.PROCESSED / "roulette_rnu4_2.tsv", sep="\t", index=False)
M.to_csv(config.RESULTS / "renu_recurrence_rates.tsv", sep="\t", index=False)
