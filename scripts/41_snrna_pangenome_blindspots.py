"""Are the snRNA loci that population short-read sequencing cannot call actually variable?

gnomAD v4.1 (short read, 76,215 genomes) reports ZERO callable alleles at RNU1-3 and RNU1-4, and only 15% of
cohort depth at RNU2-1, while the near-identical paralogs RNU1-1 and RNU1-2 are at full depth (scripts/…,
results/snrna_callability.tsv). Those zeros are consistent with either (a) genuine invariance or
(b) an inability to map short reads there.

This script distinguishes them using the HPRC Release 2 (v2.1) Minigraph-Cactus pangenome, which is built from
long-read haplotype assemblies and therefore does not depend on short-read mappability. If a locus is variable in
the assemblies but invariant in gnomAD, the gnomAD result is an ascertainment artefact and the locus is a blind
spot for disease-variant discovery.

Remote tabix range queries (snrna_vep.remote_tabix) are used so the multi-GB VCF is never downloaded.
Output: results/snrna_pangenome_blindspots.tsv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.remote_tabix import fetch, header, read_tbi

BASE = ("https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/"
        "minigraph-cactus/v2.1/hprc-v2.1-mc-grch38/hprc-v2.1-mc-grch38")
VCF, TBI = BASE + ".wave.vcf.gz", BASE + ".wave.vcf.gz.tbi"
FLANK = 0

genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t")
call = pd.read_csv(config.RESULTS / "snrna_callability.tsv", sep="\t")
G = genes.merge(call[["gene_name", "pass_variants", "AN_median", "pct_cohort_callable"]], on="gene_name", how="left")

print("reading pangenome index ...", flush=True)
idx = read_tbi(TBI)
h = header(VCF)
samples = h[-1].split("\t")[9:] if h else []
n_hap = 2 * len(samples)
print(f"  HPRC v2.1: {len(samples)} samples (~{n_hap} haplotypes)\n", flush=True)

rows = []
for r in G.itertuples():
    try:
        recs = fetch(VCF, idx, r.chrom, int(r.start) - FLANK, int(r.end) + FLANK)
    except Exception as e:
        print(f"  {r.gene_name}: query failed ({e})")
        recs = None
    if recs is None:
        continue
    n_alt_hap = 0
    for line in recs:
        f = line.split("\t")
        gts = f[9:]
        for gt in gts:
            a = gt.split(":")[0].replace("|", "/").split("/")
            n_alt_hap += sum(1 for x in a if x not in (".", "0"))
    rows.append(dict(gene_name=r.gene_name, chrom=r.chrom, length=int(r.length),
                     paralog_group=r.paralog_group, is_disease_gene=bool(r.is_disease_gene),
                     gnomad_pass=int(r.pass_variants or 0), gnomad_pct_callable=r.pct_cohort_callable,
                     hprc_variants=len(recs), hprc_alt_haplotypes=n_alt_hap))
    print(f"  {r.gene_name:10s} gnomAD {int(r.pass_variants or 0):5d} ({r.pct_cohort_callable:5.1f}% callable)"
          f"   HPRC {len(recs):4d} variants, {n_alt_hap:5d} alt haplotypes", flush=True)

D = pd.DataFrame(rows)
D["hprc_var_per_kb"] = (1000 * D.hprc_variants / D.length).round(2)
D.to_csv(config.RESULTS / "snrna_pangenome_blindspots.tsv", sep="\t", index=False)

print("\n=== loci where short-read population data is blind ===")
blind = D[D.gnomad_pct_callable < 50]
print(blind[["gene_name", "paralog_group", "length", "gnomad_pass", "gnomad_pct_callable",
             "hprc_variants", "hprc_alt_haplotypes", "hprc_var_per_kb"]].to_string(index=False))
print("\n=== callable comparators in the same families ===")
comp = D[(D.gnomad_pct_callable >= 50) & (D.paralog_group.isin(blind.paralog_group.unique()))]
print(comp[["gene_name", "paralog_group", "length", "gnomad_pass", "gnomad_pct_callable",
            "hprc_variants", "hprc_alt_haplotypes", "hprc_var_per_kb"]].to_string(index=False))
if len(blind):
    print(f"\nVERDICT: blind loci carry {blind.hprc_variants.sum()} pangenome variants across "
          f"{int(blind.hprc_alt_haplotypes.sum())} alt haplotypes. "
          f"{'REAL VARIATION HIDDEN FROM SHORT READS' if blind.hprc_variants.sum() > 0 else 'genuinely invariant'}")
