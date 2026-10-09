"""First long-read variability map of the human snRNA gene family.

Motivation. The 2026 Nature Genetics systematic snRNA screen (PMC12424890) began from 200 potentially functional
snRNA genes, EXCLUDED 43 of them (10 with predominantly low variant-allele-fraction calls, 39 overlapping the GIAB
Problematic Regions track) and analysed 81, concluding that short-read sequencing "cannot accurately resolve highly
similar snRNA loci, including most U1 and U2 canonical genes ... The analysis of these regions thus requires
long-read sequencing technologies." That analysis has not been done.

This script does it: every snRNA gene in GENCODE for the spliceosomal families is queried against the HPRC
Release 2 (v2.1) Minigraph-Cactus pangenome, built from long-read haplotype assemblies and therefore free of
short-read mappability limits. For each locus we report how many variants and how many alternate haplotypes are
observed, which separates "no variation" from "no visibility".

Output: results/snrna_family_pangenome_map.tsv
"""
import gzip
import re
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.remote_tabix import fetch, header, read_tbi

GTF = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("gencode.gtf.gz")
BASE = ("https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/"
        "minigraph-cactus/v2.1/hprc-v2.1-mc-grch38/hprc-v2.1-mc-grch38")
VCF = BASE + ".wave.vcf.gz"
FAMILY = re.compile(r"^(RNU1-|RNVU1-|RNU2-|RNU4-|RNU5[A-Z]?-|RNU6-|RNU11|RNU12|RNU4ATAC|RNU6ATAC|RNU7-)", re.I)
MAIN = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}

print("parsing GENCODE for snRNA genes ...", flush=True)
rows = []
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.split("\t")
        if f[2] != "gene" or 'gene_type "snRNA"' not in f[8]:
            continue
        name = f[8].split('gene_name "')[1].split('"')[0]
        if f[0] not in MAIN or not FAMILY.match(name):
            continue
        rows.append(dict(gene=name, chrom=f[0], start=int(f[3]), end=int(f[4]), strand=f[6]))
G = pd.DataFrame(rows).drop_duplicates(["chrom", "start", "end"]).sort_values(["chrom", "start"])
G["family"] = G.gene.str.extract(r"^(RNVU1|RNU\d+[A-Z]*|RNU[0-9]+)", expand=False).str.upper()
print(f"  spliceosomal snRNA loci in GENCODE: {len(G)}")
print(G.family.value_counts().to_string())

print("\nreading pangenome index ...", flush=True)
idx = read_tbi(VCF + ".tbi")
n_hap = 2 * len(header(VCF)[-1].split("\t")[9:])
print(f"  HPRC v2.1: {n_hap} haplotypes\n", flush=True)

out, t0 = [], time.time()
for i, r in enumerate(G.itertuples(), 1):
    try:
        recs = fetch(VCF, idx, r.chrom, r.start, r.end)
    except Exception as e:
        out.append(dict(gene=r.gene, chrom=r.chrom, start=r.start, end=r.end, family=r.family,
                        length=r.end - r.start + 1, hprc_variants=None, hprc_alt_haps=None, note=str(e)[:60]))
        continue
    alt = 0
    for line in recs:
        for gt in line.split("\t")[9:]:
            a = gt.split(":")[0].replace("|", "/").split("/")
            alt += sum(1 for x in a if x not in (".", "0"))
    out.append(dict(gene=r.gene, chrom=r.chrom, start=r.start, end=r.end, family=r.family,
                    length=r.end - r.start + 1, hprc_variants=len(recs), hprc_alt_haps=alt, note=""))
    if i % 25 == 0:
        print(f"  {i}/{len(G)} loci ({time.time() - t0:.0f}s)", flush=True)

D = pd.DataFrame(out)
D["var_per_kb"] = (1000 * D.hprc_variants / D.length).round(1)
D["pct_alt_haps"] = (100 * D.hprc_alt_haps / n_hap).round(1)
D.to_csv(config.RESULTS / "snrna_family_pangenome_map.tsv", sep="\t", index=False)

print(f"\n=== long-read variability across {len(D)} spliceosomal snRNA loci ===")
print(D.groupby("family").agg(loci=("gene", "size"), median_var_per_kb=("var_per_kb", "median"),
                              invariant=("hprc_variants", lambda v: int((v == 0).sum())),
                              max_pct_alt=("pct_alt_haps", "max")).to_string())
print(f"\ngenuinely invariant in long-read assemblies (candidate constrained loci): "
      f"{int((D.hprc_variants == 0).sum())} of {len(D)}")
print("\nmost variable loci (least likely to be disease genes):")
print(D.nlargest(12, "var_per_kb")[["gene", "chrom", "length", "hprc_variants", "pct_alt_haps", "var_per_kb"]].to_string(index=False))
print("\ninvariant loci (the ones worth sequencing in patients):")
inv = D[D.hprc_variants == 0]
print(inv[["gene", "chrom", "start", "end", "family", "length"]].head(30).to_string(index=False))
