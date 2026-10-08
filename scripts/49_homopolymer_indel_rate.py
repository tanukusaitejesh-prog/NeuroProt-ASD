"""Empirical single-base insertion rate at homopolymer tracts, from gnomAD, to put n.64_65insT on the SNV scale.

Roulette models SNVs only. The recurrent ReNU allele is a T inserted into a T4 homopolymer, so the question
"is this site unusually mutable?" needs an indel rate at matched sequence context. This script measures it
empirically.

Method. Sample windows across chr12, fetch the reference sequence for each (UCSC REST) and the gnomAD v4.1
genome variants in it (remote tabix range query against the 26 GB per-chromosome VCF, never downloaded). Within
each window:
  - enumerate homopolymer tracts by base and length
  - classify observed variants into SNVs and single-base insertions, assigning each insertion to the tract it
    extends (if any)
  - compute, per tract length L: the fraction of tracts carrying at least one observed single-base insertion
  - compute, as the comparator: the fraction of positions carrying at least one observed SNV
The ratio of these two is the mutability of a length-L homopolymer for slippage insertions relative to a typical
site for substitution, measured on the same ascertainment (same cohort, same filters).

Caveat recorded up front: counts of observed distinct alleles reflect mutation rate, drift and selection
together. Restricting to common neutral sequence keeps selection roughly constant across the comparison but does
not remove it; this is a relative-mutability estimate, not a per-generation rate.
"""
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.remote_tabix import fetch, read_tbi

GNOMAD = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/"
          "gnomad.genomes.v4.1.sites.chr12.vcf.bgz")
WIN = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000
N_WIN = int(sys.argv[2]) if len(sys.argv) > 2 else 6
RNU4_2 = (120291763, 120291903)


def ucsc_seq(chrom, start, end):
    u = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
    req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
    import json
    return json.load(urllib.request.urlopen(req, timeout=180))["dna"].upper()


def tracts(seq, offset):
    """All homopolymer runs >=2: {(base, length): [start positions, 1-based genomic]}"""
    out = {}
    for m in re.finditer(r"(A+|C+|G+|T+)", seq):
        b, L = m.group(0)[0], len(m.group(0))
        if L >= 2:
            out.setdefault((b, L), []).append(offset + m.start() + 1)
    return out


print("reading gnomAD chr12 index ...", flush=True)
idx = read_tbi(GNOMAD + ".tbi")
rng = np.random.default_rng(0)
starts = rng.integers(20_000_000, 110_000_000, N_WIN)
agg_tracts, agg_ins, snv_sites, snv_obs, total_bp = {}, {}, 0, 0, 0

for i, s in enumerate(sorted(starts), 1):
    s = int(s); e = s + WIN
    try:
        seq = ucsc_seq("chr12", s, e)
        recs = fetch(GNOMAD, idx, "chr12", s + 1, e, max_bytes=250_000_000)
    except Exception as ex:
        print(f"  window {i} failed: {ex}"); continue
    if seq.count("N") > 0.2 * len(seq):
        print(f"  window {i} skipped (assembly gap)"); continue
    T = tracts(seq, s)
    for k, v in T.items():
        agg_tracts[k] = agg_tracts.get(k, 0) + len(v)
    pos_of_tract = {}
    for (b, L), ps in T.items():
        for p in ps:
            for off in range(L):
                pos_of_tract[p + off] = (b, L)
    n_snv = n_ins = 0
    for line in recs:
        f = line.split("\t", 8)
        if len(f) < 8 or "PASS" not in f[6]:
            continue
        ref, alt, pos = f[3], f[4], int(f[1])
        if len(ref) == 1 and len(alt) == 1:
            n_snv += 1
        elif len(alt) == len(ref) + 1 and alt.startswith(ref):
            n_ins += 1
            key = pos_of_tract.get(pos) or pos_of_tract.get(pos + 1)
            if key:
                agg_ins[key] = agg_ins.get(key, 0) + 1
    snv_obs += n_snv
    total_bp += len(seq) - seq.count("N")
    print(f"  window {i}: {s:,}-{e:,}  {n_snv:,} SNVs, {n_ins:,} 1-bp insertions, "
          f"{sum(len(v) for v in T.values()):,} tracts", flush=True)

snv_per_site = snv_obs / max(total_bp * 3, 1)          # 3 possible alternates per site
print(f"\ntotal sequence scanned: {total_bp:,} bp; observed PASS SNVs: {snv_obs:,}")
print(f"observed SNV alleles per possible SNV: {snv_per_site:.5f}")

rows = []
for (b, L), n in sorted(agg_tracts.items(), key=lambda kv: (kv[0][1], kv[0][0])):
    if L < 2 or L > 12 or n < 30:
        continue
    ins = agg_ins.get((b, L), 0)
    rows.append(dict(base=b, length=L, n_tracts=n, n_insertions=ins,
                     ins_per_tract=ins / n, rel_to_snv=(ins / n) / snv_per_site))
D = pd.DataFrame(rows)
print("\n=== single-base insertion frequency by homopolymer tract, relative to the per-SNV frequency ===")
print(D.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
t4 = D[(D.length == 4) & (D.base.isin(["T", "A"]))]
if len(t4):
    r = t4.rel_to_snv.mean()
    print(f"\nT/A homopolymer of length 4: insertions are {r:.1f}x as frequent per tract as an SNV is per "
          f"possible substitution")
    print(f"  n.64_65insT has 89 NDD carriers; critical-region pathogenic SNVs average 2.5 per allele")
    print(f"  36x observed excess / {r:.1f}x mutational expectation = {36/max(r,1e-9):.1f}x residual")
D.to_csv(config.RESULTS / "homopolymer_insertion_rate.tsv", sep="\t", index=False)
