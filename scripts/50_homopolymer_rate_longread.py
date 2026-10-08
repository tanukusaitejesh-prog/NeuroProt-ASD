"""Is the homopolymer insertion rate from short reads depressed by calling sensitivity? Re-measure with long reads.

scripts/49 estimated the single-base insertion frequency at homopolymer tracts from gnomAD v4.1 (short-read),
and found tract length 4 sits at 0.26x the substitution frequency. The entire quantitative claim of the
recurrence study rests on that number, and homopolymer indel calling is the best-known failure mode of
short-read sequencing. If gnomAD under-detects these insertions by 5-10x, the true rate could be ~2x rather than
0.26x and the reported residual would shrink several-fold.

This script repeats the identical measurement against the HPRC Release 2 pangenome, whose variants are derived
from long-read haplotype assemblies and therefore do not share that failure mode. Sample sizes differ enormously
(232 assemblies vs 76,215 genomes), so absolute frequencies are not comparable; the comparable quantity is the
insertion frequency RELATIVE to the substitution frequency measured in the same dataset and the same windows,
which is what both scripts report.

Expected outcomes:
  rel_to_snv(long read) ~ rel_to_snv(short read)   -> gnomAD was not materially biased; the residual stands
  rel_to_snv(long read) >> short read              -> gnomAD depressed the estimate; the residual shrinks by
                                                      the ratio, and the manuscript must be corrected
"""
import json
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.remote_tabix import fetch, read_tbi

HPRC = ("https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/"
        "minigraph-cactus/v2.1/hprc-v2.1-mc-grch38/hprc-v2.1-mc-grch38.wave.vcf.gz")
WIN = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000
N_WIN = int(sys.argv[2]) if len(sys.argv) > 2 else 10


def ucsc_seq(chrom, start, end):
    u = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
    req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
    return json.load(urllib.request.urlopen(req, timeout=180))["dna"].upper()


def tracts(seq, offset):
    out = {}
    for m in re.finditer(r"(A+|C+|G+|T+)", seq):
        b, L = m.group(0)[0], len(m.group(0))
        if L >= 2:
            out.setdefault((b, L), []).append(offset + m.start() + 1)
    return out


print("reading HPRC index ...", flush=True)
idx = read_tbi(HPRC + ".tbi")
rng = np.random.default_rng(0)                      # SAME seed and windows as scripts/49
starts = sorted(int(s) for s in rng.integers(20_000_000, 110_000_000, N_WIN))
agg_tracts, agg_ins, snv_obs, total_bp = {}, {}, 0, 0

for i, s in enumerate(starts, 1):
    e = s + WIN
    try:
        seq = ucsc_seq("chr12", s, e)
        recs = fetch(HPRC, idx, "chr12", s + 1, e, max_bytes=250_000_000)
    except Exception as ex:
        print(f"  window {i} failed: {ex}")
        continue
    if seq.count("N") > 0.2 * len(seq):
        print(f"  window {i} skipped (assembly gap)")
        continue
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
        f = line.split("\t", 6)
        if len(f) < 5:
            continue
        ref, alt, pos = f[3], f[4].split(",")[0], int(f[1])
        if len(ref) == 1 and len(alt) == 1 and alt in "ACGT":
            n_snv += 1
        elif len(alt) == len(ref) + 1 and alt.startswith(ref):
            n_ins += 1
            key = pos_of_tract.get(pos) or pos_of_tract.get(pos + 1)
            if key:
                agg_ins[key] = agg_ins.get(key, 0) + 1
    snv_obs += n_snv
    total_bp += len(seq) - seq.count("N")
    print(f"  window {i}: {n_snv:,} SNVs, {n_ins:,} 1-bp insertions", flush=True)

snv_per_site = snv_obs / max(total_bp * 3, 1)
print(f"\nHPRC: {total_bp:,} bp scanned, {snv_obs:,} SNVs -> {snv_per_site:.6f} per possible substitution")

rows = []
for (b, L), n in sorted(agg_tracts.items(), key=lambda kv: (kv[0][1], kv[0][0])):
    if L < 2 or L > 12 or n < 30:
        continue
    ins = agg_ins.get((b, L), 0)
    rows.append(dict(base=b, length=L, n_tracts=n, n_insertions=ins,
                     ins_per_tract=ins / n, rel_to_snv=(ins / n) / snv_per_site))
LR = pd.DataFrame(rows)
LR.to_csv(config.RESULTS / "homopolymer_insertion_rate_longread.tsv", sep="\t", index=False)
print("\n=== long-read (HPRC) relative insertion frequency ===")
print(LR.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

SR = pd.read_csv(config.RESULTS / "homopolymer_insertion_rate.tsv", sep="\t")
M = LR.merge(SR, on=["base", "length"], suffixes=("_lr", "_sr"))
M["ratio_lr_over_sr"] = M.rel_to_snv_lr / M.rel_to_snv_sr.replace(0, np.nan)
print("\n=== long read vs short read, same windows ===")
print(M[["base", "length", "rel_to_snv_sr", "rel_to_snv_lr", "ratio_lr_over_sr"]]
      .to_string(index=False, float_format=lambda x: f"{x:.3f}"))
t4 = M[(M.length == 4) & (M.base.isin(["A", "T"]))]
if len(t4):
    sr, lr = t4.rel_to_snv_sr.mean(), t4.rel_to_snv_lr.mean()
    print(f"\nT/A length 4: short read {sr:.3f}, long read {lr:.3f}  -> long/short = {lr/sr:.2f}x")
    print(f"  residual in the manuscript was 7.9 / {sr:.2f} = {7.9/sr:.0f}x")
    print(f"  recomputed with the long-read rate: 7.9 / {lr:.2f} = {7.9/lr:.0f}x")
M.to_csv(config.RESULTS / "homopolymer_rate_lr_vs_sr.tsv", sep="\t", index=False)
