"""Remeasure the mutational supply for the event class n.64_65insT ACTUALLY belongs to.

Why this script exists. The manuscript described n.64_65insT as "a thymine inserted into a four-thymine
homopolymer - the canonical signature of replication slippage". That is wrong, and the reference settles it:

    chr12  120,291,835  G
           120,291,836  T  |
           120,291,837  T  | T4 tract
           120,291,838  T  |
           120,291,839  T  |  <- reference base of the T>TA call
           120,291,840  C

The variant inserts an A between the final T of the tract and the C. The run stays T4. On the transcript it is a
T inserted 5' of the A4 run, and T != A, so it does not extend that either. It is a NON-EXTENDING insertion at a
tract edge, not a slippage extension.

scripts/49 and 50 counted any single-base insertion whose anchor fell in a tract, without checking the inserted
base, so the 0.40x figure mixes two event classes with very different mechanisms, and the log-linear
tract-length scaling that was reported as validation is a property of the extending class.

This script separates them, in the same windows with the same seed, and reports:
  EXTENDING      inserted base == the tract base       (true replication slippage)
  NON-EXTENDING  inserted base != the adjacent base    (the class n.64_65insT belongs to)
each relative to the substitution frequency measured the same way, so the ratio is directly comparable with the
figure the manuscript used.

Expectation, declared before running: non-extending insertions should be RARER than extensions, because slippage
is what makes long tracts insertion-prone. If so the supply for the real event class is below 0.40 and the
reported excess was understated, not overstated. Either way the manuscript must be rebuilt on this number.
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

GNOMAD = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/"
          "gnomad.genomes.v4.1.sites.chr12.vcf.bgz")
WIN = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000
N_WIN = int(sys.argv[2]) if len(sys.argv) > 2 else 10


def ucsc_seq(chrom, start, end):
    u = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
    req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
    return json.load(urllib.request.urlopen(req, timeout=180))["dna"].upper()


def tracts(seq, offset):
    """{(base, length): [1-based genomic start positions]} for every homopolymer run >= 2."""
    out = {}
    for m in re.finditer(r"(A+|C+|G+|T+)", seq):
        b, L = m.group(0)[0], len(m.group(0))
        if L >= 2:
            out.setdefault((b, L), []).append(offset + m.start() + 1)
    return out


print("reading gnomAD chr12 index ...", flush=True)
idx = read_tbi(GNOMAD + ".tbi")
rng = np.random.default_rng(0)                                 # SAME seed and windows as scripts/49 and 50
starts = sorted(int(s) for s in rng.integers(20_000_000, 110_000_000, N_WIN))

agg_tracts, ext, non_ext, snv_obs, total_bp = {}, {}, {}, 0, 0
for i, s in enumerate(starts, 1):
    e = s + WIN
    try:
        seq = ucsc_seq("chr12", s, e)
        recs = fetch(GNOMAD, idx, "chr12", s + 1, e, max_bytes=250_000_000)
    except Exception as ex:
        print(f"  window {i} failed: {ex}", flush=True)
        continue
    if seq.count("N") > 0.2 * len(seq):
        print(f"  window {i} skipped (assembly gap)", flush=True)
        continue
    T = tracts(seq, s)
    for k, v in T.items():
        agg_tracts[k] = agg_tracts.get(k, 0) + len(v)
    # map every position inside a tract to that tract's (base, length)
    pos_of = {}
    for (b, L), ps in T.items():
        for p in ps:
            for off in range(L):
                pos_of[p + off] = (b, L)
    n_snv = ne = nx = 0
    for line in recs:
        f = line.split("\t", 8)
        if len(f) < 8 or "PASS" not in f[6]:
            continue
        ref, alt, pos = f[3], f[4], int(f[1])
        if len(ref) == 1 and len(alt) == 1:
            n_snv += 1
        elif len(alt) == len(ref) + 1 and alt.startswith(ref):
            key = pos_of.get(pos) or pos_of.get(pos + 1)
            if not key:
                continue
            inserted = alt[len(ref):]                          # the single inserted base
            if inserted == key[0]:                             # matches the tract base -> extends the run
                ext[key] = ext.get(key, 0) + 1
                nx += 1
            else:                                              # does not extend: the n.64_65insT class
                non_ext[key] = non_ext.get(key, 0) + 1
                ne += 1
    snv_obs += n_snv
    total_bp += len(seq) - seq.count("N")
    print(f"  window {i}: {n_snv:,} SNVs, {nx:,} extending, {ne:,} non-extending insertions", flush=True)

snv_per_site = snv_obs / max(total_bp * 3, 1)
print(f"\nscanned {total_bp:,} bp; {snv_obs:,} PASS SNVs -> {snv_per_site:.6f} per possible substitution")

rows = []
for (b, L), n in sorted(agg_tracts.items(), key=lambda kv: (kv[0][1], kv[0][0])):
    if L < 2 or L > 12 or n < 30:
        continue
    x, y = ext.get((b, L), 0), non_ext.get((b, L), 0)
    rows.append(dict(base=b, length=L, n_tracts=n, n_extending=x, n_non_extending=y,
                     rel_extending=(x / n) / snv_per_site, rel_non_extending=(y / n) / snv_per_site,
                     rel_any=((x + y) / n) / snv_per_site))
D = pd.DataFrame(rows)
D.to_csv(config.RESULTS / "insertion_class_rate.tsv", sep="\t", index=False)

print("\n=== insertion frequency at homopolymer tracts, split by event class ===")
print("    (relative to the per-possible-substitution frequency, same windows)")
print(D.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

t4 = D[(D.length == 4) & D.base.isin(["A", "T"])]
if len(t4):
    rx, rn, ra = t4.rel_extending.mean(), t4.rel_non_extending.mean(), t4.rel_any.mean()
    print(f"\n=== A/T tracts of length 4 - the context of n.64_65insT ===")
    print(f"  extending (slippage)        {rx:.3f}x a substitution")
    print(f"  NON-extending (this allele) {rn:.3f}x a substitution   <-- the correct expectation")
    print(f"  pooled, as scripts/49 did   {ra:.3f}x")
    print(f"\n  the manuscript used {ra:.2f} (short read, pooled). The correct class is {rn:.3f}.")
    if rn > 0:
        print(f"  the expectation changes by {ra / rn:.2f}-fold, and every residual scales with it.")
    print("\n  NOTE: this is the short-read measurement. scripts/50 showed short reads under-detect "
          "homopolymer\n  insertions ~1.7x at length 4; that correction must be re-derived for this class "
          "before use.")
