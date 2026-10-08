"""Paralog-resolved snRNA expression across human fetal brain, and what published resources cannot report.

Motivation. Variants in individual snRNA paralogs cause neurodevelopmental disease (RNU4-2, RNU2-2, RNU5B-1),
yet which paralog of each family is expressed where has never been measured systematically: near-identical copies
make read assignment hard, and the literature states that no study resolves paralog-level snRNA expression across
human brain regions.

scripts/43 established, from sequence alone, which paralogs are separable at all: U4 and U5 are fully assignable
at 101 nt, U2 is 92% assignable, and U1 and U6 are NOT assignable at any read length because their copies are
sequence-identical (U1: zero differences across 164 nt).

This script tests that prediction against real data and produces the measurement. ENCODE small RNA-seq
(101 nt, GRCh38) covers 7 fetal brain regions at 19-24 weeks plus neural stem cells, neural progenitors and
bipolar neurons. ENCODE's STAR gene quantifications count UNIQUELY mapping reads, so they are paralog-resolved
exactly where assignment is possible and structurally zero where it is not.

Two outputs:
  (a) the first paralog-resolved map of snRNA expression across fetal brain regions, for the families where it
      is possible (U4, U5);
  (b) a quantification of what is lost: the fraction of the small RNA-seq library that is discarded as
      multimapping, and the genes for which a published resource reports zero despite being among the most
      abundant RNAs in the cell.
"""
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

MANIFEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("quant_files.json")
files = json.load(open(MANIFEST))
genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t")
ids = {g.gene_id.split(".")[0]: g.gene_name for g in genes.itertuples()}
fam = {g.gene_name: g.paralog_group for g in genes.itertuples()}

rows, libs = [], []
for i, f in enumerate(files, 1):
    try:
        raw = urllib.request.urlopen(f["url"], timeout=180).read().decode("utf8", "replace")
    except Exception as e:
        print(f"  {f['file']}: {e}")
        continue
    counts, meta = {}, {}
    for line in raw.split("\n"):
        p = line.split("\t")
        if len(p) < 2:
            continue
        if p[0].startswith("N_"):
            meta[p[0]] = int(p[1])
        else:
            b = p[0].split(".")[0]
            if b in ids:
                counts[ids[b]] = int(p[1])
    assigned = sum(1 for _ in [0])  # placeholder, replaced below
    total_unique = sum(int(p[1]) for line in raw.split("\n") if (p := line.split("\t")) and len(p) > 1
                       and not p[0].startswith("N_") and p[1].isdigit())
    libs.append(dict(file=f["file"], exp=f["exp"], sample=f["sample"], rep=f["rep"],
                     unique_assigned=total_unique, multimapping=meta.get("N_multimapping", 0),
                     unmapped=meta.get("N_unmapped", 0), ambiguous=meta.get("N_ambiguous", 0)))
    for g, c in counts.items():
        rows.append(dict(file=f["file"], exp=f["exp"], sample=f["sample"], rep=f["rep"],
                         gene=g, family=fam.get(g), count=c))
    if i % 5 == 0:
        print(f"  {i}/{len(files)} files", flush=True)

E = pd.DataFrame(rows)
L = pd.DataFrame(libs)
E.to_csv(config.RESULTS / "snrna_paralog_expression_raw.tsv", sep="\t", index=False)
L.to_csv(config.RESULTS / "snrna_library_composition.tsv", sep="\t", index=False)
print(f"\nsamples: {L.file.nunique()}; snRNA gene records: {len(E)}")

# ---- (b) how much of the library is discarded
L["pct_multimapping"] = 100 * L.multimapping / (L.unique_assigned + L.multimapping + L.unmapped + L.ambiguous)
print("\n=== fraction of each small RNA-seq library discarded as multimapping ===")
print(L[["sample", "rep", "multimapping", "pct_multimapping"]].assign(
    sample=L["sample"].str.slice(0, 44)).to_string(index=False))
print(f"\nmedian across samples: {L.pct_multimapping.median():.1f}%")

# ---- which genes are reported as zero despite being abundant
Z = E.groupby(["gene", "family"]).agg(samples=("count", "size"), zero_in=("count", lambda v: int((v == 0).sum()),),
                                      median_count=("count", "median"), max_count=("count", "max")).reset_index()
print("\n=== genes reported as ZERO in every sample ===")
print(Z[Z.max_count == 0][["gene", "family", "samples"]].to_string(index=False))
print("\n=== all snRNA genes, median count across samples ===")
print(Z.sort_values("median_count", ascending=False)[["gene", "family", "median_count", "max_count", "zero_in", "samples"]].to_string(index=False))

# ---- (a) paralog ratios where assignment is possible
print("\n=== paralog-resolved expression, U4 and U5 (assignable families) ===")
piv = E[E.family.isin(["U4", "U5"])].pivot_table(index=["sample", "rep"], columns="gene", values="count", aggfunc="first")
piv = piv.loc[:, piv.max() > 0]
print(piv.astype("Int64").to_string())
if {"RNU4-2", "RNU4-1"}.issubset(piv.columns):
    r = (piv["RNU4-2"] / piv["RNU4-1"].replace(0, np.nan)).dropna()
    print(f"\nRNU4-2 : RNU4-1 ratio across samples - median {r.median():.2f}, range {r.min():.2f}-{r.max():.2f}")
    print(r.round(2).to_string())
if {"RNU5A-1", "RNU5B-1"}.issubset(piv.columns):
    r5 = (piv["RNU5A-1"] / piv["RNU5B-1"].replace(0, np.nan)).dropna()
    print(f"\nRNU5A-1 : RNU5B-1 ratio - median {r5.median():.2f}, range {r5.min():.2f}-{r5.max():.2f}")
Z.to_csv(config.RESULTS / "snrna_paralog_expression_summary.tsv", sep="\t", index=False)
piv.to_csv(config.RESULTS / "snrna_paralog_matrix.tsv", sep="\t")
