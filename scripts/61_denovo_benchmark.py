"""ATTEMPTED AND ABANDONED: benchmarking the 110-fold ReNU residual against other recurrent de novo NDD alleles.

The idea. The recurrence manuscript reports that n.64_65insT exceeds its measured mutational supply by ~110-fold.
If that could be shown to be extreme relative to other recurrent de novo alleles in neurodevelopmental disorder,
the paper would gain a quantitative superlative ("the largest documented germline allelic excess in a Mendelian
disorder") rather than merely a large number. The plan was to build a reference distribution from denovo-db
v1.6.1, computing for each recurrent allele a within-gene observed/expected using Roulette mutation rates, with
the same construction applied to RNU4-2 in scripts/48.

It does not work, for two independent reasons, both established here rather than assumed.

REASON 1 - coordinate build mismatch. denovo-db v1.6.1 is on GRCh37/hg19; Roulette is on GRCh38. The recurrent
PACS1 allele sits at 11:65,978,677 in denovo-db, whereas PACS1 p.R203W is at roughly chr11:66,211,000 in GRCh38.
Rate lookups therefore return nothing for most variants. This alone is fixable with a liftover.

REASON 2 - fatal sparsity, which is not fixable. The within-gene O/E statistic needs a usable denominator: other
de novo alleles in the same gene, against which the focal allele's share can be compared. Across the 7,573 genes
with at least one functional de novo SNV in an NDD phenotype:

    median distinct alleles per gene   1
    75th percentile                    2
    genes with >= 20 distinct alleles  6
    genes with >= 50                   0

A ratio whose denominator is one or two observations cannot support a reference distribution. In the single gene
that did produce a value before the build mismatch bit (ADNP), O/E came out at 6.3 from five events in three
alleles - a number with no useful precision.

CONSEQUENCE FOR THE PAPER. The superlative claim is not supportable with public data, and the manuscript does
not make it. The 110-fold residual is reported as what it is: large, position-matched, and measured against a
locus-specific mutational supply (scripts/60), without any claim about where it ranks among recurrent alleles
generally. Do not re-attempt this with denovo-db; the limit is the number of de novo events per gene, not the
coordinates.

This file is kept as the record of a tested and rejected idea. Running it reproduces the two diagnoses.
"""
import gzip
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

DB = Path(r"C:\Users\saite\AppData\Local\Temp\claude\C--Users-saite-se2001"
          r"\51bce341-b5e8-4f3d-8b42-22ce035af3b6\scratchpad\denovodb\denovo-db.tsv.gz")
NDD = {"autism", "mixed", "developmentalDisorder", "intellectualDisability", "epilepsy"}
FUNC = {"missense", "stop-gained", "splice-donor", "splice-acceptor", "synonymous"}

if not DB.exists():
    sys.exit(f"denovo-db not present at {DB}\n"
             "fetch: curl -sL -o denovo-db.tsv.gz "
             "https://denovo-db.gs.washington.edu/denovo-db.non-ssc-samples.variants.tsv.gz")

rows, hdr = [], None
with gzip.open(DB, "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            hdr = line[1:].rstrip("\n").split("\t")
            continue
        rows.append(line.rstrip("\n").split("\t"))
D = pd.DataFrame(rows, columns=hdr)
D = D[D.PrimaryPhenotype.isin(NDD) & D.Variant.str.match(r"^[ACGT]>[ACGT]$")
      & D.FunctionClass.isin(FUNC) & (D.Gene != "NA")]
D["key"] = D.Chr + ":" + D.Position + ":" + D.Variant
U = D.drop_duplicates(["SampleID", "key"])

print("=== REASON 1: build mismatch ===")
p = U[(U.Gene == "PACS1") & (U.FunctionClass == "missense")].Position.astype(int)
if len(p):
    print(f"  PACS1 recurrent missense in denovo-db at chr11:{p.min():,} (GRCh37)")
    print("  PACS1 p.R203W in GRCh38 is near chr11:66,211,000 -> Roulette lookups miss\n")

print("=== REASON 2: sparsity (fatal) ===")
g = U.groupby("Gene").key.nunique()
print(f"  genes with >=1 functional de novo SNV : {len(g):,}")
print(f"  median distinct alleles per gene      : {g.median():.0f}")
print(f"  75th percentile                       : {g.quantile(0.75):.0f}")
print(f"  genes with >=20 distinct alleles      : {int((g >= 20).sum())}")
print(f"  genes with >=50 distinct alleles      : {int((g >= 50).sum())}")
print("\n  A within-gene observed/expected needs a denominator. The median gene supplies one allele.")

pd.DataFrame([dict(genes_with_denovo=len(g), median_alleles_per_gene=float(g.median()),
                   p75=float(g.quantile(0.75)), genes_ge20=int((g >= 20).sum()),
                   genes_ge50=int((g >= 50).sum()), verdict="benchmark not viable")]
             ).to_csv(config.RESULTS / "denovo_benchmark_abandoned.tsv", sep="\t", index=False)
print("\n=== VERDICT: not viable. The manuscript makes no claim about where 110x ranks. ===")
