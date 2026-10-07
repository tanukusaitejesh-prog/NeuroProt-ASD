"""UK Biobank WGS homozygote controls for recessive snRNA genes (Allele Frequency Browser exports, ~490k genomes).

Exports were made by a person in a web browser (afb.ukbiobank.ac.uk; region search, filter 'Number of Homozygotes' > 0,
CSV export) and saved as data/external/ukb_afb/<GENE>_homozygotes.csv. Automated querying is not permitted by the site.
Each variant is reference-checked and restricted to the transcribed region. Genotype-quality filter: allele number
>= 50% of the gene's median allele number (sites with many failed calls give unreliable homozygote counts). Output: data/curation/ukb_afb_homozygotes.tsv
with a 'conflict' flag for variants that are also labelled pathogenic (homozygous in >= 1 UKB adult).
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.reference import fetch_base

genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t").drop_duplicates("gene_name").set_index("gene_name")
L = pd.read_csv(config.CURATION / "labels_v4.tsv", sep="\t")
path = set(L[L.label == 1].key)
rows = []
for f in sorted((config.DATA / "external/ukb_afb").glob("*_homozygotes.csv")):
    gene = f.name.replace("_homozygotes.csv", "")
    g = genes.loc[gene]
    u = pd.read_csv(f)
    low_an = u["Allele Number"] < 0.5 * u["Allele Number"].median()
    for vid in u.loc[low_an, "Variant ID"]:
        print(f"  low allele number (QC), skipped: {vid}")
    u = u[~low_an]
    for vid, nhom, af in zip(u["Variant ID"], u["Number of Homozygotes"], u["Allele Frequency"]):
        c, p, r, a = vid.split("-")
        p = int(p)
        if fetch_base(c, p) != r[0]:
            print("  reference mismatch, skipped:", vid)
            continue
        if not (g.start - 1 <= p <= g.end) or int(nhom) < 1:
            continue
        rows.append(dict(gene_name=gene, key=f"{c}:{p}:{r}:{a}", ukb_nhom=int(nhom), ukb_af=float(af)))
out = pd.DataFrame(rows)
out["conflict"] = out.key.isin(path)
out.to_csv(config.CURATION / "ukb_afb_homozygotes.tsv", sep="\t", index=False)
print(out.groupby("gene_name").agg(variants=("key", "size"), conflicts=("conflict", "sum")).to_string())
