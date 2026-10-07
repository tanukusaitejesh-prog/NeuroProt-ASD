"""Interim patient labels from ClinVar (2025-05-04 release, mirrored in the public gnomAD bucket).

Pathogenic / likely pathogenic records in the core genes become a patient table in the standard schema,
with the disease mode taken from the gene. ClinVar has no per-variant first-report date here, so this
table supports leave-one-gene-out but not the time split. The curated supplement table supersedes it.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.coords import g_to_n
from snrna_vep.genes import DISEASE_GENES

URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/resources/grch38/clinvar/"
       "clinvar_20250504.vcf.gz")
PATHO = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
MODE = {g: m.split(";")[0] for g, m, *_ in DISEASE_GENES}

path = config.PROCESSED / "clinvar_20250504_core_loci.tsv"
if not path.exists():
    import pysam
    tb = pysam.TabixFile(URL)
    genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t")
    rows = []
    for g in genes.itertuples():
        for rec in tb.fetch(g.chrom.replace("chr", ""), g.start - 1, g.end):
            f = rec.split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            rows.append(dict(gene_name=g.gene_name, chrom="chr" + f[0], pos=int(f[1]), ref=f[3], alt=f[4],
                             clinvar_id=f[2], clnsig=info.get("CLNSIG"), revstat=info.get("CLNREVSTAT"),
                             disease=info.get("CLNDN", "")[:80], origin=info.get("ORIGIN"),
                             n_pos=g_to_n(int(f[1]), g.start, g.end, g.strand)))
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

cv = pd.read_csv(path, sep="\t")
cv = cv[cv.clnsig.isin(PATHO) & cv.gene_name.isin(MODE)]
out = pd.DataFrame({
    "gene_name": cv.gene_name, "hgvs_n": "", "mode": cv.gene_name.map(MODE), "inheritance": "unknown",
    "n_probands": 1, "first_published": "", "source": "ClinVar:" + cv.clinvar_id.astype(str),
    "cohort": "ClinVar-20250504", "notes": cv.clnsig + " | " + cv.disease.fillna(""),
    "key": cv.chrom + ":" + cv.pos.astype(str) + ":" + cv.ref + ":" + cv.alt,
})
dest = config.CURATION / "patient_variants_clinvar.tsv"
out.to_csv(dest, sep="\t", index=False)
print(out.groupby(["gene_name", "mode"]).size().to_string())
print(f"-> {dest}")
