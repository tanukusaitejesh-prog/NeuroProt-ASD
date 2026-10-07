"""Remote data locations. All are public Google Cloud Storage mirrors."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
CURATION = DATA / "curation"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

GCS = "https://storage.googleapis.com"
GRCH38_FASTA = f"{GCS}/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
GENCODE_GTF_GS = "gs://gcp-public-data--broad-references/hg38/v0/gencode.v27.primary_assembly.annotation.gtf"
GNOMAD_GENOMES_VCF = f"{GCS}/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.{{chrom}}.vcf.bgz"

GENCODE_SMALL_RNA_GTF = RAW / "gencode_v27_smallRNA_genes.gtf"
