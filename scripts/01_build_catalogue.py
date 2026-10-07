"""Build the small-RNA gene catalogue from GENCODE v27 (GRCh38)."""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.genes import build_catalogue

config.RAW.mkdir(parents=True, exist_ok=True)
config.PROCESSED.mkdir(parents=True, exist_ok=True)
if not config.GENCODE_SMALL_RNA_GTF.exists():
    cmd = (f"gsutil cat {config.GENCODE_GTF_GS} | awk -F'\\t' '$3==\"gene\"' | "
           "grep -E 'gene_type \"(snRNA|scaRNA|snoRNA|misc_RNA|ribozyme)\"|gene_name \"(RNU2-2P|CHASERR)\"' "
           f"> {config.GENCODE_SMALL_RNA_GTF}")
    subprocess.run(cmd, shell=True, check=True)

cat = build_catalogue(config.GENCODE_SMALL_RNA_GTF)
out = config.PROCESSED / "small_rna_catalogue.tsv"
cat.to_csv(out, sep="\t", index=False)
print(f"{len(cat)} loci -> {out}")
print(cat.gene_type.value_counts().to_string())
print(cat[cat.paralog_group.notna() | cat.is_disease_gene][
    ["gene_name", "chrom", "start", "end", "strand", "length", "paralog_group", "disease_modes"]].to_string(index=False))
