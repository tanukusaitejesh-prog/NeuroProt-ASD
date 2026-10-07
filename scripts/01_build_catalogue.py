"""Build the small-RNA gene catalogue from GENCODE v27 (GRCh38)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.genes import build_catalogue

config.RAW.mkdir(parents=True, exist_ok=True)
config.PROCESSED.mkdir(parents=True, exist_ok=True)
if not config.GENCODE_SMALL_RNA_GTF.exists():
    # Stream the public GENCODE v27 GTF over HTTPS and keep only small-RNA gene records.
    import re
    import requests
    url = config.GENCODE_GTF_GS.replace("gs://", "https://storage.googleapis.com/")
    keep = re.compile(r'gene_type "(snRNA|scaRNA|snoRNA|misc_RNA|ribozyme)"|gene_name "(RNU2-2P|CHASERR)"')
    with requests.get(url, stream=True, timeout=60) as r, open(config.GENCODE_SMALL_RNA_GTF, "w") as out:
        r.raise_for_status()
        for line in r.iter_lines(decode_unicode=True):
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene" and keep.search(f[8]):
                out.write(line + "\n")

cat = build_catalogue(config.GENCODE_SMALL_RNA_GTF)
out = config.PROCESSED / "small_rna_catalogue.tsv"
cat.to_csv(out, sep="\t", index=False)
print(f"{len(cat)} loci -> {out}")
print(cat.gene_type.value_counts().to_string())
print(cat[cat.paralog_group.notna() | cat.is_disease_gene][
    ["gene_name", "chrom", "start", "end", "strand", "length", "paralog_group", "disease_modes"]].to_string(index=False))
