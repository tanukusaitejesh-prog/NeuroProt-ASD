"""Extract snRNA patient variants (HGVS n.) from Europe PMC full texts and supplements, for manual review.

For every HGVS n. mention: the gene is the nearest preceding snRNA gene symbol within 400 characters (else the
paper's main gene, given on the command line); a +-200 character context is kept, and words indicating benign/VUS
status are flagged. Each variant is converted to GRCh38 VCF with a reference-base check (a mismatch means the gene or
transcript assignment is wrong, and the row is rejected).
Usage: python scripts/05f_extract_literature.py <dir> PMCID:MAIN_GENE [...]   -> data/curation/literature_raw.tsv
Supplement text is read from <dir>/<PMCID>_supp/*.txt (pdftotext / xlsx->csv beforehand).
"""
import html
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.coords import hgvs_n_to_vcf
from snrna_vep.reference import fetch_base

GENES = ["RNU4ATAC", "RNU6ATAC", "RNU4-2", "RNU2-2P", "RNU2-2", "RNU5B-1", "RNU5A-1", "RNU12", "RNU11",
         "RNU6-1", "RNU6-2", "RNU6-8", "RNU6-9", "RNU6-7", "RNU4-1", "RNU1-1"]
ALIAS = {"RNU2-2": "RNU2-2P", "U4atac": "RNU4ATAC", "U6atac": "RNU6ATAC", "U12": "RNU12", "U11": "RNU11"}
GENE_RE = re.compile(r"\b(" + "|".join(map(re.escape, sorted(GENES + ["U4atac", "U6atac"], key=len, reverse=True))) + r")\b")
VAR_RE = re.compile(r"n\.\s?(\d+(?:_\d+)?)\s?([ACGTU]>[ACGTU]|del[ACGTU]*|ins[ACGTU]+|dup[ACGTU]*)")
# RefSeq accessions verified against NCBI (data/raw/refseq/*.fa). Literature numbers variants on these transcripts;
# NR_023343.1 (RNU4ATAC) starts 1 nt upstream of the GENCODE model, so n.(NR) = n.(GENCODE) + 1. The other transcripts
# start at the same base as GENCODE (checked by sequence alignment).
NR = {"NR_023343": "RNU4ATAC", "NR_023344": "RNU6ATAC", "NR_029422": "RNU12", "NR_003137": "RNU4-2",
      "NR_004394": "RNU6-1", "NR_002756": "RNU5A-1"}
NR_SHIFT = {"RNU4ATAC": 1}
NR_RE = re.compile(r"(NR_\d{6})")
FLAG = re.compile(r"benign|uncertain|VUS|polymorphism|control|gnomAD|population", re.I)
genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t").drop_duplicates("gene_name").set_index("gene_name")


def text_of(path):
    t = path.read_text(encoding="utf8", errors="ignore")
    if path.suffix == ".xml":
        t = re.sub(r"<[^>]+>", " ", t)
    return re.sub(r"\s+", " ", html.unescape(t)).replace(" ", " ")


def to_vcf(gene, pos, change):
    g = genes.loc[gene]
    change = change.replace("U", "T")
    sh = NR_SHIFT.get(gene, 0)
    pos = "_".join(str(int(x) - sh) for x in pos.split("_"))
    h = f"n.{pos}{change}"
    if change.startswith("dup"):
        a, b = (pos.split("_") + [pos])[:2]
        h = f"n.{b}_{int(b) + 1}ins{change[3:]}" if change[3:] else None
    if h is None:
        return None
    try:
        c, p, r, a = hgvs_n_to_vcf(h, g.chrom, int(g.start), int(g.end), g.strand, fetch_base)
        return f"{c}:{p}:{r}:{a}"
    except Exception:
        return None


root = Path(sys.argv[1])
rows = []
for spec in sys.argv[2:]:
    pmc, main = spec.split(":")
    files = [root / f"{pmc}.xml"] + sorted((root / f"{pmc}_supp").glob("*.txt"))
    for f in files:
        if not f.exists():
            continue
        t = text_of(f)
        for m in VAR_RE.finditer(t):
            before = t[max(0, m.start() - 400):m.start()]
            gs = GENE_RE.findall(before)
            gene = ALIAS.get(gs[-1], gs[-1]) if gs else main
            # a RefSeq accession next to the variant overrides the symbol (nearest within 150 characters)
            acc = [(abs(a.start() - m.start()), a.group(1)) for a in NR_RE.finditer(t, max(0, m.start() - 150), m.end() + 150)]
            acc = sorted(x for x in acc if x[1] in NR)
            gene_src = "accession" if acc else ("symbol" if gs else "main")
            if acc:
                gene = NR[acc[0][1]]
            ctx = t[max(0, m.start() - 200):m.end() + 200]
            key = to_vcf(gene, m.group(1), m.group(2))
            if key is None and gene != main and gene_src != "accession":   # retry with the paper's main gene
                key2 = to_vcf(main, m.group(1), m.group(2))
                if key2:
                    gene, key, gene_src = main, key2, "main-retry"
            rows.append(dict(pmcid=pmc, file=f.name, gene_name=gene, hgvs_n=f"n.{m.group(1)}{m.group(2)}", key=key, gene_src=gene_src,
                             flagged=bool(FLAG.search(t[max(0, m.start() - 120):m.end() + 120])), context=ctx))
out = pd.DataFrame(rows)
out.to_csv(config.CURATION / "literature_raw.tsv", sep="\t", index=False)
ok = out[out.key.notna()]
print(f"{len(out)} mentions, {out.key.isna().sum()} rejected (reference mismatch / unparsable)")
print(ok.drop_duplicates(["pmcid", "key"]).groupby(["pmcid", "gene_name"]).size().to_string())
