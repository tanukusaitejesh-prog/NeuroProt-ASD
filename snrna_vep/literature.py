"""Literature HGVS n. -> GRCh38 VCF keys, using the RefSeq transcript numbering that papers use.

RefSeq transcripts were checked against GENCODE by sequence alignment (data/raw/refseq/*.fa): NR_023343.1 (RNU4ATAC)
starts 1 nt upstream of the GENCODE model, so n.(NR) = n.(GENCODE) + 1; RNU6ATAC, RNU12, RNU4-2, RNU6-1 and RNU5A-1
start at the same base. Validated on RNU4ATAC n.13C>T, n.16G>A, n.46G>A, n.50G>A, n.51G>A, n.55G>A, n.111G>A (all
reproduce their ClinVar P/LP records) and on RNU2-2 n.4G>A / n.35A>G (reproduce the curated dominant variants).
Every conversion is checked against the GRCh38 reference base; failures return None.
"""
import pandas as pd

from . import config
from .coords import hgvs_n_to_vcf
from .reference import fetch_base

NR = {"NR_023343": "RNU4ATAC", "NR_023344": "RNU6ATAC", "NR_029422": "RNU12", "NR_003137": "RNU4-2",
      "NR_004394": "RNU6-1", "NR_002756": "RNU5A-1"}
NR_SHIFT = {"RNU4ATAC": 1}
_genes = None


def to_key(gene, pos, change):
    """gene symbol, position string ('51' or '13_15'), change ('G>A', 'del', 'insT', 'dupX') -> 'chrom:pos:ref:alt'."""
    global _genes
    if _genes is None:
        _genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t").drop_duplicates("gene_name").set_index("gene_name")
    if gene not in _genes.index or "*" in pos or "-" in pos or "+" in pos:
        return None
    g = _genes.loc[gene]
    change = change.replace("U", "T").strip()
    sh = NR_SHIFT.get(gene, 0)
    pos = "_".join(str(int(x) - sh) for x in pos.split("_"))
    h = f"n.{pos}{change}"
    if change.startswith("dup"):
        a, b = pos.split("_")[0], pos.split("_")[-1]
        seq = change[3:]
        if not seq:                                   # spell out the duplicated bases from the reference
            from .coords import revcomp
            from .reference import fetch_region
            plus = fetch_region(g.chrom, int(g.start), int(g.end))
            rna = plus if g.strand == "+" else revcomp(plus)
            seq = rna[int(a) - 1:int(b)]
        h = f"n.{b}_{int(b) + 1}ins{seq}"
    if h is None:
        return None
    try:
        c, p, r, a = hgvs_n_to_vcf(h, g.chrom, int(g.start), int(g.end), g.strand, fetch_base)
        return f"{c}:{p}:{r}:{a}"
    except Exception:
        return None
