"""gnomAD v4.1 genome variants for small-RNA loci, fetched by remote tabix."""
from functools import lru_cache

import pandas as pd
import pysam

from .config import GNOMAD_GENOMES_VCF
from .coords import g_to_n

INFO_FLOAT = ["AF", "cadd_phred", "phylop", "spliceai_ds_max", "pangolin_largest_ds", "AS_VQSLOD", "inbreeding_coeff"]
INFO_INT = ["AC", "AN", "nhomalt"]
INFO_FLAG = ["segdup", "lcr"]


@lru_cache(maxsize=32)
def _tabix(chrom):
    return pysam.TabixFile(GNOMAD_GENOMES_VCF.format(chrom=chrom))


def _parse_info(info):
    out = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            if k in INFO_FLOAT:
                out[k] = float(v) if v not in (".", "") else None
            elif k in INFO_INT:
                out[k] = int(v)
            elif k == "allele_type":
                out[k] = v
        elif item in INFO_FLAG:
            out[item] = True
    return out


def fetch_locus(gene, flank=0):
    """All PASS and non-PASS gnomAD genome records overlapping a gene (+/- flank)."""
    rows = []
    lo, hi = gene.start - flank, gene.end + flank
    for rec in _tabix(gene.chrom).fetch(gene.chrom, lo - 1, hi):
        f = rec.split("\t")
        pos, ref, alt, filt = int(f[1]), f[3], f[4], f[6]
        info = _parse_info(f[7])
        rows.append({
            "gene_name": gene.gene_name, "chrom": gene.chrom, "pos": pos, "ref": ref, "alt": alt,
            "filter": filt, **{k: info.get(k) for k in INFO_INT + INFO_FLOAT},
            "allele_type": info.get("allele_type"),
            "segdup": info.get("segdup", False), "lcr": info.get("lcr", False),
            "n_pos": g_to_n(pos, gene.start, gene.end, gene.strand),
        })
    return pd.DataFrame(rows)
