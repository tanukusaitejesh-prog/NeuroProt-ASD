"""GRCh38 reference access over HTTPS range requests (pysam + htslib remote FASTA)."""
from functools import lru_cache

import pysam

from .config import GRCH38_FASTA


@lru_cache(maxsize=1)
def _fasta():
    return pysam.FastaFile(GRCH38_FASTA)


@lru_cache(maxsize=4096)
def fetch_region(chrom, start, end):
    """1-based inclusive region, upper-case plus-strand sequence."""
    return _fasta().fetch(chrom, start - 1, end).upper()


def fetch_base(chrom, pos):
    # Cache whole 1-kb blocks to avoid one HTTP request per base.
    block = (pos - 1) // 1000
    seq = fetch_region(chrom, block * 1000 + 1, block * 1000 + 1000)
    return seq[(pos - 1) % 1000]
