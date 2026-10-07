"""Coordinate conversion between GRCh38 positions and non-coding transcript (HGVS n.) positions."""
import re

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def g_to_n(pos, start, end, strand):
    """1-based genomic position -> 1-based transcript position (may fall outside 1..len)."""
    return pos - start + 1 if strand == "+" else end - pos + 1


def n_to_g(n, start, end, strand):
    return start + n - 1 if strand == "+" else end - n + 1


_SNV = re.compile(r"^n\.(\d+)([ACGT])>([ACGT])$")
_INS = re.compile(r"^n\.(\d+)_(\d+)ins([ACGT]+)$")
_DEL = re.compile(r"^n\.(\d+)(?:_(\d+))?del$")


def left_normalize(chrom, pos, ref, alt, fetch):
    """Left-align a VCF-style indel. `fetch(chrom, pos)` returns the 1-based reference base."""
    while len(ref) != len(alt) and ref[-1] == alt[-1] and pos > 1:
        prev = fetch(chrom, pos - 1)
        ref = prev + ref[:-1]
        alt = prev + alt[:-1]
        pos -= 1
    return pos, ref, alt


def hgvs_n_to_vcf(hgvs, chrom, start, end, strand, fetch):
    """Convert an HGVS n. description (SNV, insertion, deletion) to a left-normalized VCF tuple."""
    m = _SNV.match(hgvs)
    if m:
        n, ref_t, alt_t = int(m.group(1)), m.group(2), m.group(3)
        pos = n_to_g(n, start, end, strand)
        ref, alt = (ref_t, alt_t) if strand == "+" else (revcomp(ref_t), revcomp(alt_t))
        if fetch(chrom, pos) != ref:
            raise ValueError(f"{hgvs}: reference mismatch at {chrom}:{pos}")
        return chrom, pos, ref, alt
    m = _INS.match(hgvs)
    if m:
        a, b, ins = int(m.group(1)), int(m.group(2)), m.group(3)
        if b != a + 1:
            raise ValueError(f"{hgvs}: insertion flanks must be adjacent")
        if strand == "+":
            anchor = n_to_g(a, start, end, strand)
            seq = ins
        else:
            anchor = n_to_g(b, start, end, strand)  # left genomic flank on minus strand
            seq = revcomp(ins)
        base = fetch(chrom, anchor)
        return (chrom, *left_normalize(chrom, anchor, base, base + seq, fetch))
    m = _DEL.match(hgvs)
    if m:
        a = int(m.group(1))
        b = int(m.group(2) or a)
        g1, g2 = sorted((n_to_g(a, start, end, strand), n_to_g(b, start, end, strand)))
        anchor = g1 - 1
        deleted = "".join(fetch(chrom, p) for p in range(g1, g2 + 1))
        base = fetch(chrom, anchor)
        return (chrom, *left_normalize(chrom, anchor, base + deleted, base, fetch))
    raise ValueError(f"unsupported HGVS: {hgvs}")
