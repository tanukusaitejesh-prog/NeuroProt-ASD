"""Minimal tabix (.tbi) reader with HTTP range fetching, so a few small loci can be pulled from a
multi-gigabyte remote bgzipped VCF without downloading it.

Used because pysam does not build on this machine and the HPRC pangenome VCFs are tens of GB.
Implements the TBI spec (samtools/htslib): BGZF container, linear + binning index, virtual offsets.
"""
import gzip
import io
import struct
import urllib.request
import zlib


def _get(url, start=None, end=None, timeout=300):
    req = urllib.request.Request(url)
    if start is not None:
        req.add_header("Range", f"bytes={start}-{'' if end is None else end}")
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def read_tbi(url, timeout=600):
    """Parse a .tbi into {chrom: {bin: [(voff_beg, voff_end), ...]}} plus the linear index."""
    raw = gzip.decompress(_get(url, timeout=timeout))
    assert raw[:4] == b"TBI\x01", "not a TBI index"
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", raw[4:36])
    names = raw[36:36 + l_nm].split(b"\x00")
    names = [n.decode() for n in names if n]
    off = 36 + l_nm
    idx = {}
    for ref in range(n_ref):
        (n_bin,) = struct.unpack("<i", raw[off:off + 4]); off += 4
        bins = {}
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack("<Ii", raw[off:off + 8]); off += 8
            chunks = []
            for _ in range(n_chunk):
                b, e = struct.unpack("<QQ", raw[off:off + 16]); off += 16
                chunks.append((b, e))
            bins[bin_id] = chunks
        (n_intv,) = struct.unpack("<i", raw[off:off + 4]); off += 4
        linear = list(struct.unpack(f"<{n_intv}Q", raw[off:off + 8 * n_intv])); off += 8 * n_intv
        idx[names[ref]] = (bins, linear)
    return idx


def reg2bins(beg, end):
    """Bins overlapping [beg, end) in the UCSC binning scheme used by tabix."""
    end -= 1
    out = [0]
    for shift, start in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        out.extend(range(start + (beg >> shift), start + (end >> shift) + 1))
    return out


def _bgzf_blocks(data):
    """Decompress a concatenation of BGZF blocks, tolerating a truncated final block."""
    out, p = io.BytesIO(), 0
    while p < len(data) - 18:
        if data[p:p + 2] != b"\x1f\x8b":
            break
        xlen = struct.unpack("<H", data[p + 10:p + 12])[0]
        bsize = None
        x, xend = p + 12, p + 12 + xlen
        while x < xend:
            si1, si2, slen = data[x], data[x + 1], struct.unpack("<H", data[x + 2:x + 4])[0]
            if si1 == 66 and si2 == 67:
                bsize = struct.unpack("<H", data[x + 4:x + 6])[0] + 1
            x += 4 + slen
        if bsize is None or p + bsize > len(data):
            break
        try:
            out.write(zlib.decompress(data[p + 12 + xlen:p + bsize - 8], -15))
        except zlib.error:
            break
        p += bsize
    return out.getvalue()


def fetch(vcf_url, idx, chrom, beg, end, max_bytes=80_000_000):
    """Yield VCF lines overlapping chrom:beg-end (1-based inclusive) from a remote bgzipped VCF."""
    if chrom not in idx:
        chrom = chrom.replace("chr", "") if chrom.startswith("chr") else "chr" + chrom
        if chrom not in idx:
            return []
    bins, linear = idx[chrom]
    b0, e0 = beg - 1, end
    want = []
    for b in reg2bins(b0, e0):
        want.extend(bins.get(b, []))
    if not want:
        return []
    # linear index gives a lower bound on the file offset for this window
    li = linear[min(b0 >> 14, len(linear) - 1)] if linear else 0
    want = [(s, e) for s, e in want if e > li]
    if not want:
        return []
    lo = min(s >> 16 for s, _ in want)
    hi = max(e >> 16 for _, e in want)
    if (hi - lo) > max_bytes:
        hi = lo + max_bytes
    blob = _get(vcf_url, lo, hi + 65536)
    text = _bgzf_blocks(blob).decode("utf8", "replace")
    out = []
    for line in text.split("\n"):
        if not line or line[0] == "#":
            continue
        f = line.split("\t", 2)
        if len(f) < 3 or f[0] != chrom:
            continue
        try:
            p = int(f[1])
        except ValueError:
            continue
        if beg <= p <= end:
            out.append(line)
    return out


def header(vcf_url, n=400_000):
    """Return the VCF header lines (and so the sample names)."""
    text = _bgzf_blocks(_get(vcf_url, 0, n)).decode("utf8", "replace")
    return [l for l in text.split("\n") if l.startswith("#")]
