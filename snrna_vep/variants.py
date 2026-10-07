"""Enumerate every possible small variant in a gene: all SNVs, 1-nt deletions and 1-nt insertions.

Each variant gets a transcript-level description (HGVS n.) and a left-normalized VCF key so that
patient variants, gnomAD records, SGE scores and external scores join on the same key.
"""
import pandas as pd

from .coords import left_normalize, n_to_g, revcomp

BASES = "ACGT"


def enumerate_variants(gene, seq_plus, fetch, indels=True):
    """`seq_plus` is the plus-strand reference of gene.start..gene.end."""
    L = gene.end - gene.start + 1
    rna = seq_plus if gene.strand == "+" else revcomp(seq_plus)
    rows = []
    for n in range(1, L + 1):
        ref_t = rna[n - 1]
        g = n_to_g(n, gene.start, gene.end, gene.strand)
        for alt_t in BASES:
            if alt_t == ref_t:
                continue
            ref, alt = (ref_t, alt_t) if gene.strand == "+" else (revcomp(ref_t), revcomp(alt_t))
            rows.append(dict(n_pos=n, hgvs=f"n.{n}{ref_t}>{alt_t}", vtype="snv",
                             chrom=gene.chrom, pos=g, ref=ref, alt=alt))
        if not indels:
            continue
        # 1-nt deletion of transcript position n
        anchor = g - 1
        base = fetch(gene.chrom, anchor)
        pos, ref, alt = left_normalize(gene.chrom, anchor, base + fetch(gene.chrom, g), base, fetch)
        rows.append(dict(n_pos=n, hgvs=f"n.{n}del", vtype="del", chrom=gene.chrom, pos=pos, ref=ref, alt=alt))
        # 1-nt insertion between n and n+1
        if n < L:
            for ins in BASES:
                anchor = g if gene.strand == "+" else n_to_g(n + 1, gene.start, gene.end, gene.strand)
                seq = ins if gene.strand == "+" else revcomp(ins)
                base = fetch(gene.chrom, anchor)
                pos, ref, alt = left_normalize(gene.chrom, anchor, base, base + seq, fetch)
                rows.append(dict(n_pos=n + 0.5, hgvs=f"n.{n}_{n + 1}ins{ins}", vtype="ins",
                                 chrom=gene.chrom, pos=pos, ref=ref, alt=alt))
    df = pd.DataFrame(rows)
    df.insert(0, "gene_name", gene.gene_name)
    df["key"] = df.chrom + ":" + df.pos.astype(str) + ":" + df.ref + ":" + df.alt
    # Equivalent indels (e.g. deleting either base of a homopolymer) share a key; keep one row.
    return df.drop_duplicates("key").reset_index(drop=True)


def apply_variant(rna, hgvs):
    """Return the mutant RNA (transcript orientation) for an n. SNV / 1-nt indel description."""
    import re
    m = re.match(r"^n\.(\d+)([ACGT])>([ACGT])$", hgvs)
    if m:
        i = int(m.group(1)) - 1
        return rna[:i] + m.group(3) + rna[i + 1:]
    m = re.match(r"^n\.(\d+)del$", hgvs)
    if m:
        i = int(m.group(1)) - 1
        return rna[:i] + rna[i + 1:]
    m = re.match(r"^n\.(\d+)_(\d+)ins([ACGT]+)$", hgvs)
    if m:
        i = int(m.group(1))
        return rna[:i] + m.group(3) + rna[i:]
    raise ValueError(hgvs)
