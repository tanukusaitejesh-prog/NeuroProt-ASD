"""Sequence and RNA-structure features per variant (ViennaRNA), in transcript orientation."""
import numpy as np
import pandas as pd
import RNA

from .variants import apply_variant

# Intermolecular partners whose base-pairing a variant can break (snRNA duplexes).
DUPLEX_PARTNER = {
    "RNU4-1": "RNU6-1", "RNU4-2": "RNU6-1",
    "RNU6-1": "RNU4-2", "RNU6-2": "RNU4-2", "RNU6-7": "RNU4-2", "RNU6-8": "RNU4-2", "RNU6-9": "RNU4-2",
    "RNU4ATAC": "RNU6ATAC", "RNU6ATAC": "RNU4ATAC",
}


def to_rna(s):
    return s.replace("T", "U")


def mfe(seq):
    return RNA.fold_compound(to_rna(seq)).mfe()[1]


def cofold_dimer_energy(a, b):
    """Free energy of the A·B heterodimer (kcal/mol) via RNAcofold."""
    return RNA.fold_compound(to_rna(a) + "&" + to_rna(b)).mfe_dimer()[1]


def unpaired_probability(seq):
    fc = RNA.fold_compound(to_rna(seq))
    fc.mfe()
    fc.exp_params_rescale(fc.mfe()[1])
    fc.pf()
    bpp = np.array(fc.bpp())[1:, 1:]
    paired = bpp.sum(axis=0) + bpp.sum(axis=1)
    return 1 - paired


def gene_features(variants, rna, partner_rna=None):
    """`variants` from enumerate_variants for one gene; `rna` is its transcript sequence (DNA alphabet)."""
    L = len(rna)
    wt_mfe = mfe(rna)
    wt_dimer = cofold_dimer_energy(rna, partner_rna) if partner_rna else None
    p_unp = unpaired_probability(rna)
    out = []
    for v in variants.itertuples():
        mut = apply_variant(rna, v.hgvs)
        i = int(np.floor(v.n_pos)) - 1
        row = {
            "key": v.key,
            "rel_pos": v.n_pos / L,
            "dist_5p": v.n_pos - 1,
            "dist_3p": L - v.n_pos,
            "is_snv": v.vtype == "snv", "is_ins": v.vtype == "ins", "is_del": v.vtype == "del",
            "transition": v.vtype == "snv" and {v.hgvs[-3], v.hgvs[-1]} in ({"A", "G"}, {"C", "T"}),
            "local_gc": sum(c in "GC" for c in rna[max(0, i - 5): i + 6]) / len(rna[max(0, i - 5): i + 6]),
            "p_unpaired_wt": float(p_unp[min(i, L - 1)]),
            "ddG_fold": mfe(mut) - wt_mfe,
        }
        if partner_rna:
            row["ddG_duplex"] = cofold_dimer_energy(mut, partner_rna) - wt_dimer
        out.append(row)
    df = pd.DataFrame(out)
    if "ddG_duplex" not in df:
        df["ddG_duplex"] = np.nan
    return df


def genomic_edit_to_rna(key, gene, seq_plus):
    """Apply a VCF-keyed variant (any length) to a gene; return (mutant RNA, first changed n_pos, vtype)
    or None if the edit is not entirely inside the transcribed region."""
    from .coords import g_to_n, revcomp
    chrom, pos, ref, alt = key.split(":")
    pos = int(pos)
    # trim the shared VCF anchor/padding so only changed bases remain
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    lo, hi = pos, pos + max(len(ref), 1) - 1
    if not ref:                      # pure insertion between pos-1 and pos
        lo, hi = pos - 1, pos
    if lo < gene.start or hi > gene.end:
        return None
    i = pos - gene.start
    mut_plus = seq_plus[:i] + alt + seq_plus[i + len(ref):]
    mut = mut_plus if gene.strand == "+" else revcomp(mut_plus)
    first = g_to_n(pos if gene.strand == "+" else pos + max(len(ref), 1) - 1, gene.start, gene.end, gene.strand)
    vtype = "snv" if len(ref) == len(alt) == 1 else ("ins" if len(alt) > len(ref) else "del") \
        if (not ref or not alt) else "complex"
    return mut, max(1, first), vtype


def featurize_keys(keys, gene, seq_plus, partner_rna=None):
    """Same features as gene_features, for arbitrary variants given as VCF keys (multi-nt indels etc.)."""
    from .coords import revcomp
    rna = seq_plus if gene.strand == "+" else revcomp(seq_plus)
    L = len(rna)
    wt_mfe = mfe(rna)
    wt_dimer = cofold_dimer_energy(rna, partner_rna) if partner_rna else None
    p_unp = unpaired_probability(rna)
    rows, skipped = [], []
    for key in keys:
        r = genomic_edit_to_rna(key, gene, seq_plus)
        if r is None:
            skipped.append(key)
            continue
        mut, n, vtype = r
        i = n - 1
        rows.append({
            "key": key, "gene_name": gene.gene_name, "n_pos": n, "vtype": vtype,
            "rel_pos": n / L, "dist_5p": n - 1, "dist_3p": L - n,
            "is_snv": vtype == "snv", "is_ins": vtype == "ins", "is_del": vtype == "del",
            "transition": vtype == "snv" and {rna[i], mut[i]} in ({"A", "G"}, {"C", "T"}),
            "local_gc": sum(c in "GC" for c in rna[max(0, i - 5): i + 6]) / len(rna[max(0, i - 5): i + 6]),
            "p_unpaired_wt": float(p_unp[min(i, L - 1)]),
            "ddG_fold": mfe(mut) - wt_mfe,
            "ddG_duplex": (cofold_dimer_energy(mut, partner_rna) - wt_dimer) if partner_rna else np.nan,
        })
    return pd.DataFrame(rows), skipped
