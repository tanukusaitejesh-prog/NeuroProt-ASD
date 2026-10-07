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
