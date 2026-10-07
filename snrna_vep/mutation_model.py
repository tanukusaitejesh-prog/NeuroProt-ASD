"""Trinucleotide-context mutability model and observed/expected depletion.

The model is the fraction of possible SNVs (by context and alt allele, strand-collapsed) that are
observed in gnomAD genomes, fitted on flanking sequence around the snRNA loci (a near-neutral
background with the same local mappability). Expected count for a site = sum over its 3 SNVs.
"""
import numpy as np
import pandas as pd

from .coords import revcomp


def context_key(tri, alt):
    """Collapse to the pyrimidine-centred strand (C or T at the middle base)."""
    if tri[1] in "AG":
        tri, alt = revcomp(tri), revcomp(alt)
    return f"{tri}>{alt}"


def site_table(chrom, start, end, seq_plus, observed_snvs):
    """All possible SNVs at positions start+1..end-1 with an `observed` flag.

    `seq_plus` covers start..end (1-based inclusive); `observed_snvs` is a set of (pos, alt).
    """
    rows = []
    for i in range(1, len(seq_plus) - 1):
        tri = seq_plus[i - 1:i + 2]
        if "N" in tri:
            continue
        pos = start + i
        for alt in "ACGT":
            if alt == tri[1]:
                continue
            rows.append((chrom, pos, tri[1], alt, context_key(tri, alt), (pos, alt) in observed_snvs))
    return pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "ctx", "observed"])


def fit(sites, pseudocount=1.0):
    g = sites.groupby("ctx")["observed"].agg(["sum", "count"])
    p = (g["sum"] + pseudocount) / (g["count"] + 2 * pseudocount)
    return p.rename("p_obs")


def observed_expected(sites, model, window=1):
    """Per-position O/E (window=1) or centred sliding-window O/E over positions."""
    s = sites.assign(exp=sites.ctx.map(model).fillna(model.mean()))
    per = s.groupby("pos").agg(obs=("observed", "sum"), exp=("exp", "sum")).sort_index()
    if window > 1:
        per = per.rolling(window, center=True, min_periods=1).sum()
    per["oe"] = per.obs / per.exp
    return per


def poisson_depletion_p(obs, exp):
    from scipy.stats import poisson
    return poisson.cdf(obs, exp)


def windows_depleted(per, window, alpha=1e-3):
    w = per[["obs", "exp"]].rolling(window, center=True, min_periods=window).sum().dropna()
    w["oe"] = w.obs / w.exp
    w["p"] = poisson_depletion_p(w.obs.values, w.exp.values)
    w["significant"] = w.p < alpha
    return w


def chi2_fit_check(sites, model):
    """Calibration of the context model on its own training sites (sanity check)."""
    s = sites.assign(exp=sites.ctx.map(model))
    return float(s.observed.sum()), float(s.exp.sum()), float(np.corrcoef(
        s.groupby("ctx").observed.mean(), model.loc[s.groupby("ctx").observed.mean().index])[0, 1])
