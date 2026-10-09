"""Figure 2 for the MuST-VEP manuscript: the dissociation between clinical labels and molecular function.

This is the paper's lead result and it needs to be visible in one figure. A pre-specified gate failed because
gnomAD population constraint was not separable from MuST-VEP on curated clinical labels. On a saturation genome
editing readout - which never consults allele frequency - the two separate decisively, and the partial correlations
are asymmetric. The explanation is that curated labels are themselves partly derived from population depletion, so
a depletion score predicts them partly for free.

  a  the dissociation: the same two scores on clinical labels (AUROC, overlapping) and on measured function
     (Spearman, separated). One panel, the whole argument.
  b  zero-shot Spearman with measured function for every score, with position-block 95% CIs
  c  partial correlations: structure given constraint, and constraint given structure
  d  scatter of MuST-VEP against the assay, with CADD on the same 485 variants for contrast
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, GREEN, AMBER, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#b9b7b2"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": SURF})
LAB = config.LABELS.replace(".tsv", "")
R = config.RESULTS


def partial_spearman(x, y, z):
    """Spearman of x with y after removing the rank-linear effect of z from both.

    Identical definition to scripts/27 (Pearson of the rank residuals), so the figure and the text report the
    same number.
    """
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    rx, ry, rz = rankdata(x[ok]), rankdata(y[ok]), rankdata(z[ok])
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    return float(np.corrcoef(ex, ey)[0, 1])


RHO = pd.read_csv(R / f"constraint_vs_function_rho_{LAB}.tsv", sep="\t").set_index("score")
V = pd.read_csv(R / f"constraint_vs_function_variants_{LAB}.tsv", sep="\t")
G2 = pd.read_csv(R / f"gate_g2_tests_{LAB}.tsv", sep="\t")
g2d = G2[G2.ctx == "dominant"].set_index("baseline")
y = -V.sge_score.to_numpy(float)

fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.0, 1.25],
                      hspace=0.52, wspace=0.34, left=0.115, right=0.955, top=0.92, bottom=0.085)

# ------------------------------------------------- a. the dissociation, as paired differences
# AUROC and Spearman are different quantities and must not share an axis. What IS comparable is the paired
# difference (MuST-VEP - constraint) within each readout, with its own interval: on clinical labels it covers
# zero, on measured function it does not. That contrast is the finding.
ax = fig.add_subplot(gs[0, 0])
PAIR = pd.read_csv(R / f"constraint_vs_function_paired_{LAB}.tsv", sep="\t")
pf = PAIR[(PAIR.structure == "snrnavep_v2") & (PAIR.population == "constraint_1_minus_oe")].iloc[0]
gl = g2d.loc["constraint_1_minus_oe"]
rows = [("measured\nfunction", "$\\Delta\\rho$, 485 assayed", pf.delta_rho, pf.lo, pf.hi, BLUE),
        ("clinical\nlabels", "$\\Delta$AUROC, held-out gene", gl.delta, gl.block_lo, gl.block_hi, ORANGE)]
for i, (lab, sub, d, lo, hi, c) in enumerate(rows):
    ax.errorbar(d, i, xerr=[[d - lo], [hi - d]], fmt="o", ms=7, color=c, ecolor=c, lw=1.8, capsize=3.5, zorder=3)
    ax.text(d, i + 0.30, f"{d:+.3f}  [{lo:+.3f}, {hi:+.3f}]", ha="center", va="bottom", fontsize=6.5, color=INK)
    ax.text(d, i - 0.32, sub, ha="center", va="top", fontsize=6.2, color=INK2)
ax.axvline(0, color=INK, lw=1.0, zorder=1)
ax.set_yticks([0, 1]); ax.set_yticklabels([r[0] for r in rows], fontsize=7.4)
ax.set_ylim(-0.75, 1.75); ax.set_xlim(-0.26, 0.48)
ax.set_xlabel("MuST-VEP $-$ gnomAD constraint", fontsize=7.8)
ax.text(-0.235, 1.52, "not separable", fontsize=6.8, color=ORANGE, weight="bold")
ax.text(-0.235, 0.52, "separable", fontsize=6.8, color=BLUE, weight="bold")
ax.set_title("The same comparison, two readouts", fontsize=8.5, color=INK, loc="left", pad=5)
ax.text(-0.30, 1.17, "a", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

# ------------------------------------------------- b. all scores vs function, with CIs
ax = fig.add_subplot(gs[0, 1])
NAME = {"nb_frac_prot": "Graph-smoothed contacts", "snrnavep_v2": "MuST-VEP",
        "frac_states_protein_contact": "Protein-contact fraction", "pop_density_w2": "gnomAD density (±2 nt)",
        "constraint_1_minus_oe": "gnomAD constraint", "pop_density_pos": "gnomAD density (position)",
        "cadd_phred": "CADD", "phylop447": "phyloP447"}
S = RHO.loc[[k for k in NAME if k in RHO.index]].sort_values("rho")
cols = [BLUE if k in ("snrnavep_v2", "nb_frac_prot", "frac_states_protein_contact") else ORANGE
        if k.startswith(("constraint", "pop_density")) else GREY for k in S.index]
ax.barh(range(len(S)), S.rho, color=cols, edgecolor=INK2, lw=0.5, height=0.68, zorder=3)
ax.errorbar(S.rho, range(len(S)), xerr=[S.rho - S.lo, S.hi - S.rho], fmt="none", ecolor=INK2, lw=0.8,
            capsize=2.2, zorder=4)
ax.set_yticks(range(len(S))); ax.set_yticklabels([NAME[k] for k in S.index], fontsize=7.2)
ax.axvline(0, color=INK2, lw=0.7)
ax.set_xlabel("Spearman $\\rho$ with measured damage (485 variants, zero-shot)", fontsize=7.6)
ax.set_xlim(-0.2, 0.68)
ax.set_title("Structure tracks function", fontsize=8.5, color=INK, loc="left", pad=5)
ax.text(-0.30, 1.17, "b", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

# ------------------------------------------------- c. partial correlations
ax = fig.add_subplot(gs[1, 0])
p_sc = partial_spearman(V.snrnavep_v2.to_numpy(float), y, V.constraint_1_minus_oe.to_numpy(float))
p_cs = partial_spearman(V.constraint_1_minus_oe.to_numpy(float), y, V.snrnavep_v2.to_numpy(float))
p_gc = partial_spearman(V.nb_frac_prot.to_numpy(float), y, V.constraint_1_minus_oe.to_numpy(float))
p_cg = partial_spearman(V.constraint_1_minus_oe.to_numpy(float), y, V.nb_frac_prot.to_numpy(float))
bars = [("MuST-VEP\ngiven constraint", p_sc, BLUE), ("constraint\ngiven MuST-VEP", p_cs, ORANGE),
        ("contacts\ngiven constraint", p_gc, BLUE), ("constraint\ngiven contacts", p_cg, ORANGE)]
ax.barh(range(4), [b[1] for b in bars], color=[b[2] for b in bars], edgecolor=INK2, lw=0.5, height=0.66, zorder=3)
for i, (_, v, _) in enumerate(bars):
    ax.text(v + 0.012, i, f"{v:+.3f}", va="center", fontsize=7.2, color=INK)
ax.set_yticks(range(4)); ax.set_yticklabels([b[0] for b in bars], fontsize=6.9)
ax.axvline(0, color=INK2, lw=0.7); ax.set_xlim(-0.03, 0.52)
ax.set_xlabel("partial Spearman $\\rho$ with measured damage", fontsize=7.6)
ax.set_title("The information is not shared", fontsize=8.5, color=INK, loc="left", pad=5)
ax.text(-0.30, 1.17, "c", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

# ------------------------------------------------- d. scatter
ax = fig.add_subplot(gs[1, 1])
for col, c, lab, rho in [("snrnavep_v2", BLUE, "MuST-VEP", RHO.loc["snrnavep_v2", "rho"]),
                         ("cadd_phred", GREY, "CADD", RHO.loc["cadd_phred", "rho"])]:
    m = V[col].notna() & np.isfinite(y)
    xr = rankdata(V.loc[m, col]) / m.sum()
    ax.scatter(xr, y[m.to_numpy()], s=5, alpha=0.4, color=c, lw=0, zorder=3 if col == "snrnavep_v2" else 2,
               label=f"{lab}  $\\rho$ = {rho:.3f}")
    z = np.polyfit(xr, y[m.to_numpy()], 1)
    xx = np.linspace(0, 1, 50)
    ax.plot(xx, np.polyval(z, xx), color=c, lw=1.6, zorder=5)
ax.set_xlabel("within-gene percentile rank of score", fontsize=7.6)
ax.set_ylabel("measured damage  (–SGE function score)", fontsize=7.6)
ax.legend(fontsize=6.9, frameon=False, loc="upper left", markerscale=2.2)
ax.set_title("RNU4-2 SGE, U4 family withheld", fontsize=8.5, color=INK, loc="left", pad=5)
ax.text(-0.14, 1.17, "d", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

out = config.FIGURES / "fig2_dissociation.png"
fig.savefig(out, dpi=400); fig.savefig(str(out).replace(".png", ".pdf"))
print(f"wrote {out} and .pdf")
print(f"  a: clinical labels dAUROC {gl.delta:+.3f} [{gl.block_lo:+.3f}, {gl.block_hi:+.3f}] (covers 0); "
      f"measured function drho {pf.delta_rho:+.3f} [{pf.lo:+.3f}, {pf.hi:+.3f}] (excludes 0)")
print(f"  c: structure|constraint {p_sc:+.3f}, constraint|structure {p_cs:+.3f}, "
      f"contacts|constraint {p_gc:+.3f}, constraint|contacts {p_cg:+.3f}")

pd.DataFrame([dict(structure=s, population=p, structure_given_population=a, population_given_structure=b)
              for s, p, a, b in [("snrnavep_v2", "constraint_1_minus_oe", p_sc, p_cs),
                                 ("nb_frac_prot", "constraint_1_minus_oe", p_gc, p_cg)]]
             ).to_csv(R / f"partial_correlations_{LAB}.tsv", sep="\t", index=False)
