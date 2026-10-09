"""Figures 1-3 for the ReNU recurrence manuscript.

  Fig 1  mutation rate across RNU4-2: per-nucleotide Roulette rate with the critical region shaded and patient
         variants overlaid; inside versus outside; rate against patient and population carriers (both flat)
  Fig 2  homopolymer insertion frequency: log-linear scaling with tract length, split by base, and the
         short-read versus long-read comparison that forced the 34x -> 20x correction
  Fig 3  the residual and its parental origin: carriers per allele against mutational expectation; the
         RNU4-2 / RNU4-1 alignment showing conversion cannot produce the allele; parental origin in both cohorts
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, GREEN, AMBER, PINK, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#b9b7b2"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": SURF})
F, R = config.FIGURES, config.RESULTS
CR = (120291825, 120291842)
END = 120291903                     # RNU4-2 is minus strand; n.pos = END - genomic + 1

RO = pd.read_csv(config.PROCESSED / "roulette_rnu4_2.tsv", sep="\t")
C = pd.read_csv(config.CURATION / "renu_variant_counts.tsv", sep="\t", comment="#")
C["ndd"] = C.gel_ndd + C.nongel_ndd
C["pop"] = C.gnomad + C.ukb + C.aou
C["gpos"] = C.key.str.split(":").str[1].astype(int)
C["npos"] = END - C.gpos + 1
SR = pd.read_csv(R / "homopolymer_insertion_rate.tsv", sep="\t")
M = pd.read_csv(R / "homopolymer_rate_lr_vs_sr.tsv", sep="\t")

# ============================================================ FIGURE 1
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw=dict(width_ratios=[2.0, 0.7, 1.1], wspace=0.40))
fig.subplots_adjust(left=0.075, right=0.985, top=0.84, bottom=0.235)

ax = axes[0]
per = RO.groupby("pos").MR.mean().reset_index()
per["npos"] = END - per.pos + 1
per = per.sort_values("npos")
ax.fill_between([END - CR[1] + 1, END - CR[0] + 1], 0, per.MR.max() * 1.12, color=AMBER, alpha=0.18, lw=0,
                zorder=1, label="critical region")
ax.plot(per.npos, per.MR, color=INK2, lw=0.9, zorder=3)
pv = C[C.ndd > 0]
ax.scatter(pv.npos, [per.MR.max() * 1.05] * len(pv), s=np.clip(pv.ndd * 1.6, 6, 70), color=ORANGE,
           zorder=4, lw=0, alpha=0.85, label="patient variants")
ax.set_xlabel("position in RNU4-2 (n.)", fontsize=7.6)
ax.set_ylabel("Roulette relative rate", fontsize=7.6)
ax.set_ylim(0, per.MR.max() * 1.12)
ax.legend(fontsize=6.5, frameon=False, loc="upper left")
ax.set_title("Mutation rate does not peak where variants cluster", fontsize=8.2, color=INK, loc="left", pad=4)
ax.text(-0.10, 1.22, "a", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

ax = axes[1]
ins_, out_ = RO[RO.in_CR].MR, RO[~RO.in_CR].MR
ax.bar([0, 1], [ins_.mean(), out_.mean()], color=[AMBER, GREY], edgecolor=INK2, lw=0.6, width=0.62, zorder=3)
for i, v in enumerate([ins_, out_]):
    ax.errorbar(i, v.mean(), yerr=v.std() / np.sqrt(len(v)), fmt="none", ecolor=INK2, lw=0.9, capsize=3, zorder=4)
    ax.text(i, v.mean() + v.std() / np.sqrt(len(v)) + 0.004, f"{v.mean():.3f}", ha="center",
            fontsize=7, color=INK)
ax.set_xticks([0, 1]); ax.set_xticklabels(["inside\nCR", "outside\nCR"], fontsize=7.2)
ax.set_ylabel("mean relative rate", fontsize=7.6)
ax.set_title("Less mutable", fontsize=8.2, color=INK, loc="left", pad=4)
ax.text(-0.42, 1.22, "b", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

ax = axes[2]
S = C.merge(RO.groupby("key").MR.mean().rename("rate"), on="key", how="left").dropna(subset=["rate"])
S = S[(S.ndd + S["pop"]) > 0]
r1 = spearmanr(S.rate, S.ndd); r2 = spearmanr(S.rate, S["pop"])
ax.scatter(S.rate, S.ndd, s=18, color=ORANGE, lw=0, alpha=0.8, label=f"patients  $\\rho$={r1.statistic:+.2f}")
ax.scatter(S.rate, S["pop"], s=18, color=BLUE, lw=0, alpha=0.8, marker="s",
           label=f"population  $\\rho$={r2.statistic:+.2f}")
ax.set_yscale("symlog", linthresh=1)
ax.set_xlabel("Roulette relative rate", fontsize=7.6)
ax.set_ylabel("carriers", fontsize=7.6)
ax.legend(fontsize=6.4, frameon=False, loc="upper right")
ax.set_title("Rate predicts neither", fontsize=8.2, color=INK, loc="left", pad=4)
ax.text(-0.26, 1.22, "c", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")
fig.savefig(F / "recur_fig1_rate.png", dpi=400); fig.savefig(F / "recur_fig1_rate.pdf")
print(f"fig1: rate in CR {ins_.mean():.3f} vs out {out_.mean():.3f}; "
      f"rho patients {r1.statistic:+.2f} (p {r1.pvalue:.2f}), population {r2.statistic:+.2f} (p {r2.pvalue:.2f})")

# ============================================================ FIGURE 2
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw=dict(wspace=0.38))
fig.subplots_adjust(left=0.085, right=0.985, top=0.84, bottom=0.235)

# fit on every A/T (base, length) row, as reported in the text - not on per-length means, which gives 1.86
AT = SR[SR.base.isin(["A", "T"])]
ax = axes[0]
for b, mk in [("A", "o"), ("T", "s")]:
    d = AT[AT.base == b]
    ax.scatter(d.length, d.rel_to_snv, s=24, color=BLUE, marker=mk, zorder=4, lw=0, alpha=0.8, label=f"{b} tracts")
k = np.polyfit(AT.length, np.log10(AT.rel_to_snv), 1)
_r2 = np.corrcoef(AT.length, np.log10(AT.rel_to_snv))[0, 1] ** 2
xx = np.linspace(AT.length.min(), AT.length.max(), 50)
ax.plot(xx, 10 ** np.polyval(k, xx), color=INK2, lw=1.0, ls="--", zorder=3,
        label=f"{10 ** k[0]:.2f}× per base, R²={_r2:.2f}")
ax.axhline(1.0, color=INK2, lw=0.7, zorder=2)
t4 = AT[AT.length == 4].rel_to_snv.mean()
ax.scatter([4], [t4], s=80, facecolor="none", edgecolor=ORANGE, lw=1.6, zorder=5)
ax.annotate("T₄: the ReNU context", xy=(4.25, t4), xytext=(6.2, 0.115), fontsize=6.8, color=ORANGE,
            arrowprops=dict(arrowstyle="-", lw=0.7, color=ORANGE))
ax.set_yscale("log"); ax.set_xlabel("homopolymer tract length", fontsize=7.6)
ax.set_ylabel("insertion freq. ÷ substitution freq.", fontsize=7.6)
ax.legend(fontsize=6.6, frameon=False, loc="upper left")
ax.set_title("Log-linear slippage scaling", fontsize=8.2, color=INK, loc="left", pad=4)
ax.text(-0.22, 1.22, "a", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

ax = axes[1]
for b, c in [("A", GREEN), ("T", BLUE), ("C", AMBER), ("G", PINK)]:
    d = SR[SR.base == b]
    ax.plot(d.length, d.rel_to_snv, "-o", ms=3, lw=1.0, color=c, label=b)
ax.axhline(1.0, color=INK2, lw=0.7)
ax.set_yscale("log"); ax.set_xlabel("tract length", fontsize=7.6)
ax.set_ylabel("relative insertion freq.", fontsize=7.6)
ax.legend(fontsize=6.6, frameon=False, ncol=2, loc="upper left")
ax.set_title("By base", fontsize=8.2, color=INK, loc="left", pad=4)
ax.text(-0.26, 1.22, "b", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

ax = axes[2]
MT = M[M.base.isin(["A", "T"])].groupby("length")[["rel_to_snv_sr", "rel_to_snv_lr"]].mean().reset_index()
w = 0.36
ax.bar(MT.length - w / 2, MT.rel_to_snv_sr, w, color=GREY, edgecolor=INK2, lw=0.5, label="short read (gnomAD)",
       zorder=3)
ax.bar(MT.length + w / 2, MT.rel_to_snv_lr, w, color=GREEN, edgecolor=INK2, lw=0.5, label="long read (HPRC)",
       zorder=3)
ax.axhline(1.0, color=INK2, lw=0.7, zorder=2)
ax.set_yscale("log"); ax.set_xlabel("tract length", fontsize=7.6)
ax.set_ylabel("relative insertion freq.", fontsize=7.6)
ax.legend(fontsize=6.4, frameon=False, loc="upper left")
ax.set_title("Short-read under-detection", fontsize=8.2, color=INK, loc="left", pad=4)
ax.text(-0.26, 1.22, "c", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")
fig.savefig(F / "recur_fig2_homopolymer.png", dpi=400); fig.savefig(F / "recur_fig2_homopolymer.pdf")
r4 = MT[MT.length == 4].iloc[0]
print(f"fig2: slope {10 ** k[0]:.2f}x/base; T/A length 4 short {r4.rel_to_snv_sr:.3f} long {r4.rel_to_snv_lr:.3f}")

# ============================================================ FIGURE 3
fig = plt.figure(figsize=(7.2, 2.7))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.45, 0.8], wspace=0.46,
                      left=0.085, right=0.985, top=0.83, bottom=0.21)

# panel a: the position-matched comparison. All four variants sit at one nucleotide, so assay coverage, calling
# behaviour and ascertainment are held constant; the two other insertions are the calling-bias control.
ax = fig.add_subplot(gs[0])
exp = 0.403
FOCAL = 120291839
f = C[C.gpos == FOCAL].sort_values("ndd", ascending=False)
labs = ["insT", "insG", "A>G", "insC"]
vals = [int(f[f.hgvs == h].ndd.iloc[0]) for h in ["n.64_65insT", "n.64_65insG", "n.65A>G", "n.64_65insC"]]
cols = [ORANGE, GREY, GREY, GREY]
ax.bar(range(4), vals, color=cols, edgecolor=INK2, lw=0.6, width=0.62, zorder=3)
for i, v in enumerate(vals):
    ax.text(i, v + 2.2, str(v), ha="center", fontsize=7.4, color=INK, weight="bold" if i == 0 else "normal")
ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=7.0)
ax.set_xlabel("variant at n.64_65 / n.65", fontsize=7.0)
ax.set_ylabel("patient carriers", fontsize=7.6)
ax.set_ylim(0, max(vals) * 1.34)
ax.text(0.68, 0.76, "all four variants at\nchr12:120,291,839", transform=ax.transAxes, ha="center", va="center",
        fontsize=6.6, color=INK2)
ax.text(0.68, 0.48, "vs the substitution\nat this base:\n89/91 obs, 28.7% exp\np = 1.3×10$^{-45}$",
        transform=ax.transAxes, ha="center", va="center", fontsize=6.4, color=ORANGE, weight="bold")
ax.set_title("Position-matched", fontsize=8.2, color=INK, loc="left", pad=4)
ax.text(-0.28, 1.19, "a", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

ax = fig.add_subplot(gs[1]); ax.axis("off")
COMP = str.maketrans("ACGTN", "TGCAN")
RF = pd.read_csv(config.PROCESSED / "reference_core_loci.tsv.gz", sep="\t").set_index("gene_name")


def mature(g):
    s = RF.loc[g].seq.upper()
    return s[2000:len(s) - 2000].translate(COMP)[::-1]


a2, a1 = mature("RNU4-2"), mature("RNU4-1")
lo, hi = 59, 80                                   # n.60-80, spanning the critical region and the tract
ax.set_xlim(-0.6, hi - lo + 0.6); ax.set_ylim(-1.6, 3.4)
for j, i in enumerate(range(lo, hi)):
    same = a2[i] == a1[i]
    intract = 64 <= i <= 67                       # n.65-68 = the A4 run (T4 on the genomic plus strand)
    for row, seq, lab in [(2, a2, "RNU4-2"), (1, a1, "RNU4-1")]:
        ax.text(j, row, seq[i], ha="center", va="center", fontsize=7.4, family="monospace",
                color=ORANGE if intract else (INK if same else BLUE),
                weight="bold" if intract else "normal")
    ax.text(j, 0.25, "|" if same else "×", ha="center", va="center", fontsize=6.4,
            color=GREY if same else BLUE)
    if i % 5 == 4:
        ax.text(j, -0.55, str(i + 1), ha="center", va="center", fontsize=5.8, color=INK2)
ax.text(-0.9, 2, "RNU4-2", ha="right", va="center", fontsize=7.2, color=INK)
ax.text(-0.9, 1, "RNU4-1", ha="right", va="center", fontsize=7.2, color=INK)
ax.plot([64 - lo - 0.45, 67 - lo + 0.45], [2.72, 2.72], color=ORANGE, lw=2.2, solid_capstyle="butt")
ax.text((64 + 67) / 2 - lo, 2.95, "the tract n.64_65insT extends — identical in both", ha="center",
        fontsize=6.6, color=ORANGE, weight="bold")
ax.text((hi + lo) / 2 - lo, -1.25, "97.2% identical; all 4 differences lie OUTSIDE the critical region\n"
                                   "→ conversion explains 0/11 alleles, 0/114 carriers",
        ha="center", fontsize=6.8, color=INK)
ax.set_title("Conversion from RNU4-1 cannot produce it", fontsize=8.2, color=INK, loc="left", pad=4)
ax.text(-0.14, 1.19, "b", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

ax = fig.add_subplot(gs[2])
coh = [("discovery\n(Chen 2024)", 54, 0), ("French\n(2025)", 47, 4)]
for i, (lab, mat, pat) in enumerate(coh):
    ax.bar(i, mat, color=PINK, edgecolor=INK2, lw=0.6, width=0.6, zorder=3)
    ax.bar(i, pat, bottom=mat, color=BLUE, edgecolor=INK2, lw=0.6, width=0.6, zorder=3)
    ax.text(i, mat / 2, f"{mat}", ha="center", va="center", fontsize=7.4, color="white", weight="bold")
    if pat:
        ax.text(i, mat + pat + 1.6, f"{pat} pat.\n(none insT)", ha="center", fontsize=6.2, color=BLUE)
ax.set_xticks(range(len(coh))); ax.set_xticklabels([c[0] for c in coh], fontsize=6.9)
ax.set_ylabel("informative cases", fontsize=7.6); ax.set_ylim(0, 64)
ax.text(0.02, 0.96, "maternal", transform=ax.transAxes, fontsize=7, color=PINK, weight="bold", va="top")
ax.set_title("Parental origin", fontsize=8.2, color=INK, loc="center", pad=4)
ax.text(-0.34, 1.19, "c", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")
fig.savefig(F / "recur_fig3_residual.png", dpi=400); fig.savefig(F / "recur_fig3_residual.pdf")
print(f"fig3: position-matched {vals[0]}/{sum(vals)} at chr12:{FOCAL:,}; others {vals[1:]}")
print(f"wrote recur_fig1_rate, recur_fig2_homopolymer, recur_fig3_residual (.png and .pdf) to {F}")
