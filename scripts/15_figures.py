"""Main figures for snRNA-VEP (from results/ tables; no model fitting here). Palette: validated categorical order
(dataviz reference palette, light mode); every figure has a legend or direct labels and a source table.
Usage: SNRNA_LABELS=labels_v5.tsv python scripts/15_figures.py"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
SER = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": SURF})
config.FIGURES.mkdir(exist_ok=True)
UL = {"RNU2-2P": "RNU2-2", "RNU6": "RNU6 (RP)"}

# ---- Fig 1: held-out AUROC by gene and mechanism
res = pd.read_csv(config.result("v2_heldout.tsv"), sep="\t")
nested = pd.read_csv(config.result("v2_nested.tsv"), sep="\t")
res = pd.concat([res, nested.rename(columns={"chosen": "_c"}).assign(score="final")[["unit", "ctx", "score", "auroc", "n_pos"]]])
SC = [("final", "snRNA-VEP (nested)"), ("combo_untrained", "contacts + phyloP"), ("frac_states_protein_contact", "protein contacts"),
      ("cadd_phred", "CADD"), ("phylop447", "phyloP447")]
fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
for ax, ctx in zip(axes, ["dominant", "recessive"]):
    r = res[res.ctx == ctx]
    units = [u for u in ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6", "RNU4ATAC", "RNU12"] if u in set(r.unit)]
    w = 0.8 / len(SC)
    for j, (s, lab) in enumerate(SC):
        v = [r[(r.unit == u) & (r.score == s)].auroc.mean() for u in units] + [r[r.score == s].groupby("unit").auroc.mean().mean()]
        ax.bar(np.arange(len(units) + 1) + (j - 2) * w, v, w * 0.88, color=SER[j], label=lab, zorder=3)
    ax.axhline(0.5, color=INK2, lw=0.8, ls=":", zorder=2)
    ax.yaxis.grid(True, color=GRID, lw=0.6, zorder=0)
    npos = [int(r[(r.unit == u)].n_pos.max()) for u in units]
    ax.set_xticks(np.arange(len(units) + 1), [f"{UL.get(u, u)}\n(n={n})" for u, n in zip(units, npos)] + ["mean"])
    ax.set_ylim(0.3, 1.0)
    ax.set_title(f"{ctx.capitalize()} mechanism (held-out genes)", color=INK, loc="left", fontsize=10)
axes[0].set_ylabel("AUROC")
axes[1].legend(frameon=False, fontsize=8, loc="upper right", ncol=1)
fig.tight_layout(); fig.savefig(config.FIGURES / "fig1_heldout_auroc.png", dpi=300); plt.close(fig)

# ---- Fig 2: SGE zero-shot
zf = config.result("v2_sge_zeroshot.tsv")
z = pd.read_csv(zf if zf.exists() else config.RESULTS / "v2_sge_zeroshot.tsv", sep="\t")
lab = {"v2 trained without RNU4-2": "snRNA-VEP (no RNU4-2 labels)", "v2 trained without U4 family": "snRNA-VEP (no U4-family labels)",
       "nb_frac_prot": "graph-smoothed contacts", "frac_states_protein_contact": "protein contacts", "cadd_phred": "CADD",
       "phylop447": "phyloP447"}
z = z[z.score.isin(lab)].iloc[::-1]
fig, ax = plt.subplots(figsize=(6.2, 2.8))
ax.barh([lab[s] for s in z.score], z.spearman, 0.6, color=SER[0], zorder=3)
for y, v in enumerate(z.spearman):
    ax.text(v + (0.01 if v >= 0 else -0.01), y, f"{v:.2f}", va="center", ha="left" if v >= 0 else "right", color=INK, fontsize=8)
ax.axvline(0, color=INK2, lw=0.8); ax.xaxis.grid(True, color=GRID, lw=0.6, zorder=0)
ax.set_xlabel("Spearman ρ with RNU4-2 saturation-editing damage (485 variants)")
ax.set_xlim(-0.08, 0.55)
fig.tight_layout(); fig.savefig(config.FIGURES / "fig2_sge_zeroshot.png", dpi=300); plt.close(fig)

# ---- Fig 3: mechanism, protein contacts (position level)
P = pd.read_csv(config.result("mechanism_proteins.tsv"), sep="\t")
P = P[P.level == "position"]
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharex=True)
for ax, ctx, col in zip(axes, ["dominant", "recessive"], [SER[0], SER[1]]):
    g = P[P.ctx == ctx].sort_values("p").head(10).iloc[::-1]
    y = np.arange(len(g))
    sig = g.fdr < 0.05
    ax.scatter(g.OR_MH, y, s=np.where(sig, 46, 22), color=np.where(sig, col, "#b9b8b2"), zorder=3, edgecolor=SURF, linewidth=1.5)
    ax.set_yticks(y, [f"{p} ({n})" for p, n in zip(g.protein, g.n_path_contact)])
    ax.axvline(1, color=INK2, lw=0.8, ls=":")
    ax.set_xscale("log"); ax.xaxis.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_title(f"{ctx.capitalize()} variants: contacted proteins", loc="left", fontsize=10, color=INK)
    ax.set_xlabel("Mantel–Haenszel odds ratio (pathogenic vs control)")
fig.text(0.01, 0.01, "Filled: FDR < 0.05; grey: not significant. (n) = pathogenic positions contacting the protein.",
         fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(config.FIGURES / "fig3_mechanism_proteins.png", dpi=300); plt.close(fig)

# ---- Fig 4: RNU4-2 profile
f = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t", low_memory=False,
                usecols=["gene_name", "n_pos", "frac_states_protein_contact", "sge_score", "phylop447"])
g = f[f.gene_name == "RNU4-2"].groupby("n_pos").agg(contact=("frac_states_protein_contact", "mean"),
                                                   sge=("sge_score", "mean"), phylop=("phylop447", "mean"))
sc = lambda s: (s - s.min()) / (s.max() - s.min())
fig, ax = plt.subplots(figsize=(10, 2.8))
ax.axvspan(62, 79, color="#eda100", alpha=0.15, lw=0)
ax.text(70.5, 1.08, "ReNU critical region", ha="center", fontsize=8, color=INK2)
gs = g.dropna(subset=["sge"])
ax.scatter(gs.index, sc(-gs.sge), s=14, color=SER[1], zorder=4, edgecolor=SURF, linewidth=0.8,
           label="SGE damage (scaled; positions measured in the screen)")
ax.plot(g.index, g.contact, color=SER[0], lw=2, label="protein-contact fraction (11 cryo-EM states)")
ax.plot(g.index, sc(g.phylop), color="#9a9993", lw=1.2, label="phyloP447 (scaled)")
ax.set_xlabel("RNU4-2 position (n.)"); ax.set_ylim(0, 1.15); ax.set_xlim(1, g.index.max())
ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.28))
fig.tight_layout(); fig.savefig(config.FIGURES / "fig4_rnu4-2_profile.png", dpi=300); plt.close(fig)
print("written:", sorted(p.name for p in config.FIGURES.glob("fig*.png")))
