"""Figure 1 for the MuST-VEP manuscript: workflow schematic, state x family coverage, and the G1 gate.

Both comparable papers in the target venue lead with a workflow figure, and it carries real weight with reviewers.
Panel a is the schematic; panels b and c are data, so the figure argues rather than merely illustrates.

  a  pipeline: 11 cryo-EM structures + GENCODE genes -> contact extraction -> family projection -> multi-state
     aggregation and graph smoothing -> mechanism-aware model -> position score, allele score, calibration
  b  which snRNA family is resolved in which state, with nucleotide counts. This is the argument for the
     multi-state representation: U4 is absent from every structure after activation, so no single state can
     serve RNU4-2 and RNU2-2 at once.
  c  gate G1: mean held-out AUROC for all 11 structures, each single state, the nested single-state selection,
     and the degree-preserving shuffled-contact null (200 replicates).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, GREEN, AMBER, PINK = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": SURF})

# canonical cycle order; tri-snRNP is the free particle that docks to form pre-B
ORDER = ["U2-snRNP", "tri-snRNP", "pre-B", "B", "Bact", "Bact-mature", "C*", "P", "minor-preB", "minor-Bact"]
PDB_OF = {"U2-snRNP": "6Y5Q", "tri-snRNP": "3JCR/6QW6", "pre-B": "6QX9", "B": "5O9Z", "Bact": "6FF7",
          "Bact-mature": "5Z56", "C*": "5XJC", "P": "6QDV", "minor-preB": "8Y6O", "minor-Bact": "7DVQ"}

fig = plt.figure(figsize=(7.2, 8.4))
gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 1.15, 1.0], hspace=0.42,
                      left=0.105, right=0.975, top=0.965, bottom=0.085)

# ----------------------------------------------------------------------------------------------- a. workflow
ax = fig.add_subplot(gs[0]); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 34)


def box(x, y, w, h, head, body, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=0.9",
                                fc=fc, ec=INK2, lw=0.7, zorder=2))
    ax.text(x + w / 2, y + h - 1.9, head, ha="center", va="center", fontsize=6.9, zorder=3,
            color=INK, weight="bold")
    ax.text(x + w / 2, y + (h - 3.4) / 2, body, ha="center", va="center", fontsize=6.0, zorder=3,
            color=INK2, linespacing=1.45)


def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=7,
                                 lw=0.75, color=INK2, zorder=1, shrinkA=0, shrinkB=0))


W, XS = 20.5, [0.5, 27.0, 53.0, 79.0]
box(XS[0], 23.0, W, 9.5, "Structures", "11 cryo-EM states\nof the splicing cycle\n4.5 Å heavy-atom cutoff", "#dceafb")
box(XS[0], 12.2, W, 9.0, "Genes", "17 snRNA genes\nGENCODE v27\nparalogue alignment", "#dceafb")
box(XS[0], 1.5, W, 9.0, "Labels", "266 pathogenic variants\n14 sources + ClinVar\nmechanism-aware controls", "#dceafb")

box(XS[1], 18.0, W, 11.0, "Contact extraction", "per nucleotide, per state:\nprotein residues\nmin. distance\n"
                          "snRNA and other RNA", "#fdeadf")
box(XS[1], 4.0, W, 10.0, "Family projection", "paralogues and minor-\nspliceosome analogues\n"
                         "onto shared coordinates", "#fdeadf")

box(XS[2], 10.5, W, 15.5, "Multi-state aggregation", "fraction of states\nwith contact\nmax protein residues\n"
                          "snRNA contact fraction\n1-hop graph smoothing", "#dcf5ec")

box(XS[3], 19.0, W, 10.5, "Model", "L2 logistic regression\n[feat, feat × recessive]\nscorer by nested LOGO",
    "#fdf0d4")
box(XS[3], 2.5, W, 13.5, "Outputs", "position score\n13,280 variants\nallele score where paired\n"
                         "ClinGen calibration\ncallability annotation", "#fce7ef")

arrow(XS[0] + W, 27.7, XS[1], 25.5)
arrow(XS[0] + W, 16.7, XS[1], 21.5)
arrow(XS[0] + W, 6.0, XS[1], 9.0)
arrow(XS[1] + W, 23.5, XS[2], 20.5)
arrow(XS[1] + W, 9.0, XS[2], 16.0)
arrow(XS[2] + W, 18.2, XS[3], 24.2)
arrow(XS[3] + W / 2, 19.0, XS[3] + W / 2, 16.0)
ax.text(0.5, 33.8, "a", fontsize=11, weight="bold", color=INK, va="bottom")

# ------------------------------------------------------------------------- b. state x family coverage (data)
ax = fig.add_subplot(gs[1])
C = pd.read_csv(config.PROCESSED / "contacts_per_structure.tsv", sep="\t")
M = C.pivot_table(index="paralog_family", columns="state", values="ref_n_pos", aggfunc="size").fillna(0)
M = M.reindex(index=["U1", "U2", "U4", "U5", "U6"], columns=[s for s in ORDER if s in M.columns])
im = ax.imshow(M.to_numpy(), cmap="Blues", aspect="auto", vmin=0, vmax=np.nanmax(M.to_numpy()))
ax.set_xticks(range(M.shape[1]))
ax.set_xticklabels([f"{c}\n{PDB_OF.get(c, '')}" for c in M.columns], fontsize=6.6)
ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index, fontsize=8)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M.to_numpy()[i, j]
        ax.text(j, i, f"{int(v)}" if v else "–", ha="center", va="center", fontsize=6.4,
                color="white" if v > 0.6 * np.nanmax(M.to_numpy()) else INK)
ax.axvline(M.columns.get_loc("Bact") - 0.5, color=ORANGE, lw=1.6, ls="--")
ax.text(M.columns.get_loc("Bact") - 0.42, -0.78, "activation: U4 released", color=ORANGE, fontsize=6.8,
        va="bottom", ha="left", weight="bold")
ax.set_title("Nucleotides resolved per snRNA family in each state", fontsize=8.5, color=INK, pad=16, loc="left")
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)
ax.text(-0.085, 1.20, "b", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

# --------------------------------------------------------------------------------- c. gate G1 (data)
ax = fig.add_subplot(gs[2])
S = pd.read_csv(config.RESULTS / "gate_g1_summary_labels_v5c.tsv", sep="\t", index_col=0)
N = pd.read_csv(config.RESULTS / "gate_g1_null_labels_v5c.tsv", sep="\t")
states = [(i.replace("state:", ""), S.loc[i, "dominant"]) for i in S.index if i.startswith("state:")]
states.sort(key=lambda t: t[1])
labels = [s for s, _ in states] + ["nested\nsingle state", "ALL 11\nstructures"]
vals = [v for _, v in states] + [S.loc["NESTED single state", "dominant"], S.loc["ALL", "dominant"]]
cols = [GRID] * len(states) + [AMBER, BLUE]
ax.bar(range(len(vals)), vals, color=cols, edgecolor=INK2, lw=0.5, width=0.74, zorder=3)
v = ax.violinplot([N.dominant.to_numpy()], positions=[len(vals) + 0.9], widths=1.5, showextrema=False)
for b in v["bodies"]:
    b.set_facecolor(ORANGE); b.set_alpha(0.32); b.set_edgecolor(ORANGE); b.set_linewidth(0.7)
ax.plot([len(vals) + 0.9], [N.dominant.mean()], "o", ms=3.4, color=ORANGE, zorder=4)
ax.axhline(S.loc["ALL", "dominant"], color=BLUE, lw=0.8, ls=":", zorder=1)
ax.axhline(0.5, color=INK2, lw=0.6, zorder=1)
ax.set_xticks(list(range(len(vals))) + [len(vals) + 0.9])
ax.set_xticklabels(labels + ["degree-preserving\nnull (200×)"], fontsize=6.3, rotation=38, ha="right")
ax.set_ylabel("mean held-out AUROC, dominant", fontsize=8)
ax.set_ylim(0.40, 0.90)
ax.annotate(f"ALL = {S.loc['ALL', 'dominant']:.3f}\nbeats 200/200 nulls, p = 0.005",
            xy=(len(vals) - 0.3, S.loc["ALL", "dominant"]), xytext=(len(vals) - 4.6, 0.865),
            fontsize=6.8, color=INK, ha="left",
            arrowprops=dict(arrowstyle="-", lw=0.6, color=INK2))
ax.set_title("Gate G1: no single splicing state substitutes for the ensemble", fontsize=8.5, color=INK,
             pad=6, loc="left")
ax.text(-0.085, 1.15, "c", transform=ax.transAxes, fontsize=11, weight="bold", color=INK, va="top")

out = config.FIGURES / "fig1_workflow.png"
fig.savefig(out, dpi=400)
fig.savefig(str(out).replace(".png", ".pdf"))
print(f"wrote {out} and .pdf")
print(f"  panel b: {M.shape[0]} families x {M.shape[1]} states")
print(f"  panel c: {len(states)} single states, null n={len(N)}")
