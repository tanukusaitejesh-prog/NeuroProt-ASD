"""
Generate Publication Figure: Stratified Performance Across Constraint Regimes
=============================================================================
Compares PLM vs forecASD vs Genomic Baselines in:
1. Highly Constrained Genes (pLI >= 0.9, n=3,186)
2. Intermediate / Unconstrained Genes (pLI < 0.9, n=14,769, 82.2% of genome)
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

FIGURES_DIR = Path("figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Data from our verified stratified evaluation
categories = ["Constrained Subspace\n(pLI ≥ 0.9, n=3,186 genes)", "Unconstrained / Intermediate Subspace\n(pLI < 0.9, n=14,769 genes, 82.2% of Proteome)"]

# AUROCs
models = {
    "forecASD (Network + BrainSpan)": [0.8975, 0.7750],
    "Genomic Baseline (Length / Mutation)": [0.7087, 0.5840],
    "Our Paralog-Aware PLM Model": [0.6322, 0.8252]
}

x = np.arange(len(categories))
width = 0.25

fig, ax = plt.subplots(figsize=(8, 5))

colors = ["#2b5c8f", "#7f7f7f", "#d95f02"]

for i, (name, aurocs) in enumerate(models.items()):
    pos = x + (i - 1) * width
    rects = ax.bar(pos, aurocs, width, label=name, color=colors[i], alpha=0.9, edgecolor="black", linewidth=0.8)
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f"{height:.3f}",
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center", va="bottom", fontsize=9, fontweight="bold")

ax.set_ylabel("Out-of-Fold AUROC", fontsize=11, fontweight="bold")
ax.set_title("Stratified Prioritization: Where Protein Language Models Add Value\n(Evaluating Beyond Saturated Constraint)", fontsize=12, fontweight="bold", pad=12)
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=10, fontweight="bold")
ax.set_ylim(0.4, 1.0)
ax.axhline(0.5, color="black", ls="--", lw=1, alpha=0.5, label="Random Guess (0.50)")
ax.legend(loc="upper left", framealpha=0.9, fontsize=9)
ax.grid(axis="y", linestyle=":", alpha=0.5)

# Annotation box highlighting the key finding
ax.annotate(
    "PLM Outperforms forecASD by +5.0 AUROC points\nacross the 82% of genome where pLI is uninformative",
    xy=(1.12, 0.83), xytext=(0.85, 0.91),
    arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.2),
    fontsize=8.5, bbox=dict(boxstyle="round,pad=0.4", fc="yellow", alpha=0.3, ec="orange")
)

fig.tight_layout()
pdf_path = FIGURES_DIR / "stratified_performance_pli.pdf"
png_path = FIGURES_DIR / "stratified_performance_pli.png"
fig.savefig(pdf_path, bbox_inches="tight")
fig.savefig(png_path, dpi=200, bbox_inches="tight")
print(f"[INFO] Saved stratified plot to {pdf_path} and {png_path}")
