"""
§7.3 — Head-to-Head Against forecASD
=====================================

Downloads forecASD genome-wide gene scores, merges with your model's scores,
computes rank correlation, and defines the "novel candidate" gene set:
genes your model ranks highly that forecASD does NOT, excluding known SFARI genes.

Inputs:
    - forecASD_table.csv  (auto-downloaded from GitHub if missing)
    - genome_wide_scores.csv  (YOUR model's output from §6; columns: gene_symbol, score)

Outputs:
    - results/forecasd_comparison.csv     (merged table)
    - results/novel_candidates.csv        (genes in ours_only set)
    - figures/forecasd_scatter.pdf        (quadrant scatter plot)
    - Printed summary stats

Usage:
    python 01_forecasd_comparison.py --our-scores genome_wide_scores.csv
"""

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# ── Configuration ──────────────────────────────────────────────────
FORECASD_URL = (
    "https://raw.githubusercontent.com/LeoBman/forecASD/master/"
    "suppl_data/forecASD_table.csv"
)
FORECASD_LOCAL = "forecASD_table.csv"

# Actual column names (verified from the file):
#   symbol     = gene symbol (human, e.g. "SCN2A")
#   forecASD   = forecASD probability score (0–1, higher = more likely ASD)
#   SFARI_listed = TRUE/FALSE
#   SFARI_score  = 1, 2, 3, S, or NA
FORECASD_GENE_COL = "symbol"
FORECASD_SCORE_COL = "forecASD"

# How many top genes define each model's "high-confidence" set
TOP_N = 500


def download_forecasd(dest: str = FORECASD_LOCAL) -> str:
    """Download forecASD scores if not already present."""
    if os.path.exists(dest):
        print(f"[INFO] Using cached {dest}")
        return dest
    print(f"[INFO] Downloading forecASD scores from GitHub...")
    import urllib.request
    urllib.request.urlretrieve(FORECASD_URL, dest)
    print(f"[INFO] Saved to {dest}")
    return dest


def load_forecasd(path: str) -> pd.DataFrame:
    """Load and validate forecASD table."""
    df = pd.read_csv(path)
    assert FORECASD_GENE_COL in df.columns, (
        f"Expected column '{FORECASD_GENE_COL}', got: {list(df.columns)}"
    )
    assert FORECASD_SCORE_COL in df.columns, (
        f"Expected column '{FORECASD_SCORE_COL}', got: {list(df.columns)}"
    )
    # Convert forecASD score to float (should already be, but be safe)
    df[FORECASD_SCORE_COL] = pd.to_numeric(df[FORECASD_SCORE_COL], errors="coerce")
    print(f"[INFO] Loaded {len(df)} genes from forecASD")
    print(f"       Score range: [{df[FORECASD_SCORE_COL].min():.4f}, "
          f"{df[FORECASD_SCORE_COL].max():.4f}]")
    print(f"       SFARI-listed: {(df['SFARI_listed'] == True).sum()}")
    return df


def load_our_scores(path: str) -> pd.DataFrame:
    """Load your model's genome-wide scores.

    Expected columns: gene_symbol, score
    Adjust column names here if your output uses different headers.
    """
    df = pd.read_csv(path)
    # Normalize column names — adapt these if your output differs
    col_map = {}
    for c in df.columns:
        cl = c.lower().strip()
        if cl in ("gene_symbol", "gene", "symbol", "gene_name"):
            col_map[c] = "gene_symbol"
        elif cl in ("score", "prob", "probability", "asd_score", "our_score", "prediction") or cl.endswith("_score"):
            col_map[c] = "score"
    df = df.rename(columns=col_map)
    assert "gene_symbol" in df.columns, (
        f"Could not find gene symbol column. Columns: {list(df.columns)}"
    )
    assert "score" in df.columns, (
        f"Could not find score column. Columns: {list(df.columns)}"
    )
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score"])
    print(f"[INFO] Loaded {len(df)} genes from your model")
    return df


def main(our_scores_path: str, top_n: int = TOP_N):
    os.makedirs("results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)

    # ── 1. Load both ranked lists ──────────────────────────────────
    forecasd_path = download_forecasd()
    fc = load_forecasd(forecasd_path)
    ours = load_our_scores(our_scores_path)

    # ── 2. Merge on gene symbol ────────────────────────────────────
    merged = fc.merge(
        ours,
        left_on=FORECASD_GENE_COL,
        right_on="gene_symbol",
        how="inner",
    )
    print(f"\n[RESULT] Overlapping genes: {len(merged)}")

    # ── 3. Genome-wide Spearman correlation ────────────────────────
    mask = merged[FORECASD_SCORE_COL].notna() & merged["score"].notna()
    rho, p_rho = spearmanr(
        merged.loc[mask, FORECASD_SCORE_COL],
        merged.loc[mask, "score"],
    )
    # Ensure safe stdout on Windows
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    print(f"[RESULT] Spearman rho = {rho:.4f} (p = {p_rho:.2e})")

    # ── 4. Define gene sets ────────────────────────────────────────
    # Known ASD genes: SFARI-listed genes in the forecASD table itself
    known_asd = set(
        fc.loc[fc["SFARI_listed"] == True, FORECASD_GENE_COL]
    )
    # If you have a separate SFARI gene list, load it here instead:
    # known_asd = set(pd.read_csv("sfari_genes.csv")["gene_symbol"])
    print(f"[INFO] Known SFARI genes in forecASD table: {len(known_asd)}")

    our_top = set(merged.nlargest(top_n, "score")[FORECASD_GENE_COL])
    fc_top = set(merged.nlargest(top_n, FORECASD_SCORE_COL)[FORECASD_GENE_COL])

    both_high = our_top & fc_top
    ours_only = our_top - fc_top - known_asd  # NOVEL CANDIDATES
    fc_only = fc_top - our_top - known_asd
    neither = set(merged[FORECASD_GENE_COL]) - our_top - fc_top

    print(f"\n[RESULT] Top-{top_n} gene set comparison:")
    print(f"  Both high (agreement):     {len(both_high)}")
    print(f"  Ours only (novel):         {len(ours_only)}  ← candidates for IMPC validation")
    print(f"  forecASD only:             {len(fc_only)}")
    print(f"  Neither:                   {len(neither)}")

    # ── 5. Save results ───────────────────────────────────────────
    merged_out = merged[[FORECASD_GENE_COL, FORECASD_SCORE_COL, "score"]].copy()
    merged_out.columns = ["gene_symbol", "forecasd_score", "our_score"]
    merged_out["in_our_top"] = merged_out["gene_symbol"].isin(our_top)
    merged_out["in_forecasd_top"] = merged_out["gene_symbol"].isin(fc_top)
    merged_out["known_sfari"] = merged_out["gene_symbol"].isin(known_asd)
    merged_out["category"] = "neither"
    merged_out.loc[merged_out["gene_symbol"].isin(both_high), "category"] = "both_high"
    merged_out.loc[merged_out["gene_symbol"].isin(ours_only), "category"] = "ours_only_novel"
    merged_out.loc[merged_out["gene_symbol"].isin(fc_only), "category"] = "forecasd_only"
    merged_out = merged_out.sort_values("our_score", ascending=False)
    merged_out.to_csv("results/forecasd_comparison.csv", index=False)

    novel_df = merged_out[merged_out["category"] == "ours_only_novel"].copy()
    novel_df.to_csv("results/novel_candidates.csv", index=False)
    print(f"\n[INFO] Saved results/forecasd_comparison.csv ({len(merged_out)} genes)")
    print(f"[INFO] Saved results/novel_candidates.csv ({len(novel_df)} novel candidates)")

    # ── 6. Scatter plot ────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 8))

    # Background
    bg = merged_out[merged_out["category"] == "neither"]
    ax.scatter(bg["forecasd_score"], bg["our_score"],
               s=2, alpha=0.15, c="grey", label=f"Neither top-{top_n}", zorder=1)

    # Known SFARI
    sfari_mask = merged_out["known_sfari"]
    ax.scatter(merged_out.loc[sfari_mask, "forecasd_score"],
               merged_out.loc[sfari_mask, "our_score"],
               s=18, alpha=0.6, c="crimson", edgecolors="darkred",
               linewidths=0.3, label=f"Known SFARI genes", zorder=3)

    # Novel (ours only)
    novel_mask = merged_out["category"] == "ours_only_novel"
    ax.scatter(merged_out.loc[novel_mask, "forecasd_score"],
               merged_out.loc[novel_mask, "our_score"],
               s=18, alpha=0.7, c="dodgerblue", edgecolors="navy",
               linewidths=0.3, label=f"Novel (ours only, n={novel_mask.sum()})", zorder=4)

    # Both high
    both_mask = merged_out["category"] == "both_high"
    ax.scatter(merged_out.loc[both_mask, "forecasd_score"],
               merged_out.loc[both_mask, "our_score"],
               s=18, alpha=0.6, c="gold", edgecolors="darkorange",
               linewidths=0.3, label=f"Both top-{top_n} (n={both_mask.sum()})", zorder=2)

    ax.set_xlabel("forecASD score", fontsize=12)
    ax.set_ylabel("Our model score", fontsize=12)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax.set_title(f"Genome-wide comparison  (Spearman ρ = {rho:.3f}, p = {p_rho:.1e})",
                 fontsize=13)

    fig.tight_layout()
    fig.savefig("figures/forecasd_scatter.pdf", dpi=300, bbox_inches="tight")
    fig.savefig("figures/forecasd_scatter.png", dpi=150, bbox_inches="tight")
    print(f"[INFO] Saved figures/forecasd_scatter.pdf")

    return novel_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="§7.3: Head-to-head vs forecASD")
    parser.add_argument(
        "--our-scores", required=True,
        help="Path to your genome-wide scores CSV (columns: gene_symbol, score)"
    )
    parser.add_argument(
        "--top-n", type=int, default=TOP_N,
        help=f"Top N genes per method to compare (default: {TOP_N})"
    )
    args = parser.parse_args()
    main(args.our_scores, args.top_n)
