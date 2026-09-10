"""
§7.5 — MAGMA Common-Variant GWAS Enrichment for Novel Candidates
================================================================

Tests whether the "novel candidates" uniquely prioritized by your model
(and missed by forecASD) carry significantly elevated common-variant
association signal in the largest published ASD GWAS (Grove et al. 2019).

Solves the ID mapping problem:
    MAGMA results typically use Entrez IDs or Ensembl IDs (ENSG...).
    This script provides robust HGNC mapping supporting:
    - Entrez Gene ID -> Approved Symbol
    - Ensembl Gene ID -> Approved Symbol
    - Gene Symbol aliases -> Approved Symbol
    Using local forecASD table cache + HGNC complete dataset download.

Inputs:
    - --magma-results: Path to MAGMA gene-level results file
      (e.g., asd.genes.out from PGC Grove et al. 2019 GWAS summary stats,
       or a downloaded table with columns: GENE, P, and optional ZSTAT)
    - --novel-candidates: Path to novel candidates CSV (results/novel_candidates.csv)

Outputs:
    - results/magma_enrichment_summary.json (Mann-Whitney U, odds ratios, CDF stats)
    - figures/magma_pvalue_distribution.pdf (P-value CDF comparison plot)
    - results/novel_candidates_magma_hits.csv (Novel candidates with MAGMA stats)

Usage:
    python 03_magma_enrichment.py --magma-results asd_magma.genes.out --novel-candidates results/novel_candidates.csv
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy.stats import fisher_exact, mannwhitneyu

HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HGNC_LOCAL = "hgnc_complete_set.txt"
FORECASD_LOCAL = "forecASD_table.csv"


class GeneSymbolMapper:
    """Robust gene identifier mapper for Entrez, Ensembl, and Symbol aliases."""

    def __init__(self):
        self.entrez_to_symbol: Dict[str, str] = {}
        self.ensembl_to_symbol: Dict[str, str] = {}
        self.alias_to_symbol: Dict[str, str] = {}
        self._load_mappings()

    def _load_mappings(self):
        # 1. Fast bootstrap from local forecASD table if available
        if os.path.exists(FORECASD_LOCAL):
            try:
                fc = pd.read_csv(FORECASD_LOCAL, usecols=["entrez", "ensembl_string", "symbol"])
                for _, row in fc.iterrows():
                    sym = str(row["symbol"]).strip()
                    if pd.notna(row["entrez"]):
                        self.entrez_to_symbol[str(int(row["entrez"]))] = sym
                    if pd.notna(row["ensembl_string"]):
                        self.ensembl_to_symbol[str(row["ensembl_string"]).strip()] = sym
                print(f"[INFO] Initialized ID mapping with {len(self.entrez_to_symbol)} genes from {FORECASD_LOCAL}")
            except Exception as e:
                print(f"[WARNING] Could not read local {FORECASD_LOCAL}: {e}")

        # 2. Comprehensive HGNC mapping
        if not os.path.exists(HGNC_LOCAL):
            print(f"[INFO] Downloading official HGNC complete set from {HGNC_URL}...")
            try:
                resp = requests.get(HGNC_URL, timeout=60)
                resp.raise_for_status()
                with open(HGNC_LOCAL, "wb") as f:
                    f.write(resp.content)
                print(f"[INFO] Saved {HGNC_LOCAL}")
            except Exception as e:
                print(f"[WARNING] Could not download HGNC set ({e}). Using local forecASD mappings.")
                return

        if os.path.exists(HGNC_LOCAL):
            df = pd.read_csv(HGNC_LOCAL, sep="\t", low_memory=False)
            for _, row in df.iterrows():
                sym = str(row["symbol"]).strip()
                # Entrez
                if pd.notna(row.get("entrez_id")):
                    self.entrez_to_symbol[str(int(row["entrez_id"]))] = sym
                # Ensembl
                if pd.notna(row.get("ensembl_gene_id")):
                    self.ensembl_to_symbol[str(row["ensembl_gene_id"]).strip()] = sym
                # Aliases
                if pd.notna(row.get("alias_symbol")):
                    for alias in str(row["alias_symbol"]).split("|"):
                        self.alias_to_symbol[alias.strip().upper()] = sym
                if pd.notna(row.get("prev_symbol")):
                    for prev in str(row["prev_symbol"]).split("|"):
                        self.alias_to_symbol[prev.strip().upper()] = sym

            print(f"[INFO] Total ID mappings active: {len(self.entrez_to_symbol)} Entrez, {len(self.ensembl_to_symbol)} Ensembl")

    def map_id(self, val) -> Optional[str]:
        """Maps an unknown identifier (Entrez, Ensembl, or Symbol) to approved HGNC symbol."""
        if pd.isna(val):
            return None
        val_str = str(val).strip()

        # Check Entrez (digits)
        if val_str.isdigit():
            return self.entrez_to_symbol.get(val_str, val_str)

        # Check Ensembl
        if val_str.startswith("ENS"):
            base_ens = val_str.split(".")[0]  # strip version if present
            return self.ensembl_to_symbol.get(base_ens, self.ensembl_to_symbol.get(val_str, val_str))

        # Check alias / uppercase symbol
        val_upper = val_str.upper()
        return self.alias_to_symbol.get(val_upper, val_str)


def load_magma_results(path: str, mapper: GeneSymbolMapper) -> pd.DataFrame:
    """Loads MAGMA .genes.out or tabular results file and normalizes columns."""
    # MAGMA standard format is whitespace-delimited
    try:
        df = pd.read_csv(path, sep=r"\s+", low_memory=False)
    except Exception:
        df = pd.read_csv(path, low_memory=False)

    # Normalize column names
    col_map = {}
    for c in df.columns:
        cl = c.strip().upper()
        if cl in ("GENE", "GENE_ID", "ENTREZ", "SYMBOL", "GENE_NAME"):
            col_map[c] = "GENE"
        elif cl in ("P", "PVAL", "P_VALUE", "PVALUE"):
            col_map[c] = "P"
        elif cl in ("ZSTAT", "Z", "Z_SCORE", "ZSCORE"):
            col_map[c] = "ZSTAT"

    df = df.rename(columns=col_map)
    assert "GENE" in df.columns, f"Could not identify gene identifier column in {list(df.columns)}"
    assert "P" in df.columns, f"Could not identify p-value column in {list(df.columns)}"

    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df.dropna(subset=["P"])

    # Map gene IDs to approved human symbols
    print("[INFO] Mapping MAGMA gene identifiers to approved symbols...")
    df["gene_symbol"] = df["GENE"].apply(mapper.map_id)
    df = df.drop_duplicates(subset=["gene_symbol"]).dropna(subset=["gene_symbol"])
    print(f"[INFO] Loaded MAGMA statistics for {len(df)} unique genes")

    return df


def run_magma_enrichment(
    magma_path: str,
    novel_candidates_path: str = "results/novel_candidates.csv",
    nominal_cutoff: float = 0.05,
):
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    os.makedirs("results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)

    mapper = GeneSymbolMapper()
    magma_df = load_magma_results(magma_path, mapper)

    # Load novel candidates
    cand_df = pd.read_csv(novel_candidates_path)
    cand_col = "gene_symbol"
    for c in cand_df.columns:
        if c.lower() in ("gene", "symbol", "gene_symbol", "gene_name"):
            cand_col = c
            break
    novel_genes = set(cand_df[cand_col].astype(str).str.strip())
    print(f"[INFO] Loaded {len(novel_genes)} novel candidate genes from {novel_candidates_path}")

    # Partition MAGMA results
    novel_mask = magma_df["gene_symbol"].isin(novel_genes)
    novel_magma = magma_df[novel_mask].copy()
    bg_magma = magma_df[~novel_mask].copy()

    print(f"[INFO] Novel candidates present in GWAS MAGMA results: {len(novel_magma)} / {len(novel_genes)}")

    if len(novel_magma) == 0:
        print("[ERROR] None of the novel candidates were found in the MAGMA results. Check ID formatting.")
        return

    # 1. Mann-Whitney U test on P-values (one-sided: novel candidates have LOWER p-values)
    mwu_stat, mwu_p = mannwhitneyu(novel_magma["P"], bg_magma["P"], alternative="less")

    # 2. Fisher's exact test at nominal cutoff (P < 0.05)
    a = (novel_magma["P"] < nominal_cutoff).sum()
    b = len(novel_magma) - a
    c = (bg_magma["P"] < nominal_cutoff).sum()
    d = len(bg_magma) - c

    table = np.array([[a, b], [c, d]])
    fisher_or, fisher_p = fisher_exact(table, alternative="greater")

    novel_sig_rate = a / len(novel_magma)
    bg_sig_rate = c / len(bg_magma)

    print("\n" + "=" * 65)
    print("§7.5 MAGMA COMMON-VARIANT GWAS ENRICHMENT TEST (GROVE ET AL. 2019)")
    print("=" * 65)
    print(f"Total genome genes with MAGMA stats:   {len(magma_df)}")
    print(f"Novel candidates with MAGMA stats:     {len(novel_magma)}")
    print(f"Median P-value (Novel Candidates):     {novel_magma['P'].median():.4e}")
    print(f"Median P-value (Genome Background):    {bg_magma['P'].median():.4e}")
    print("-" * 65)
    print(f"Mann-Whitney U Test (Lower P-value):   U = {mwu_stat:.1f}, p = {mwu_p:.4e}")
    print("-" * 65)
    print(f"Enrichment at Nominal P < {nominal_cutoff}:")
    print(f"  Novel Candidate Rate:                {novel_sig_rate*100:.2f}% ({a}/{len(novel_magma)})")
    print(f"  Background Rate:                     {bg_sig_rate*100:.2f}% ({c}/{len(bg_magma)})")
    print(f"  Odds Ratio (OR):                     {fisher_or:.3f}")
    print(f"  Fisher's Exact p-value:              {fisher_p:.4e}")
    print("=" * 65)

    # Save hits
    novel_magma = novel_magma.sort_values("P")
    hits_path = "results/novel_candidates_magma_hits.csv"
    novel_magma.to_csv(hits_path, index=False)
    print(f"[INFO] Saved novel candidate MAGMA statistics to {hits_path}")

    # Plot CDF of -log10(P)
    fig, ax = plt.subplots(figsize=(7, 5))
    sorted_novel = np.sort(-np.log10(novel_magma["P"]))
    sorted_bg = np.sort(-np.log10(bg_magma["P"]))

    ax.plot(sorted_bg, np.linspace(0, 1, len(sorted_bg)),
            label=f"Genome Background (n={len(sorted_bg)})", color="grey", lw=1.5)
    ax.plot(sorted_novel, np.linspace(0, 1, len(sorted_novel)),
            label=f"Novel Candidates (n={len(sorted_novel)})", color="dodgerblue", lw=2.5)

    ax.set_xlabel(r"$-\log_{10}(P_{\mathrm{MAGMA}})$", fontsize=11)
    ax.set_ylabel("Cumulative Fraction", fontsize=11)
    ax.set_title(f"MAGMA Common-Variant Enrichment (MWU p = {mwu_p:.2e})", fontsize=12)
    ax.legend(loc="lower right", framealpha=0.9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig("figures/magma_pvalue_distribution.pdf", bbox_inches="tight")
    fig.savefig("figures/magma_pvalue_distribution.png", dpi=150, bbox_inches="tight")
    print(f"[INFO] Saved figures/magma_pvalue_distribution.pdf")

    summary = {
        "n_novel_candidates_with_magma": int(len(novel_magma)),
        "n_background_genes": int(len(bg_magma)),
        "median_p_novel": float(novel_magma["P"].median()),
        "median_p_background": float(bg_magma["P"].median()),
        "mann_whitney_u_stat": float(mwu_stat),
        "mann_whitney_u_p_value": float(mwu_p),
        "nominal_p_threshold": nominal_cutoff,
        "n_novel_significant": int(a),
        "n_background_significant": int(c),
        "fisher_odds_ratio": float(fisher_or),
        "fisher_p_value": float(fisher_p),
    }
    with open("results/magma_enrichment_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Saved results/magma_enrichment_summary.json")

    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="§7.5: MAGMA Common-Variant GWAS Enrichment")
    parser.add_argument("--magma-results", required=True, help="Path to MAGMA .genes.out or tabular results file")
    parser.add_argument("--novel-candidates", default="results/novel_candidates.csv", help="Path to novel candidates CSV")
    parser.add_argument("--nominal-cutoff", type=float, default=0.05, help="Nominal p-value threshold (default: 0.05)")
    args = parser.parse_args()
    run_magma_enrichment(args.magma_results, args.novel_candidates, args.nominal_cutoff)
