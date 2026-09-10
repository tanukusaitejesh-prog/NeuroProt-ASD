"""
§7.4 — IMPC Mouse Phenotype Enrichment for Novel Candidates
============================================================

Tests whether the "novel candidates" uniquely prioritized by your model
(and missed by forecASD) are statistically enriched for abnormal
nervous system and behavioral phenotypes in standardized knockout mice
from the International Mouse Phenotyping Consortium (IMPC).

Data Sources (Live & Automated):
    1. IMPC Solr API:
       - Tested universe: 'statistical-result' core (all ~9,900+ genes tested)
       - Neuro/behavior phenotypes: 'genotype-phenotype' core
         (MP:0005386 'behavior/neurological phenotype' OR
          MP:0003631 'nervous system phenotype')
    2. MGI Orthology:
       - The Jackson Laboratory / MGI HOM_MouseHumanSequence.rpt
         (maps mouse marker symbols to approved human gene symbols)

Inputs:
    - results/novel_candidates.csv (from 01_forecasd_comparison.py)
      or any CSV with a 'gene_symbol' column.

Outputs:
    - results/impc_enrichment_summary.json  (stats, contingency table, OR, p-value)
    - results/novel_candidates_impc_hits.csv (novel genes with specific mouse phenotypes)

Usage:
    python 02_impc_enrichment.py --novel-candidates results/novel_candidates.csv
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, Set, Tuple

import numpy as np
import pandas as pd
import requests
from scipy.stats import fisher_exact

# ── Endpoints ──────────────────────────────────────────────────────
IMPC_SOLR_STAT = "https://www.ebi.ac.uk/mi/impc/solr/statistical-result/select"
IMPC_SOLR_GENO = "https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select"
MGI_ORTHOLOG_URL = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
MGI_LOCAL_CACHE = "HOM_MouseHumanSequence.rpt"


def get_mouse_to_human_mapping() -> Tuple[Dict[str, str], Dict[str, Set[str]]]:
    """Downloads MGI orthology report and builds mouse-to-human orthology maps.

    Returns:
        mouse_to_human: dict mapping UPPERCASE mouse symbol -> human symbol
        human_to_mouse: dict mapping UPPERCASE human symbol -> set of mouse symbols
    """
    if not os.path.exists(MGI_LOCAL_CACHE):
        print(f"[INFO] Downloading MGI ortholog table from {MGI_ORTHOLOG_URL}...")
        headers = {"User-Agent": "Mozilla/5.0"}
        resp = requests.get(MGI_ORTHOLOG_URL, headers=headers, timeout=60)
        resp.raise_for_status()
        with open(MGI_LOCAL_CACHE, "wb") as f:
            f.write(resp.content)
        print(f"[INFO] Saved {MGI_LOCAL_CACHE}")
    else:
        print(f"[INFO] Using cached {MGI_LOCAL_CACHE}")

    df = pd.read_csv(MGI_LOCAL_CACHE, sep="\t", low_memory=False)
    human_df = df[df["Common Organism Name"] == "human"][["DB Class Key", "Symbol"]].rename(
        columns={"Symbol": "human_symbol"}
    )
    mouse_df = df[df["Common Organism Name"] == "mouse, laboratory"][["DB Class Key", "Symbol"]].rename(
        columns={"Symbol": "mouse_symbol"}
    )

    ortho = mouse_df.merge(human_df, on="DB Class Key")

    mouse_to_human = {}
    human_to_mouse = {}
    for _, row in ortho.iterrows():
        m_sym = str(row["mouse_symbol"]).strip().upper()
        h_sym = str(row["human_symbol"]).strip()
        mouse_to_human[m_sym] = h_sym

        h_upper = h_sym.upper()
        if h_upper not in human_to_mouse:
            human_to_mouse[h_upper] = set()
        human_to_mouse[h_upper].add(row["mouse_symbol"])

    print(f"[INFO] Loaded ortholog mappings for {len(mouse_to_human)} mouse genes -> human")
    return mouse_to_human, human_to_mouse


def fetch_impc_genes() -> Tuple[Set[str], Set[str], pd.DataFrame]:
    """Queries IMPC Solr API for:

    1. All genes tested in IMPC experiments (true background)
    2. All genes with statistically significant nervous system or behavior phenotypes
    3. Detailed phenotype hit annotations for reporting

    Returns:
        all_tested_mouse: set of mouse symbols
        neuro_mouse: set of mouse symbols with neuro/behavior phenotype
        neuro_details_df: DataFrame of individual phenotype assertions
    """
    print("[INFO] Fetching tested genes universe from IMPC Solr (statistical-result core)...")
    r1 = requests.get(
        IMPC_SOLR_STAT,
        params={
            "q": "*:*",
            "rows": 0,
            "facet": "true",
            "facet.field": "marker_symbol",
            "facet.limit": -1,
            "facet.mincount": 1,
            "wt": "json",
        },
        timeout=60,
    )
    r1.raise_for_status()
    f1 = r1.json()["facet_counts"]["facet_fields"]["marker_symbol"]
    # Solr facet array alternates [val1, count1, val2, count2, ...]
    all_tested_mouse = {f1[i].strip() for i in range(0, len(f1), 2)}
    print(f"[INFO] Total IMPC tested mouse genes: {len(all_tested_mouse)}")

    print("[INFO] Fetching neuro/behavior phenotype hits from IMPC Solr (genotype-phenotype core)...")
    neuro_query = 'top_level_mp_term_id:"MP:0005386" OR top_level_mp_term_id:"MP:0003631"'
    r2 = requests.get(
        IMPC_SOLR_GENO,
        params={
            "q": neuro_query,
            "rows": 0,
            "facet": "true",
            "facet.field": "marker_symbol",
            "facet.limit": -1,
            "facet.mincount": 1,
            "wt": "json",
        },
        timeout=60,
    )
    r2.raise_for_status()
    f2 = r2.json()["facet_counts"]["facet_fields"]["marker_symbol"]
    neuro_mouse = {f2[i].strip() for i in range(0, len(f2), 2)}
    print(f"[INFO] IMPC mouse genes with nervous system / behavioral phenotypes: {len(neuro_mouse)}")

    # Also fetch detailed phenotypes for hit annotations
    print("[INFO] Fetching detailed phenotype terms for neuro hits...")
    r3 = requests.get(
        IMPC_SOLR_GENO,
        params={
            "q": neuro_query,
            "fl": "marker_symbol,mp_term_id,mp_term_name,p_value,top_level_mp_term_name",
            "rows": 100000,
            "wt": "json",
        },
        timeout=90,
    )
    r3.raise_for_status()
    docs = r3.json()["response"]["docs"]
    neuro_details_df = pd.DataFrame(docs)
    print(f"[INFO] Retrieved {len(neuro_details_df)} phenotype assertion records")

    return all_tested_mouse, neuro_mouse, neuro_details_df


def run_enrichment_test(
    novel_candidates_csv: str,
    gene_col: str = "gene_symbol",
):
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    os.makedirs("results", exist_ok=True)

    # 1. Load novel candidates
    cand_df = pd.read_csv(novel_candidates_csv)
    if gene_col not in cand_df.columns:
        # Try finding standard gene symbol column
        for c in cand_df.columns:
            if c.lower() in ("gene", "symbol", "gene_symbol", "gene_name"):
                gene_col = c
                break
    novel_genes = set(cand_df[gene_col].astype(str).str.strip())
    print(f"\n[INFO] Loaded {len(novel_genes)} novel candidate genes from {novel_candidates_csv}")

    # 2. Get orthology
    mouse_to_human, human_to_mouse = get_mouse_to_human_mapping()

    # 3. Get IMPC data
    all_tested_mouse, neuro_mouse, neuro_details_df = fetch_impc_genes()

    # Map mouse sets to human gene symbols
    all_impc_tested_human = set()
    for m in all_tested_mouse:
        h = mouse_to_human.get(m.upper())
        if h:
            all_impc_tested_human.add(h)

    impc_neuro_human = set()
    for m in neuro_mouse:
        h = mouse_to_human.get(m.upper())
        if h:
            impc_neuro_human.add(h)

    print(f"\n[INFO] Mapped to Human:")
    print(f"  Total IMPC tested genes (human orthologs): {len(all_impc_tested_human)}")
    print(f"  IMPC neuro/behavior phenotype genes (human orthologs): {len(impc_neuro_human)}")

    # 4. Construct 2x2 Contingency Table
    # Fair comparison: restricted to genes that were actually tested by IMPC
    novel_tested = novel_genes & all_impc_tested_human
    novel_neuro = novel_tested & impc_neuro_human

    non_novel_tested = all_impc_tested_human - novel_genes
    non_novel_neuro = non_novel_tested & impc_neuro_human

    a = len(novel_neuro)                         # Novel & Neuro Hit
    b = len(novel_tested) - a                    # Novel & No Neuro Hit
    c = len(non_novel_neuro)                     # Non-novel & Neuro Hit
    d = len(non_novel_tested) - c                # Non-novel & No Neuro Hit

    table = np.array([[a, b], [c, d]])
    odds_ratio, p_value = fisher_exact(table, alternative="greater")

    novel_rate = (a / len(novel_tested)) if len(novel_tested) > 0 else 0.0
    bg_rate = (c / len(non_novel_tested)) if len(non_novel_tested) > 0 else 0.0

    print("\n" + "=" * 65)
    print("§7.4 IMPC MOUSE PHENOTYPE ENRICHMENT TEST (ONE-SIDED FISHER'S EXACT)")
    print("=" * 65)
    print(f"Novel candidate genes submitted:         {len(novel_genes)}")
    print(f"Novel candidate genes tested by IMPC:     {len(novel_tested)} ({len(novel_tested)/len(novel_genes)*100:.1f}%)")
    print(f"Novel candidates with neuro/behavior hit: {a}")
    print("-" * 65)
    print(f"Contingency Table:")
    print(f"                          Neuro Phenotype (+)   Neuro Phenotype (-)")
    print(f"  Novel Candidates:       {a:<21} {b:<21}")
    print(f"  Genome Background:      {c:<21} {d:<21}")
    print("-" * 65)
    print(f"Novel Candidate Hit Rate:  {novel_rate*100:.2f}% ({a}/{len(novel_tested)})")
    print(f"Background Hit Rate:       {bg_rate*100:.2f}% ({c}/{len(non_novel_tested)})")
    print(f"Odds Ratio (OR):           {odds_ratio:.3f}")
    print(f"Fisher's Exact p-value:    {p_value:.4e}")
    print("=" * 65)

    # 5. Extract specific novel candidate hits and their mouse phenotypes
    # Map back to mouse symbols
    novel_hits_records = []
    for h_gene in sorted(novel_neuro):
        mouse_symbols = human_to_mouse.get(h_gene.upper(), set())
        for m_sym in mouse_symbols:
            sub = neuro_details_df[neuro_details_df["marker_symbol"] == m_sym]
            if not sub.empty:
                for _, row in sub.iterrows():
                    novel_hits_records.append({
                        "human_gene_symbol": h_gene,
                        "mouse_marker_symbol": m_sym,
                        "mp_term_id": row.get("mp_term_id", ""),
                        "mp_term_name": row.get("mp_term_name", ""),
                        "p_value": row.get("p_value", np.nan),
                        "top_level_category": str(row.get("top_level_mp_term_name", "")),
                    })

    hits_df = pd.DataFrame(novel_hits_records)
    if not hits_df.empty:
        hits_df = hits_df.sort_values(["human_gene_symbol", "p_value"])
        hits_path = "results/novel_candidates_impc_hits.csv"
        hits_df.to_csv(hits_path, index=False)
        print(f"\n[INFO] Saved {len(hits_df)} phenotype records for {len(novel_neuro)} novel genes to {hits_path}")
        print("\nTop 10 Phenotype Hits for Novel Candidates:")
        print(hits_df[["human_gene_symbol", "mp_term_name", "p_value"]].head(10).to_string(index=False))

    # 6. Save summary JSON
    summary = {
        "n_novel_candidates_total": len(novel_genes),
        "n_novel_candidates_impc_tested": len(novel_tested),
        "n_novel_candidates_neuro_hit": a,
        "n_novel_candidates_no_hit": b,
        "n_background_impc_tested": len(non_novel_tested),
        "n_background_neuro_hit": c,
        "n_background_no_hit": d,
        "novel_hit_rate": novel_rate,
        "background_hit_rate": bg_rate,
        "odds_ratio": float(odds_ratio),
        "fisher_p_value": float(p_value),
        "significant_at_05": bool(p_value < 0.05),
        "tested_mp_categories": ["MP:0005386 (behavior/neurological)", "MP:0003631 (nervous system)"],
        "validated_novel_genes": sorted(list(novel_neuro)),
    }
    with open("results/impc_enrichment_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Saved results/impc_enrichment_summary.json")

    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="§7.4: IMPC Mouse Phenotype Enrichment")
    parser.add_argument(
        "--novel-candidates", default="results/novel_candidates.csv",
        help="Path to novel candidates CSV (default: results/novel_candidates.csv)"
    )
    parser.add_argument(
        "--gene-col", default="gene_symbol",
        help="Gene symbol column name (default: gene_symbol)"
    )
    args = parser.parse_args()
    run_enrichment_test(args.novel_candidates, args.gene_col)
