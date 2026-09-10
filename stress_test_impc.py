"""
Stress-Testing the 393-Gene IMPC Enrichment & Candidate Selection Pipeline
==========================================================================

Executes 4 comprehensive sensitivity tests for hidden ascertainment/selection bias:

1. THE CONFOUNDER COUNTERFACTUAL:
   Ranks genes using the Confounder-Only model (AUROC = 0.8584), extracts its
   novel candidates vs forecASD, and tests if it achieves comparable or weaker
   IMPC neurobehavioral enrichment than the PLM combined model.

2. PHENOTYPIC SPECIFICITY (NEGATIVE CONTROLS):
   Tests if the 393 candidates are specifically enriched for Neuro/Behavior
   (MP:0005386, MP:0003631) versus non-neurological control organ systems:
   - Cardiovascular (MP:0005385)
   - Immune System (MP:0005387)
   - Integument/Skin (MP:0010771)

3. MULTIVARIABLE BIAS ADJUSTMENT:
   Fits a logistic regression on all IMPC-tested genes:
   NeuroHit ~ IsCandidate + log(Length) + pLI + MutationRate
   to test if candidate status remains significant after controlling for confounders.

4. TOP-N CUTOFF SENSITIVITY:
   Evaluates Odds Ratio and Fisher p-value across N in [250, 400, 500, 750, 1000].
"""

import json
import os
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy.stats import fisher_exact, norm
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

# Safe stdout on Windows
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from config import BASE_DIR, FIGURES_DIR, RANDOM_SEED, RESULTS_DIR
from src.data import extract_gene_family, get_representative_sequence, load_curated_dataset
from src.features import extract_confound_features, extract_conservation_features

IMPC_SOLR_GENO = "https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select"
IMPC_SOLR_STAT = "https://www.ebi.ac.uk/mi/impc/solr/statistical-result/select"
MGI_LOCAL = BASE_DIR / "HOM_MouseHumanSequence.rpt"


def get_mouse_to_human():
    df = pd.read_csv(MGI_LOCAL, sep="\t", low_memory=False)
    human_df = df[df["Common Organism Name"] == "human"][["DB Class Key", "Symbol"]].rename(columns={"Symbol": "human_symbol"})
    mouse_df = df[df["Common Organism Name"] == "mouse, laboratory"][["DB Class Key", "Symbol"]].rename(columns={"Symbol": "mouse_symbol"})
    ortho = mouse_df.merge(human_df, on="DB Class Key")
    return dict(zip(ortho["mouse_symbol"].str.upper(), ortho["human_symbol"]))


def query_impc_facet(query_str: str) -> set:
    r = requests.get(
        IMPC_SOLR_GENO,
        params={
            "q": query_str,
            "rows": 0,
            "facet": "true",
            "facet.field": "marker_symbol",
            "facet.limit": -1,
            "facet.mincount": 1,
            "wt": "json",
        },
        timeout=60,
    )
    r.raise_for_status()
    f = r.json()["facet_counts"]["facet_fields"]["marker_symbol"]
    return {f[i].strip() for i in range(0, len(f), 2)}


def fetch_tested_impc_genes() -> set:
    r = requests.get(
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
    r.raise_for_status()
    f = r.json()["facet_counts"]["facet_fields"]["marker_symbol"]
    return {f[i].strip() for i in range(0, len(f), 2)}


def run_fisher_test(candidate_genes: set, target_human_genes: set, background_human_genes: set):
    cand_tested = candidate_genes & background_human_genes
    cand_hit = cand_tested & target_human_genes
    bg_tested = background_human_genes - candidate_genes
    bg_hit = bg_tested & target_human_genes

    a = len(cand_hit)
    b = len(cand_tested) - a
    c = len(bg_hit)
    d = len(bg_tested) - c

    table = np.array([[a, b], [c, d]])
    or_val, p_val = fisher_exact(table, alternative="greater")
    cand_rate = a / len(cand_tested) if len(cand_tested) > 0 else 0.0
    bg_rate = c / len(bg_tested) if len(bg_tested) > 0 else 0.0

    return {
        "n_tested": len(cand_tested),
        "n_hits": a,
        "cand_rate": cand_rate,
        "bg_rate": bg_rate,
        "odds_ratio": float(or_val),
        "p_value": float(p_val),
        "table": [[a, b], [c, d]]
    }


def main():
    print("=" * 70)
    print("  STRESS-TESTING IMPC NOVEL CANDIDATE ENRICHMENT FOR SELECTION BIAS")
    print("=" * 70)

    # 1. Load Base Data & Orthologs
    m2h = get_mouse_to_human()
    all_tested_mouse = fetch_tested_impc_genes()
    all_tested_human = {m2h[m.upper()] for m in all_tested_mouse if m.upper() in m2h}
    print(f"[INFO] IMPC tested universe: {len(all_tested_human)} human orthologs")

    # Load our calibrated genome scores and forecASD comparison
    calib_df = pd.read_csv(RESULTS_DIR / "calibrated_genome_wide_scores.csv")
    forecasd_df = pd.read_csv(RESULTS_DIR / "forecasd_comparison.csv")
    novel_df = pd.read_csv(RESULTS_DIR / "novel_candidates.csv")
    our_novel_genes = set(novel_df["gene_symbol"])
    print(f"[INFO] Our model uniquely nominated candidates: {len(our_novel_genes)}")

    # ── STRESS TEST 1: The Confounder-Only Model Counterfactual ────────
    print("\n" + "-" * 70)
    print("STRESS TEST 1: Confounder-Only Model Counterfactual (AUROC Paradox Test)")
    print("-" * 70)
    # Train confounder-only model on curated dataset
    df_curated = load_curated_dataset(balanced_negatives=True)
    conf_train = extract_confound_features(df_curated).values
    y_train = df_curated["label_ASD"].values

    scaler = StandardScaler()
    conf_train_scaled = scaler.fit_transform(conf_train)
    conf_clf = LogisticRegression(class_weight="balanced", max_iter=2000, random_state=RANDOM_SEED)
    conf_clf.fit(conf_train_scaled, y_train)

    # Score full genome using confounder-only model
    df_genome = load_curated_dataset(balanced_negatives=False)
    conf_genome = extract_confound_features(df_genome).values
    conf_genome_scaled = scaler.transform(conf_genome)
    df_genome["confound_score"] = conf_clf.predict_proba(conf_genome_scaled)[:, 1]

    # Find top 500 from confounder-only model
    fc_top500 = set(forecasd_df.nlargest(500, "forecasd_score")["gene_symbol"])
    known_sfari = set(df_genome.loc[df_genome["label_ASD"] == 1, "gene_symbol"])
    conf_top500 = set(df_genome.nlargest(500, "confound_score")["gene_symbol"])

    confound_novel_genes = conf_top500 - fc_top500 - known_sfari
    print(f"Confounder-only model uniquely nominated novel candidates: {len(confound_novel_genes)}")

    # Query Neuro/Behavior
    neuro_mouse = query_impc_facet('top_level_mp_term_id:"MP:0005386" OR top_level_mp_term_id:"MP:0003631"')
    neuro_human = {m2h[m.upper()] for m in neuro_mouse if m.upper() in m2h}

    # Test our model vs confounder-only model
    res_our_model = run_fisher_test(our_novel_genes, neuro_human, all_tested_human)
    res_conf_model = run_fisher_test(confound_novel_genes, neuro_human, all_tested_human)

    print("\nComparison of Novel Candidates in IMPC In-Vivo Neuro Knockouts:")
    print(f"  PLM Combined Model (AUROC 0.702):")
    print(f"    Hit Rate:    {res_our_model['cand_rate']*100:.2f}% ({res_our_model['n_hits']}/{res_our_model['n_tested']})")
    print(f"    Odds Ratio:  {res_our_model['odds_ratio']:.3f}")
    print(f"    Fisher p:    {res_our_model['p_value']:.4e}  {'*** SIGNIFICANT' if res_our_model['p_value'] < 0.05 else ''}")
    print(f"  Confounder-Only Model (AUROC 0.858):")
    print(f"    Hit Rate:    {res_conf_model['cand_rate']*100:.2f}% ({res_conf_model['n_hits']}/{res_conf_model['n_tested']})")
    print(f"    Odds Ratio:  {res_conf_model['odds_ratio']:.3f}")
    print(f"    Fisher p:    {res_conf_model['p_value']:.4e}  {'*** SIGNIFICANT' if res_conf_model['p_value'] < 0.05 else ''}")

    # ── STRESS TEST 2: Phenotypic Specificity (Negative Controls) ──────
    print("\n" + "-" * 70)
    print("STRESS TEST 2: Phenotypic Specificity Across Physiological Systems")
    print("-" * 70)
    systems = {
        "Nervous System / Behavior (Target)": 'top_level_mp_term_id:"MP:0005386" OR top_level_mp_term_id:"MP:0003631"',
        "Cardiovascular System (Control)": 'top_level_mp_term_id:"MP:0005385"',
        "Immune System (Control)": 'top_level_mp_term_id:"MP:0005387"',
        "Integument / Skin (Control)": 'top_level_mp_term_id:"MP:0010771"',
        "Auditory / Ear (Control)": 'top_level_mp_term_id:"MP:0005377"',
    }

    spec_results = {}
    for sys_name, q in systems.items():
        m_set = query_impc_facet(q)
        h_set = {m2h[m.upper()] for m in m_set if m.upper() in m2h}
        res = run_fisher_test(our_novel_genes, h_set, all_tested_human)
        spec_results[sys_name] = res
        print(f"  {sys_name:<38} Hit Rate: {res['cand_rate']*100:.1f}% (vs bg {res['bg_rate']*100:.1f}%) | OR: {res['odds_ratio']:.2f} | p = {res['p_value']:.4f}")

    # ── STRESS TEST 3: Multivariable Confounder Adjustment ────────────
    print("\n" + "-" * 70)
    print("STRESS TEST 3: Multivariable Logistic Regression (Adjusting for Length & pLI)")
    print("-" * 70)
    # Assemble dataframe of all tested IMPC genes
    df_impc_tested = df_genome[df_genome["gene_symbol"].isin(all_tested_human)].copy()
    df_impc_tested["is_candidate"] = df_impc_tested["gene_symbol"].isin(our_novel_genes).astype(int)
    df_impc_tested["has_neuro_phenotype"] = df_impc_tested["gene_symbol"].isin(neuro_human).astype(int)

    # Compute length
    lens = [len(get_representative_sequence(g)) for g in df_impc_tested["gene_symbol"]]
    df_impc_tested["log_length"] = np.log1p(lens)
    df_impc_tested["pLI"] = pd.to_numeric(df_impc_tested["pLI"], errors="coerce").fillna(0.0)
    df_impc_tested["mutation_rate"] = pd.to_numeric(df_impc_tested["mutation_rate"], errors="coerce").fillna(0.0)

    # Fit Logistic Regression and compute Wald test statistics via Hessian
    X_mat = df_impc_tested[["is_candidate", "log_length", "pLI", "mutation_rate"]].values
    X_design = np.hstack([np.ones((len(X_mat), 1)), X_mat])  # add intercept
    y_pheno = df_impc_tested["has_neuro_phenotype"].values

    # Unregularized / weak regularization fit to obtain MLE
    logit_clf = LogisticRegression(penalty=None, max_iter=2000, random_state=RANDOM_SEED)
    logit_clf.fit(X_mat, y_pheno)

    betas = np.array([logit_clf.intercept_[0]] + list(logit_clf.coef_[0]))
    p_pred = logit_clf.predict_proba(X_mat)[:, 1]
    
    # Compute information matrix W and Hessian X^T W X
    w = p_pred * (1.0 - p_pred)
    H = X_design.T @ (X_design * w[:, None])
    cov_beta = np.linalg.inv(H)
    se_beta = np.sqrt(np.diag(cov_beta))
    z_scores = betas / se_beta
    p_values = 2.0 * (1.0 - norm.cdf(np.abs(z_scores)))

    var_names = ["Intercept", "is_candidate", "log_length", "pLI", "mutation_rate"]
    print(f"{'Variable':<15} {'Coef':>10} {'Std.Err':>10} {'z':>10} {'P>|z|':>10} {'OR':>10}")
    print("-" * 70)
    for name, b, se, z, p in zip(var_names, betas, se_beta, z_scores, p_values):
        print(f"{name:<15} {b:>10.4f} {se:>10.4f} {z:>10.4f} {p:>10.4e} {np.exp(b):>10.3f}")

    cand_coef = betas[1]
    cand_p = p_values[1]
    cand_or = np.exp(cand_coef)
    print(f"\nAdjusted Odds Ratio for Candidate Status: OR = {cand_or:.3f} (p = {cand_p:.4e})")

    # ── STRESS TEST 4: Top-N Cutoff Sensitivity ───────────────────────
    print("\n" + "-" * 70)
    print("STRESS TEST 4: Sensitivity to Top-N Candidate Threshold")
    print("-" * 70)
    cutoffs = [250, 400, 500, 750, 1000]
    cutoff_res = []

    for N in cutoffs:
        top_ours = set(calib_df.nlargest(N, "our_score")["gene_symbol"])
        top_fc = set(forecasd_df.nlargest(N, "forecasd_score")["gene_symbol"])
        novel_set = top_ours - top_fc - known_sfari
        res_n = run_fisher_test(novel_set, neuro_human, all_tested_human)
        cutoff_res.append({
            "N": N,
            "n_novel": len(novel_set),
            "n_tested": res_n["n_tested"],
            "n_hits": res_n["n_hits"],
            "odds_ratio": res_n["odds_ratio"],
            "p_value": res_n["p_value"]
        })
        print(f"  Top-{N:<4}: Novel={len(novel_set):<4} Tested={res_n['n_tested']:<4} Hits={res_n['n_hits']:<3} OR={res_n['odds_ratio']:.3f} p={res_n['p_value']:.4e}")

    # Plot sensitivity curve
    fig, ax1 = plt.subplots(figsize=(7, 4.5))
    df_cut = pd.DataFrame(cutoff_res)
    color = "#1f77b4"
    ax1.set_xlabel("Top-N Candidate Threshold", fontsize=11)
    ax1.set_ylabel("Odds Ratio (OR)", color=color, fontsize=11)
    ax1.plot(df_cut["N"], df_cut["odds_ratio"], marker="o", color=color, lw=2, label="Odds Ratio")
    ax1.tick_params(axis="y", labelcolor=color)
    ax1.axhline(1.0, color="grey", ls="--", alpha=0.5)

    ax2 = ax1.twinx()
    color = "crimson"
    ax2.set_ylabel(r"$-\log_{10}(P_{\mathrm{Fisher}})$", color=color, fontsize=11)
    ax2.plot(df_cut["N"], -np.log10(df_cut["p_value"]), marker="s", color=color, lw=2, label=r"$-\log_{10}(P)$")
    ax2.axhline(-np.log10(0.05), color="crimson", ls=":", label="p = 0.05 cutoff")
    ax2.tick_params(axis="y", labelcolor=color)

    fig.tight_layout()
    plot_path = FIGURES_DIR / "impc_sensitivity_thresholds.pdf"
    fig.savefig(plot_path, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "impc_sensitivity_thresholds.png", dpi=150, bbox_inches="tight")
    print(f"\n[INFO] Saved threshold sensitivity plot to {plot_path}")

    # Save summary JSON
    stress_summary = {
        "confounder_counterfactual": {
            "plm_combined_or": res_our_model["odds_ratio"],
            "plm_combined_p": res_our_model["p_value"],
            "confound_only_or": res_conf_model["odds_ratio"],
            "confound_only_p": res_conf_model["p_value"],
        },
        "phenotypic_specificity": {k: {"or": v["odds_ratio"], "p": v["p_value"]} for k, v in spec_results.items()},
        "multivariable_regression": {
            "candidate_adjusted_or": float(cand_or),
            "candidate_adjusted_p": float(cand_p),
        },
        "threshold_sensitivity": cutoff_res,
    }
    with open(RESULTS_DIR / "impc_stress_test_summary.json", "w") as f:
        json.dump(stress_summary, f, indent=2)
    print(f"[INFO] Saved results/impc_stress_test_summary.json")

    print("\n" + "=" * 70)
    print("  STRESS TESTS COMPLETED")
    print("=" * 70)


if __name__ == "__main__":
    main()
