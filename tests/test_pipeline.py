"""
Test Suite for ASD Gene Prioritization Pipeline
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.conformal import fit_conformal_calibrator, score_genome_wide
from src.data import extract_gene_family, get_representative_sequence, load_curated_dataset
from src.features import (
    SaProtExtractor,
    extract_confound_features,
    extract_conservation_features
)
from src.models import evaluate_oof_auroc, run_joint_multidisorder_model
from src.mutagenesis import compute_mutational_profile, validate_mutagenesis_against_variants


def test_gene_family_extraction():
    """Verify paralog family extraction prevents paralog leakage."""
    assert extract_gene_family("SCN1A") == "SCN"
    assert extract_gene_family("SCN2A") == "SCN"
    assert extract_gene_family("SCN8A") == "SCN"
    assert extract_gene_family("KCNQ2") == "KCNQ"
    assert extract_gene_family("GRIN2B") == "GRIN"


def test_data_loading_and_labels():
    """Verify curated dataset contains valid multi-disorder labels."""
    df = load_curated_dataset(sample_size=40, balanced_negatives=True)
    assert len(df) == 40
    assert "label_ASD" in df.columns
    assert "label_ID" in df.columns
    assert "label_SCZ" in df.columns
    assert "label_EPI" in df.columns
    assert "gene_family" in df.columns
    assert df["label_ASD"].isin([0, 1]).all()


def test_confound_and_conservation_features():
    """Verify confound and conservation feature calculation."""
    df = load_curated_dataset(sample_size=10, balanced_negatives=True)
    conf = extract_confound_features(df)
    cons = extract_conservation_features(df)

    assert conf.shape == (10, 4)
    assert cons.shape == (10, 3)
    assert not conf.isna().any().any()
    assert not cons.isna().any().any()


def test_saprot_tokenization():
    """Verify SaProt 3Di sequence mapping and SA-token formatting."""
    extractor = SaProtExtractor()
    aa_seq = "ACDEFGHIKLMNPQRSTVWY"
    di_seq = extractor.generate_3di_tokens(aa_seq)
    assert len(di_seq) == len(aa_seq)
    sa_tokens = extractor.to_sa_tokens(aa_seq, di_seq)
    assert len(sa_tokens) == len(aa_seq) * 2
    assert sa_tokens[0] == aa_seq[0]
    assert sa_tokens[1] == di_seq[0]


def test_paralog_safe_groupkfold():
    """Verify GroupKFold out-of-fold evaluation."""
    np.random.seed(42)
    N = 50
    X = np.random.randn(N, 5)
    y = np.random.randint(0, 2, size=N)
    # 5 paralog family groups
    groups = np.array([f"FAM_{i % 5}" for i in range(N)])

    auroc = evaluate_oof_auroc(X, y, groups=groups, n_splits=3)
    assert 0.0 <= auroc <= 1.0


def test_conformal_calibration():
    """Verify split-conformal calibration generates valid qhat and coverage."""
    np.random.seed(42)
    N = 60
    X = np.random.randn(N, 6)
    y = np.random.choice([0, 1], size=N, p=[0.7, 0.3])
    groups = np.array([f"GRP_{i % 6}" for i in range(N)])

    clf, scaler, qhat, metrics = fit_conformal_calibrator(
        X, y, groups=groups, test_size=0.25, target_coverage=0.90
    )
    assert 0.0 <= qhat <= 1.0
    assert "empirical_calibration_coverage" in metrics


def test_mutagenesis_profile():
    """Verify in-silico mutagenesis scanner evaluates all positions."""
    test_seq = "MAQSVLVPPGP"

    def mock_predict(seqs):
        return np.array([len(s) * 0.1 + (s.count("P") * 0.2) for s in seqs])

    profile, wt_score = compute_mutational_profile(test_seq, mock_predict)
    assert len(profile) == len(test_seq)
    assert (profile >= 0.0).all()

    # Test variant validation
    stat, p_val, metrics = validate_mutagenesis_against_variants(
        gene_symbol="TEST_GENE",
        profile=profile,
        observed_variant_positions=[3, 7]
    )
    assert 0.0 <= p_val <= 1.0
