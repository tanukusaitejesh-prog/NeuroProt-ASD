"""
Split-Conformal Prediction & Calibrated Candidate Prioritization (§6)
"""

from typing import Tuple
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import StandardScaler

from config import CONFORMAL_CALIBRATION_FRACTION, CONFORMAL_TARGET_COVERAGE, RANDOM_SEED


def fit_conformal_calibrator(
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    test_size: float = CONFORMAL_CALIBRATION_FRACTION,
    target_coverage: float = CONFORMAL_TARGET_COVERAGE
) -> Tuple[LogisticRegression, StandardScaler, float, dict]:
    """Fits model with split-conformal calibration on a disjoint, paralog-grouped fold.

    Returns:
        clf: Fitted classifier
        scaler: Fitted scaler
        qhat: Conformal nonconformity threshold at target coverage
        metrics: Diagnostic calibration metrics
    """
    gss = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=RANDOM_SEED)
    train_idx, calib_idx = next(gss.split(X, y, groups=groups))

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X[train_idx])
    X_calib = scaler.transform(X[calib_idx])

    clf = LogisticRegression(class_weight="balanced", max_iter=3000, random_state=RANDOM_SEED)
    clf.fit(X_train, y[train_idx])

    # Compute nonconformity on calibration set
    calib_probs = clf.predict_proba(X_calib)[:, 1]
    y_calib = y[calib_idx]

    # Nonconformity score: 1 - p for positives, p for negatives
    nonconformity = np.where(y_calib == 1, 1.0 - calib_probs, calib_probs)

    # Conformal quantile with finite-sample correction
    n_calib = len(calib_idx)
    # ceiling((n + 1) * (1 - alpha)) / n
    alpha = 1.0 - target_coverage
    q_level = min(1.0, np.ceil((n_calib + 1) * (1.0 - alpha)) / n_calib)
    qhat = float(np.quantile(nonconformity, q_level))

    # Empirical coverage on calibration set
    covered = (1.0 - calib_probs[y_calib == 1]) <= qhat
    empirical_coverage = float(np.mean(covered)) if len(covered) > 0 else 0.0

    metrics = {
        "n_train": len(train_idx),
        "n_calibration": n_calib,
        "target_coverage": target_coverage,
        "empirical_calibration_coverage": empirical_coverage,
        "qhat_threshold": qhat
    }

    print("\n" + "=" * 65)
    print("§6 SPLIT-CONFORMAL PREDICTION CALIBRATION")
    print("=" * 65)
    print(f"Training genes:        {len(train_idx)}")
    print(f"Calibration genes:     {n_calib} (disjoint, paralog-grouped)")
    print(f"Target Coverage:       {target_coverage * 100:.1f}%")
    print(f"Empirical Coverage:    {empirical_coverage * 100:.1f}%")
    print(f"Nonconformity Cutoff:  q_hat = {qhat:.4f}")
    print("=" * 65)

    return clf, scaler, qhat, metrics


def score_genome_wide(
    clf: LogisticRegression,
    scaler: StandardScaler,
    qhat: float,
    X_all: np.ndarray,
    df_all: pd.DataFrame
) -> pd.DataFrame:
    """Computes calibrated probability scores and confidence tiers for all genes."""
    X_scaled = scaler.transform(X_all)
    probs = clf.predict_proba(X_scaled)[:, 1]

    df_scored = df_all.copy()
    df_scored["our_score"] = probs
    # Included in 90% coverage tier
    df_scored["in_90pct_conformal_set"] = (1.0 - probs) <= qhat
    df_scored["conformal_tier"] = np.where(
        df_scored["in_90pct_conformal_set"],
        "High-Confidence (90% Conformal Tier)",
        "Standard"
    )

    df_scored = df_scored.sort_values("our_score", ascending=False).reset_index(drop=True)
    df_scored["our_rank"] = np.arange(1, len(df_scored) + 1)

    return df_scored
