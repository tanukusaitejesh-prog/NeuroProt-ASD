"""
Model Training, Six-Number Confound Ablation (§5.2), and Multi-Disorder Classifier (§5.3)
"""

from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.multioutput import MultiOutputClassifier
from sklearn.preprocessing import StandardScaler

from config import RANDOM_SEED, TARGET_DISORDERS


def evaluate_oof_auroc(
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    n_splits: int = 5,
    C: float = 1.0
) -> float:
    """Computes out-of-fold AUROC using paralog-safe GroupKFold cross-validation."""
    gkf = GroupKFold(n_splits=n_splits)
    oof_preds = np.zeros(len(y))

    # Identify unique groups
    for train_idx, val_idx in gkf.split(X, y, groups=groups):
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_idx])
        X_val = scaler.transform(X[val_idx])

        clf = LogisticRegression(
            C=C,
            class_weight="balanced",
            max_iter=2000,
            random_state=RANDOM_SEED
        )
        clf.fit(X_train, y[train_idx])
        oof_preds[val_idx] = clf.predict_proba(X_val)[:, 1]

    return float(roc_auc_score(y, oof_preds))


def run_six_number_comparison(
    X_confound: np.ndarray,
    X_conservation: np.ndarray,
    X_esm2: np.ndarray,
    X_saprot: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    n_splits: int = 5
) -> pd.DataFrame:
    """Executes Revised §5.2 — Six Numbers, Not Four.

    1. Simple confound-only AUROC
    2. Conservation-augmented confound-only AUROC
    3. ESM2 embedding-only AUROC
    4. SaProt (structure-aware) embedding-only AUROC
    5. ESM2 + full confound + conservation combined AUROC
    6. ESM2 + SaProt + full confound + conservation, all combined AUROC
    """
    print("\n" + "=" * 65)
    print("REVISED §5.2 — SIX-NUMBER CONFOUND & EMBEDDING COMPARISON")
    print("=" * 65)

    results = []

    # 1. Confound-only
    auc1 = evaluate_oof_auroc(X_confound, y, groups, n_splits=n_splits)
    results.append({"Model": "1. Simple confound-only", "AUROC": auc1, "Feature Space": "Confounders (length, GC, family, mut_rate)"})

    # 2. Conservation-augmented confound-only
    X_conf_cons = np.hstack([X_confound, X_conservation])
    auc2 = evaluate_oof_auroc(X_conf_cons, y, groups, n_splits=n_splits)
    results.append({"Model": "2. Conservation-augmented confound-only", "AUROC": auc2, "Feature Space": "Confounders + 1D Conservation"})

    # 3. ESM2 embedding-only
    auc3 = evaluate_oof_auroc(X_esm2, y, groups, n_splits=n_splits)
    results.append({"Model": "3. ESM2 embedding-only", "AUROC": auc3, "Feature Space": "ESM-2 sequence embeddings"})

    # 4. SaProt structure-aware embedding-only
    auc4 = evaluate_oof_auroc(X_saprot, y, groups, n_splits=n_splits)
    results.append({"Model": "4. SaProt structure-aware embedding-only", "AUROC": auc4, "Feature Space": "SaProt 3D structure-aware embeddings"})

    # 5. ESM2 + confound + conservation
    X_comb5 = np.hstack([X_esm2, X_confound, X_conservation])
    auc5 = evaluate_oof_auroc(X_comb5, y, groups, n_splits=n_splits)
    results.append({"Model": "5. ESM2 + confound + conservation combined", "AUROC": auc5, "Feature Space": "ESM2 + Confounders + Conservation"})

    # 6. ESM2 + SaProt + confound + conservation, all combined
    X_comb6 = np.hstack([X_esm2, X_saprot, X_confound, X_conservation])
    auc6 = evaluate_oof_auroc(X_comb6, y, groups, n_splits=n_splits)
    results.append({"Model": "6. All Combined (ESM2 + SaProt + Confounders + Conservation)", "AUROC": auc6, "Feature Space": "Full sequence, structure, and genomic signal"})

    df_results = pd.DataFrame(results)
    for _, row in df_results.iterrows():
        print(f"  {row['Model']:<55} AUROC = {row['AUROC']:.4f}")
    print("=" * 65)

    return df_results


def run_joint_multidisorder_model(
    X: np.ndarray,
    Y_multidisorder: np.ndarray,
    disorder_names: List[str] = TARGET_DISORDERS
) -> Tuple[MultiOutputClassifier, Dict[str, float]]:
    """Executes New §5.3 — Joint Multi-Disorder Modeling.

    Trains a multi-output classifier sharing the embedding representation.
    Computes cosine similarity between ASD's learned weight vector and each
    other disorder's weight vector (ID, SCZ, EPI).
    """
    print("\n" + "=" * 65)
    print("NEW §5.3 — JOINT MULTI-DISORDER MODELING (COSINE SIMILARITY)")
    print("=" * 65)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    clf = MultiOutputClassifier(
        LogisticRegression(class_weight="balanced", max_iter=3000, random_state=RANDOM_SEED)
    )
    clf.fit(X_scaled, Y_multidisorder)

    # First estimator corresponds to ASD
    asd_coef = clf.estimators_[0].coef_.ravel()
    asd_norm = np.linalg.norm(asd_coef)

    similarities = {}
    print(f"Learned disorder weight vector alignment vs ASD:")
    for i, disorder in enumerate(disorder_names[1:], start=1):
        other_coef = clf.estimators_[i].coef_.ravel()
        other_norm = np.linalg.norm(other_coef)
        if asd_norm > 0 and other_norm > 0:
            cos_sim = float(np.dot(asd_coef, other_coef) / (asd_norm * other_norm))
        else:
            cos_sim = 0.0
        similarities[disorder] = cos_sim
        print(f"  ASD vs {disorder:<6} Weight Cosine Similarity = {cos_sim:+.4f}")

    print("=" * 65)
    return clf, similarities
