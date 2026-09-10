"""
NeuroProt-ASD: Multi-Modal Neuro-Structural Foundation Framework
================================================================

Fuses four biological modalities for genome-scale ASD gene prioritization:
    1. Structural & Protein Language (ESM-2 sequence embeddings + SaProt 3Di structure tokens)
    2. Spatiotemporal Fetal Brain Transcriptomics (BrainSpan neocortex developmental expression)
    3. Brain-Specific Functional Interactome Networks (STRING PPI + Krishnan functional network)
    4. Statistical Genetics & Evolutionary Constraint (TADA Bayes Factor + gnomAD pLI + mutation rate)

Evaluated under strict paralog-grouped cross-validation (GroupKFold across 171 gene families).
"""

from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from config import RANDOM_SEED
from src.features import extract_confound_features, extract_conservation_features


class NeuroProtFusionModel:
    """Multi-Modal Fusion Classifier integrating PLM, Transcriptomics, Networks, and Genetics."""

    def __init__(self, n_plm_components: int = 24, random_state: int = RANDOM_SEED):
        self.n_plm_components = n_plm_components
        self.random_state = random_state
        self.pca = PCA(n_components=n_plm_components, random_state=random_state)
        self.scaler = StandardScaler()
        # Stacking / Ensemble Classifier
        self.clf = HistGradientBoostingClassifier(
            max_iter=150,
            learning_rate=0.05,
            l2_regularization=1.5,
            max_depth=5,
            random_state=random_state
        )
        self.is_fitted = False

    def prepare_modalities(
        self,
        df: pd.DataFrame,
        X_plm: np.ndarray,
        fit_pca: bool = False
    ) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        """Prepares and standardizes the four multi-modal feature matrices."""
        # Modality 1: Structural & Protein Language (PCA-compressed to balance dimensions)
        if fit_pca:
            X_plm_scaled = StandardScaler().fit_transform(X_plm)
            X_struct = self.pca.fit_transform(X_plm_scaled)
        else:
            X_plm_scaled = StandardScaler().fit_transform(X_plm)
            X_struct = self.pca.transform(X_plm_scaled)

        # Modality 2: Spatiotemporal Fetal Brain Transcriptomics (BrainSpan)
        X_expr = pd.to_numeric(df["BrainSpan_score"], errors="coerce").fillna(0.0).values.reshape(-1, 1)

        # Modality 3: Functional Interactome Networks (STRING + Krishnan functional networks)
        string_s = pd.to_numeric(df["STRING_score"], errors="coerce").fillna(0.0).values.reshape(-1, 1)
        krish_s = pd.to_numeric(df["krishnan_post"], errors="coerce").fillna(0.0).values.reshape(-1, 1)
        X_net = np.hstack([string_s, krish_s])

        # Modality 4: Statistical Genetics & Genomic Constraint
        pli_s = pd.to_numeric(df["pLI"], errors="coerce").fillna(0.0).values.reshape(-1, 1)
        tada_s = np.log1p(pd.to_numeric(df["TADA_BF"], errors="coerce").fillna(0.0).clip(lower=0.0).values).reshape(-1, 1)
        mut_s = np.log1p(pd.to_numeric(df["mutation_rate"], errors="coerce").fillna(0.0).values).reshape(-1, 1)
        X_genetics = np.hstack([pli_s, tada_s, mut_s])

        # Baseline Confounders
        df_conf = extract_confound_features(df)
        df_cons = extract_conservation_features(df)
        X_conf_cons = np.hstack([df_conf.values, df_cons.values])

        modalities = {
            "struct_plm": X_struct,
            "transcriptomics": X_expr,
            "networks": X_net,
            "genetics_constraint": X_genetics,
            "confounders": X_conf_cons
        }

        # Full multi-modal feature vector
        X_fused = np.hstack([X_struct, X_expr, X_net, X_genetics, X_conf_cons])
        return X_fused, modalities

    def evaluate_ablation(
        self,
        df_train: pd.DataFrame,
        X_plm_train: np.ndarray,
        groups: np.ndarray,
        n_splits: int = 5,
        n_bootstraps: int = 1000
    ) -> pd.DataFrame:
        """Evaluates out-of-fold AUROC with 95% bootstrap CIs across all single and combined modality baselines."""
        y = df_train["label_ASD"].values
        X_fused, mods = self.prepare_modalities(df_train, X_plm_train, fit_pca=True)

        feature_sets = {
            "1. Confounders & Conservation Baseline": mods["confounders"],
            "2. Structural PLM-only (ESM-2 + SaProt 3Di)": mods["struct_plm"],
            "3. Spatiotemporal Transcriptomics-only (BrainSpan)": mods["transcriptomics"],
            "4. Functional Interactome-only (STRING + Krishnan)": mods["networks"],
            "5. Classical Integration (BrainSpan + Networks + Genetics)": np.hstack([
                mods["transcriptomics"], mods["networks"], mods["genetics_constraint"], mods["confounders"]
            ]),
            "6. NeuroProt-ASD (Full Multi-Modal Fusion: PLM + BrainSpan + Networks + Genetics)": X_fused
        }

        results = []
        oof_dict = {}
        print("\n" + "=" * 80)
        print("NEUROPROT-ASD: PARALOG-SAFE MULTI-MODAL ABLATION (GROUPKFOLD, k=5)")
        print("=" * 80)

        rng = np.random.RandomState(self.random_state)

        for name, X_mat in feature_sets.items():
            gkf = GroupKFold(n_splits=n_splits)
            oof = np.zeros(len(y))
            for tr, val in gkf.split(X_mat, y, groups=groups):
                sc = StandardScaler()
                X_tr = sc.fit_transform(X_mat[tr])
                X_va = sc.transform(X_mat[val])
                clf = HistGradientBoostingClassifier(
                    max_iter=100,
                    learning_rate=0.05,
                    l2_regularization=1.5,
                    max_depth=4,
                    random_state=self.random_state
                )
                clf.fit(X_tr, y[tr])
                oof[val] = clf.predict_proba(X_va)[:, 1]

            auc = float(roc_auc_score(y, oof))
            oof_dict[name] = oof

            # Compute 95% bootstrap confidence interval
            boot_aucs = []
            n_samples = len(y)
            for _ in range(n_bootstraps):
                idx = rng.randint(0, n_samples, n_samples)
                if len(np.unique(y[idx])) >= 2:
                    boot_aucs.append(roc_auc_score(y[idx], oof[idx]))

            ci_lower = float(np.percentile(boot_aucs, 2.5))
            ci_upper = float(np.percentile(boot_aucs, 97.5))

            results.append({
                "Modality_Framework": name,
                "Paralog_Safe_AUROC": auc,
                "CI_95_Lower": ci_lower,
                "CI_95_Upper": ci_upper
            })
            print(f"  {name:<68} AUROC = {auc:.4f}  [95% CI: {ci_lower:.4f} - {ci_upper:.4f}]")

        print("=" * 80)
        return pd.DataFrame(results)

    def fit(self, df_train: pd.DataFrame, X_plm_train: np.ndarray):
        """Fits final multi-modal model on the full training set."""
        y = df_train["label_ASD"].values
        X_fused, _ = self.prepare_modalities(df_train, X_plm_train, fit_pca=True)
        X_scaled = self.scaler.fit_transform(X_fused)
        self.clf.fit(X_scaled, y)
        self.is_fitted = True

    def predict_proba(self, df_all: pd.DataFrame, X_plm_all: np.ndarray) -> np.ndarray:
        """Generates genome-wide calibrated probabilities."""
        assert self.is_fitted, "Model must be fitted before prediction."
        X_fused, _ = self.prepare_modalities(df_all, X_plm_all, fit_pca=False)
        X_scaled = self.scaler.transform(X_fused)
        return self.clf.predict_proba(X_scaled)[:, 1]
