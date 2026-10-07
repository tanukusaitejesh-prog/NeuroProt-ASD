"""snRNA-VEP model: gene-agnostic features -> pathogenicity, evaluated leave-one-gene-out / time-split.

Feature groups can be switched off for ablations. The `population` group (gnomAD depletion) is OFF by
default because the benign labels come from gnomAD; turning it on is only legitimate for SGE or
time-split tests where labels don't come from gnomAD.
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

FEATURE_GROUPS = {
    "position": ["rel_pos", "dist_5p", "dist_3p", "is_snv", "is_ins", "is_del", "transition", "local_gc"],
    "structure2d": ["p_unpaired_wt", "ddG_fold", "ddG_duplex"],
    "contacts": ["n_states_seen", "frac_states_protein_contact", "max_protein_res", "mean_protein_res",
                 "max_protein_chains", "min_protein_dist", "frac_states_snrna_contact",
                 "max_other_rna_contacts"],
    "conservation": ["phylop"],
    "lm": ["rnalm_dll", "evo2_dll"],
    "paralog": ["paralog_patho_density"],
    "population": ["rel_oe_w10"],
}
DEFAULT_GROUPS = ["position", "structure2d", "contacts", "conservation", "lm", "paralog"]


def paralog_density(train, test, window=2):
    """For each test variant: share of training pathogenic variants (other genes) at the same
    family-reference position +/- window. Computed inside each fold, so no label leakage."""
    patho = train[train.label == 1].dropna(subset=["family_ref_pos"])
    out = np.zeros(len(test))
    for i, (fam, pos, gene) in enumerate(zip(test.paralog_family, test.family_ref_pos, test.gene_name)):
        if pd.isna(pos) or pd.isna(fam):
            continue
        p = patho[(patho.paralog_family == fam) & (patho.gene_name != gene)]
        if len(p):
            out[i] = (np.abs(p.family_ref_pos - pos) <= window).mean()
    return out


def make_model(kind="gbm"):
    if kind == "gbm":
        return HistGradientBoostingClassifier(max_depth=3, learning_rate=0.05, max_iter=300,
                                              l2_regularization=1.0, class_weight="balanced",
                                              random_state=0)
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                         LogisticRegression(C=0.3, class_weight="balanced", max_iter=2000))


def _columns(df, groups):
    cols = [c for g in groups for c in FEATURE_GROUPS[g] if c in df.columns]
    return [c for c in cols if df[c].notna().any()]


def fit_predict(train, test, groups=DEFAULT_GROUPS, kind="gbm"):
    train, test = train.copy(), test.copy()
    if "paralog" in groups:
        train["paralog_patho_density"] = 0.0  # within-train leave-gene-out version below
        for g in train.gene_name.unique():
            m = train.gene_name == g
            train.loc[m, "paralog_patho_density"] = paralog_density(train[~m], train[m])
        test["paralog_patho_density"] = paralog_density(train, test)
    cols = _columns(train, groups)
    model = make_model(kind)
    model.fit(train[cols].astype(float), train.label)
    return model.predict_proba(test[cols].astype(float))[:, 1], cols, model


def leave_one_gene_out(data, groups=DEFAULT_GROUPS, kind="gbm", min_pos=3):
    preds = []
    for g in sorted(data.gene_name.unique()):
        test = data[data.gene_name == g]
        train = data[data.gene_name != g]
        if test.label.sum() < min_pos or train.label.sum() < min_pos:
            continue
        p, _, _ = fit_predict(train, test, groups, kind)
        preds.append(test[["gene_name", "key", "label"]].assign(pred=p, fold=g))
    return pd.concat(preds, ignore_index=True) if preds else pd.DataFrame()


def time_split(data, cutoff="2025-01", groups=DEFAULT_GROUPS, kind="gbm"):
    """Train on pathogenic variants first published before `cutoff` (+ all controls of those
    genes); test on pathogenic variants published after, against controls of the same genes."""
    cutoff = pd.Timestamp(cutoff)
    patho = data[data.label == 1]
    old = patho[patho.first_published < cutoff]
    new = patho[patho.first_published >= cutoff]
    ctrl = data[data.label == 0]
    rng = np.random.default_rng(0)
    test_genes = set(new.gene_name)
    ctrl_test = ctrl[ctrl.gene_name.isin(test_genes)]
    ctrl_test = ctrl_test[rng.random(len(ctrl_test)) < 0.5]
    ctrl_train = ctrl.drop(ctrl_test.index)
    train = pd.concat([old, ctrl_train])
    test = pd.concat([new, ctrl_test])
    p, cols, _ = fit_predict(train, test, groups, kind)
    return test[["gene_name", "key", "label"]].assign(pred=p), cols
