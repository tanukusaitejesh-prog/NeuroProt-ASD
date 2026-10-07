"""Metrics with bootstrap confidence intervals."""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve


def sens_at_spec(y, s, spec=0.95):
    fpr, tpr, _ = roc_curve(y, s)
    ok = fpr <= 1 - spec
    return float(tpr[ok].max()) if ok.any() else 0.0


def metrics(y, s):
    y, s = np.asarray(y), np.asarray(s, dtype=float)
    keep = ~np.isnan(s)
    y, s = y[keep], s[keep]
    if len(np.unique(y)) < 2:
        return dict(n=len(y), n_pos=int(y.sum()), auroc=np.nan, auprc=np.nan, sens_at_95spec=np.nan)
    return dict(n=len(y), n_pos=int(y.sum()), auroc=roc_auc_score(y, s),
                auprc=average_precision_score(y, s), sens_at_95spec=sens_at_spec(y, s))


def bootstrap_auroc(y, s, n=1000, seed=0):
    rng = np.random.default_rng(seed)
    y, s = np.asarray(y), np.asarray(s, dtype=float)
    keep = ~np.isnan(s)
    y, s = y[keep], s[keep]
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    vals = []
    for _ in range(n):
        idx = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        vals.append(roc_auc_score(y[idx], s[idx]))
    return np.percentile(vals, [2.5, 97.5])


def paired_bootstrap_delta(y, s1, s2, n=2000, seed=0):
    """AUROC(s1) - AUROC(s2) on the same variants; returns (delta, ci_lo, ci_hi, p_two_sided)."""
    rng = np.random.default_rng(seed)
    y, s1, s2 = map(np.asarray, (y, s1, s2))
    keep = ~(np.isnan(s1.astype(float)) | np.isnan(s2.astype(float)))
    y, s1, s2 = y[keep], s1[keep], s2[keep]
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    d = []
    for _ in range(n):
        idx = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        d.append(roc_auc_score(y[idx], s1[idx]) - roc_auc_score(y[idx], s2[idx]))
    d = np.array(d)
    obs = roc_auc_score(y, s1) - roc_auc_score(y, s2)
    p = 2 * min((d <= 0).mean(), (d >= 0).mean())
    return obs, *np.percentile(d, [2.5, 97.5]), min(p, 1.0)


def score_table(df, label_col, score_cols, group_col=None):
    """Metrics per score (and per group if given). Scores must be oriented: higher = more damaging."""
    rows = []
    groups = [("all", df)] + (list(df.groupby(group_col)) if group_col else [])
    for gname, d in groups:
        for c in score_cols:
            if c not in d:
                continue
            m = metrics(d[label_col], d[c])
            rows.append(dict(group=gname, score=c, **m))
    return pd.DataFrame(rows)


def spearman_vs_sge(df, score_cols, sge_col="sge_score"):
    """SGE: lower score = more depleted = more damaging, so report rho against -sge."""
    rows = []
    for c in score_cols:
        d = df[[c, sge_col]].dropna()
        if len(d) > 5:
            rho, p = spearmanr(d[c], -d[sge_col])
            rows.append(dict(score=c, n=len(d), spearman=rho, p=p))
    return pd.DataFrame(rows)
