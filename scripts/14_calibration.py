"""Clinical calibration of held-out snRNA-VEP scores into ACMG/ClinGen PP3/BP4 evidence strengths.

Uses only held-out predictions (each variant scored by a model that never saw its gene). For score thresholds t,
positive likelihood ratio LR+ = P(score >= t | pathogenic) / P(score >= t | control) and LR- = P(score < t | path) /
P(score < t | control). Bayesian ACMG point system (Tavtigian 2018/2020): OddsPath for supporting/moderate/strong
pathogenic = 2.08 / 4.33 / 18.7; benign supporting/moderate = 1/2.08 (0.48) / 1/4.33 (0.23). A threshold is reported at a
strength only if the one-sided 95% bootstrap lower bound of LR+ (upper bound of LR- for benign) meets it.
Calibration is per mechanism context (scores are not comparable across contexts).
Usage: python scripts/14_calibration.py [predictions.tsv.gz] [score column] [context]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

f = sys.argv[1] if len(sys.argv) > 1 else str(config.RESULTS / "v2_heldout_predictions.tsv.gz")
col = sys.argv[2] if len(sys.argv) > 2 else "snrnavep_v2"
P = pd.read_csv(f, sep="\t", low_memory=False)
ctx = sys.argv[3] if len(sys.argv) > 3 else "dominant"
P = P[P[col].notna() & (P.ctx == ctx)]
# a control appears once per context; for calibration each (variant, context) row is one assessment
y, s = P.y_path.values.astype(int), P[col].values
rng = np.random.default_rng(0)
PATH = {"supporting": 2.08, "moderate": 4.33, "strong": 18.7}
BEN = {"supporting": 1 / 2.08, "moderate": 1 / 4.33}


def lrs(y, s, t):
    tp, fp = ((s >= t) & (y == 1)).sum(), ((s >= t) & (y == 0)).sum()
    fn, tn = ((s < t) & (y == 1)).sum(), ((s < t) & (y == 0)).sum()
    npos, nneg = max(1, y.sum()), max(1, (1 - y).sum())
    lrp = (tp / npos) / max(fp / nneg, 0.5 / nneg)          # continuity: at least half a false positive
    lrn = max(fn / npos, 0.5 / npos) / max(tn / nneg, 1e-9)
    return lrp, lrn


grid = np.unique(np.quantile(s, np.linspace(0.02, 0.98, 97)))
boot = []
for _ in range(1000):
    i = np.concatenate([rng.choice(np.where(y == 1)[0], y.sum()), rng.choice(np.where(y == 0)[0], (1 - y).sum())])
    boot.append([lrs(y[i], s[i], t) for t in grid])
boot = np.array(boot)                                         # (B, T, 2)
point = np.array([lrs(y, s, t) for t in grid])
lo_p = np.quantile(boot[:, :, 0], 0.05, 0)
hi_n = np.quantile(boot[:, :, 1], 0.95, 0)
rows = []
for name, thr in PATH.items():
    ok = np.where(lo_p >= thr)[0]
    if len(ok):
        k = ok[0]
        rows.append(dict(evidence=f"PP3_{name}", threshold=f">= {grid[k]:.3f}", LR=point[k, 0], bound=lo_p[k],
                         frac_path=(s[y == 1] >= grid[k]).mean(), frac_ctrl=(s[y == 0] >= grid[k]).mean()))
    else:
        rows.append(dict(evidence=f"PP3_{name}", threshold="not reached"))
for name, thr in BEN.items():
    ok = np.where(hi_n <= thr)[0]
    if len(ok):
        k = ok[-1]
        rows.append(dict(evidence=f"BP4_{name}", threshold=f"< {grid[k]:.3f}", LR=point[k, 1], bound=hi_n[k],
                         frac_path=(s[y == 1] < grid[k]).mean(), frac_ctrl=(s[y == 0] < grid[k]).mean()))
    else:
        rows.append(dict(evidence=f"BP4_{name}", threshold="not reached"))
out = pd.DataFrame(rows)
print(f"{col} [{ctx}]: {len(P)} held-out assessments ({y.sum()} pathogenic)")
print(out.round(3).to_string(index=False))
out.to_csv(config.RESULTS / f"calibration_{col}_{ctx}.tsv", sep="\t", index=False)
