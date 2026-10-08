"""Decisive screen for the GRIN idea: does state-resolved structure predict the DIRECTION of functional change
(gain vs loss of function), and does it beat the two things that could explain it trivially?

Labels   : Supplemental Table S4 of Hamed et al. HMG 2025 (ddae156) - 127 GluN1/2A/2B missense variants whose
           GoF/LoF analysis by the Myers et al. (2023) criteria is complete. Direction = "Call for this paper"
           (Increase = GoF, Decrease = LoF). STRICT subset = formal call contains "Likely" (drops Possible /
           Indeterminant), the label-quality sensitivity analysis.
Features : data/processed/nmdar_residue_features.tsv.gz (scripts/31), aggregated over the gating-state series.
Comparison (all under leave-one-GENE-out, so no gene's own variants train its model):
  DOMAIN    one-hot domain from Supplemental Table S9 boundaries  <- the control that could kill the idea
  LIGAND    min distance to agonist and to any ligand             <- the published feature set (static, MCC 0.523)
  STATE     frac states with inter-subunit contact, SD across states, open-minus-non-active delta, mean/max contacts
  LIGAND+STATE
Metric   : AUROC for GoF vs LoF, pooled over held-out genes and per gene, with stratified bootstrap CIs; plus the
           mean per-gene AUROC. A positive result requires STATE to beat DOMAIN and LIGAND with a CI excluding 0.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

rng = np.random.default_rng(0)
GENE_ACC = {"GRIN1": "Q05586", "GRIN2A": "Q12879", "GRIN2B": "Q13224"}

# ---- labels
S4 = pd.read_excel(config.DATA / "external" / "grin_gofl_s4.xlsx", sheet_name="GoF LoF", header=1)
S4 = S4[S4.Gene.isin(GENE_ACC)].copy()
S4["y"] = (S4["Call for this paper"].astype(str).str.strip() == "Increase").astype(int)   # 1 = GoF
S4["strict"] = S4["Formal Call by Myers et al criteria"].astype(str).str.contains("Likely", case=False, na=False)
S4["uniprot"] = S4.Gene.map(GENE_ACC)
S4["pos"] = S4["Residue Position"].astype(int)

# ---- domains from S9
S9 = pd.read_excel(config.DATA / "external" / "grin_domains_s9.xlsx", header=2)
dom_rows = []
for gene, (c0, c1) in {"GRIN1": (1, 2), "GRIN2A": (4, 5), "GRIN2B": (7, 8)}.items():
    for _, r in S9.iterrows():
        name = str(r.iloc[0]).strip()
        try:
            lo, hi = int(r.iloc[c0]), int(r.iloc[c1])
        except (ValueError, TypeError):
            continue
        dom_rows.append(dict(gene=gene, domain=name, lo=lo, hi=hi))
DOM = pd.DataFrame(dom_rows)


def domain_of(gene, pos):
    d = DOM[(DOM.gene == gene) & (DOM.lo <= pos) & (DOM.hi >= pos)]
    return d.domain.iloc[0] if len(d) else "other"


S4["domain"] = [domain_of(g, p) for g, p in zip(S4.Gene, S4.pos)]

# ---- features
F = pd.read_csv(config.PROCESSED / "nmdar_residue_features.tsv.gz", sep="\t")
D = S4.merge(F, on=["uniprot", "pos"], how="left")
print(f"labels {len(S4)}; with structural features {D.frac_states_subunit_contact.notna().sum()}")
print(D.groupby("Gene").apply(lambda d: pd.Series({
    "n": len(d), "GoF": int(d.y.sum()), "LoF": int((1 - d.y).sum()),
    "mapped": int(d.frac_states_subunit_contact.notna().sum())}), include_groups=False).to_string())
D = D[D.frac_states_subunit_contact.notna()].copy()
print(f"\ndomain distribution:\n{D.groupby(['domain', 'y']).size().unstack(fill_value=0).to_string()}")

STATE = ["frac_states_subunit_contact", "sd_subunit_contact", "delta_open_nonactive", "mean_subunit_contact", "max_subunit_contact"]
LIGAND = ["min_dist_agonist", "min_dist_ligand"]
D["delta_open_nonactive"] = D.delta_open_nonactive.fillna(0.0)
DUM = pd.get_dummies(D.domain, prefix="dom").astype(float)
SETS = {"DOMAIN": list(DUM.columns), "LIGAND": LIGAND, "STATE": STATE, "LIGAND+STATE": LIGAND + STATE,
        "ALL": LIGAND + STATE + list(DUM.columns)}
X = pd.concat([D[LIGAND + STATE].reset_index(drop=True), DUM.reset_index(drop=True)], axis=1)
y, gene = D.y.to_numpy(), D.Gene.to_numpy()


def logo_scores(cols, mask=None):
    """Leave-one-gene-out predictions for a feature set."""
    m = np.ones(len(D), bool) if mask is None else mask
    out = np.full(len(D), np.nan)
    for g in sorted(set(gene[m])):
        tr, te = m & (gene != g), m & (gene == g)
        if tr.sum() < 10 or len(set(y[tr])) < 2 or te.sum() < 5 or len(set(y[te])) < 2:
            continue
        mod = make_pipeline(StandardScaler(), LogisticRegression(C=0.5, class_weight="balanced", max_iter=5000))
        mod.fit(X.loc[tr, cols], y[tr])
        out[te] = mod.predict_proba(X.loc[te, cols])[:, 1]
    return out


def auc_ci(yy, ss, n=2000):
    ok = ~np.isnan(ss)
    yy, ss = yy[ok], ss[ok]
    if len(set(yy)) < 2:
        return np.nan, np.nan, np.nan
    pos, neg = np.where(yy == 1)[0], np.where(yy == 0)[0]
    b = [roc_auc_score(yy[i], ss[i]) for i in (np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]) for _ in range(n))]
    return roc_auc_score(yy, ss), np.quantile(b, .025), np.quantile(b, .975)


for label, mask in [("ALL CALLS", np.ones(len(D), bool)), ("STRICT (Likely only)", D.strict.to_numpy())]:
    print(f"\n=== {label}: n={int(mask.sum())} (GoF {int(y[mask].sum())}, LoF {int((1 - y[mask]).sum())}) ===")
    preds, rows = {}, []
    for name, cols in SETS.items():
        s = logo_scores(cols, mask)
        preds[name] = s
        a, lo, hi = auc_ci(y[mask], s[mask])
        per = {}
        for g in sorted(set(gene[mask])):
            gm = mask & (gene == g)
            ok = ~np.isnan(s) & gm
            per[g] = roc_auc_score(y[ok], s[ok]) if ok.sum() > 4 and len(set(y[ok])) > 1 else np.nan
        rows.append(dict(features=name, pooled_auroc=a, lo=lo, hi=hi,
                         mean_per_gene=np.nanmean(list(per.values())), **{f"auc_{g}": v for g, v in per.items()}))
    R = pd.DataFrame(rows)
    print(R.round(3).to_string(index=False))
    # paired deltas vs the two threats
    print("\n  paired deltas (bootstrap, same variants):")
    for base in ["DOMAIN", "LIGAND"]:
        for test in ["STATE", "LIGAND+STATE"]:
            ok = mask & ~np.isnan(preds[test]) & ~np.isnan(preds[base])
            yy, s1, s2 = y[ok], preds[test][ok], preds[base][ok]
            if len(set(yy)) < 2:
                continue
            pos, neg = np.where(yy == 1)[0], np.where(yy == 0)[0]
            d = np.array([roc_auc_score(yy[i], s1[i]) - roc_auc_score(yy[i], s2[i])
                          for i in (np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]) for _ in range(2000))])
            obs = roc_auc_score(yy, s1) - roc_auc_score(yy, s2)
            p = 2 * min((d <= 0).mean(), (d >= 0).mean())
            print(f"    {test:12s} - {base:7s}: {obs:+.3f} [{np.quantile(d, .025):+.3f}, {np.quantile(d, .975):+.3f}] p {p:.3f}")
    R.to_csv(config.RESULTS / f"grin_direction_{'all' if 'ALL' in label else 'strict'}.tsv", sep="\t", index=False)

# single-feature direction signal (unsupervised, within gene)
print("\n=== single features: AUROC for GoF vs LoF, per gene (no training) ===")
sf = []
for f in STATE + LIGAND:
    r = {"feature": f}
    for g in sorted(set(gene)):
        gm = gene == g
        v = D.loc[gm, f].to_numpy(float)
        ok = ~np.isnan(v)
        r[g] = roc_auc_score(y[gm][ok], v[ok]) if len(set(y[gm][ok])) > 1 else np.nan
    sf.append(r)
print(pd.DataFrame(sf).round(3).to_string(index=False))
D.to_csv(config.RESULTS / "grin_direction_variants.tsv", sep="\t", index=False)
