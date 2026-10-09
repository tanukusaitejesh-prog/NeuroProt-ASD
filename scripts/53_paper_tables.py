"""Assemble the manuscript tables for MuSE (Multi-State Ensemble), in the format the target venue expects.

Briefings in Bioinformatics "Problem Solving Protocol" papers are judged on a fixed skeleton: a baseline table, an
ablation table, several evaluation splits, and a named method. This script produces those four tables from results
already computed under the pre-registered held-out-unit protocol (ANALYSIS_PLAN.md, committed before any analysis),
and adds the one baseline that was missing - RNA-FM, the RNA foundation model, which is the first thing a reviewer
will ask for.

  Table 1  baselines          13 competitors, mean per-unit held-out AUROC by disease mechanism, paired
                              position-block and hierarchical bootstrap of the difference
  Table 2  ablations          every single structure, the nested and oracle single-structure selections, the graph
                              term, and the two pre-gated negatives from scripts/51 and 52
  Table 3  evaluation regimes four protocols of increasing difficulty, including the timestamped prospective set
  Table 4  functional transfer zero-shot correlation with saturation genome editing - the ascertainment-free readout

Outputs results/TABLE_[1-4].md and the machine-readable .tsv beside each.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

RNG = np.random.default_rng(0)
RES = config.RESULTS
LAB = config.LABELS.replace(".tsv", "")


def auc_safe(y, s):
    y, s = np.asarray(y), np.asarray(s, dtype=float)
    ok = np.isfinite(s)
    if ok.sum() < 5 or not 0 < y[ok].sum() < ok.sum():
        return np.nan
    return roc_auc_score(y[ok], s[ok])


def block_boot(Q, col, n=2000, hierarchical=False):
    """Paired bootstrap of (MuSE - baseline) AUROC, resampling position blocks within each held-out unit."""
    units = list(Q.unit.unique())
    cache = {}
    for u in units:
        g = Q[Q.unit == u]
        pos = g.node_name.to_numpy()
        names = np.unique(pos)
        cache[u] = (g.y_path.to_numpy(), g.snrnavep_v2.to_numpy(float), g[col].to_numpy(float),
                    [np.where(pos == p)[0] for p in names])
    out = []
    for _ in range(n):
        us = [units[i] for i in RNG.integers(0, len(units), len(units))] if hierarchical else units
        ds = []
        for u in us:
            y, f, s, groups = cache[u]
            ii = np.concatenate([groups[j] for j in RNG.integers(0, len(groups), len(groups))])
            ds.append(auc_safe(y[ii], f[ii]) - auc_safe(y[ii], s[ii]))
        out.append(np.nanmean(ds))
    return np.array(out)


def md(df, path, title, note=""):
    with open(path, "w", encoding="utf8") as fh:
        fh.write(f"**{title}**\n\n")
        fh.write(df.to_markdown(index=False))
        if note:
            fh.write(f"\n\n{note}\n")
    print(f"\n{'=' * 100}\n{title}\n{'=' * 100}")
    print(df.to_string(index=False))


# ---------------------------------------------------------------------------------------------------------------
# RNA-FM: the missing foundation-model baseline, scored under the identical held-out-unit protocol
# ---------------------------------------------------------------------------------------------------------------
P = pd.read_csv(RES / f"v2_heldout_predictions_{LAB}.tsv.gz", sep="\t", low_memory=False)
FM = pd.read_csv(config.DATA / "external" / "rnafm.tsv", sep="\t")
P = P.merge(FM, on="key", how="left")
# standard zero-shot LM convention: damage = -(log P(alt) - log P(ref)); lower likelihood under the model = worse
P["rna_fm"] = -P.rnafm_dll
cov = P.rna_fm.notna().mean()
print(f"RNA-FM coverage of held-out rows: {cov:.1%} ({P.rna_fm.notna().sum()} of {len(P)})")

fm_rows = []
for ctx in ["dominant", "recessive"]:
    Q = P[P.ctx == ctx]
    per = Q.groupby("unit").apply(
        lambda g: pd.Series(dict(final=auc_safe(g.y_path, g.snrnavep_v2), rna_fm=auc_safe(g.y_path, g.rna_fm))),
        include_groups=False).dropna()
    if per.empty:
        continue
    dl = per.final - per.rna_fm
    r = dict(ctx=ctx, baseline="rna_fm", final=per.final.mean(), base=per.rna_fm.mean(), delta=dl.mean(),
             units_better=f"{(dl > 0).sum()}/{len(dl)}")
    Qf = Q[Q.rna_fm.notna()]
    for h, lab in [(False, "block"), (True, "hier")]:
        bs = block_boot(Qf, "rna_fm", hierarchical=h)
        r[lab + "_lo"], r[lab + "_hi"] = np.nanquantile(bs, .025), np.nanquantile(bs, .975)
    fm_rows.append(r)

# ---------------------------------------------------------------------------------------------------------------
# TABLE 1 - baselines
# ---------------------------------------------------------------------------------------------------------------
G2 = pd.read_csv(RES / f"gate_g2_tests_{LAB}.tsv", sep="\t")
T1raw = pd.concat([G2, pd.DataFrame(fm_rows)], ignore_index=True)

NAMES = {
    "constraint_1_minus_oe": "gnomAD constraint (1 - o/e, 10-nt window)",
    "pop_density_w2": "gnomAD allele density (+/-2 nt)",
    "pop_density_pos": "gnomAD allele density (position)",
    "trained_secondary_structure": "Trained secondary-structure model",
    "paralog_density_mitotip": "Paralogue pathogenic density (MitoTIP-style)",
    "dist_critical": "Distance to nearest known pathogenic position",
    "rfam_element_enrichment": "Rfam element-type enrichment",
    "cadd_phred": "CADD (PHRED)",
    "phylop447": "phyloP (447-way)",
    "rna_fm": "RNA-FM (foundation model, zero-shot)",
    "vienna_duplex": "ViennaRNA duplex ddG",
    "vienna_paired": "ViennaRNA pairing probability",
    "vienna_destab": "ViennaRNA fold ddG",
}
T1raw["name"] = T1raw.baseline.map(NAMES).fillna(T1raw.baseline)

rows = []
for b, g in T1raw.groupby("name", sort=False):
    r = {"Baseline": b}
    for ctx, tag in [("dominant", "Dominant"), ("recessive", "Recessive")]:
        h = g[g.ctx == ctx]
        if h.empty:
            r[f"{tag} AUROC"] = "-"
            r[f"{tag} dAUROC [95% CI]"] = "-"
        else:
            h = h.iloc[0]
            r[f"{tag} AUROC"] = f"{h.base:.3f}"
            star = "*" if h.block_lo > 0 else ""
            r[f"{tag} dAUROC [95% CI]"] = f"{h.delta:+.3f} [{h.block_lo:+.3f}, {h.block_hi:+.3f}]{star}"
    rows.append(r)
T1 = pd.DataFrame(rows)
order = T1raw[T1raw.ctx == "dominant"].sort_values("base", ascending=False).name.tolist()
T1 = T1.set_index("Baseline").reindex([o for o in order if o in T1.Baseline.values]).reset_index()
muse_dom = T1raw[T1raw.ctx == "dominant"].final.iloc[0]
muse_rec = T1raw[T1raw.ctx == "recessive"].final.iloc[0]
T1 = pd.concat([pd.DataFrame([{"Baseline": "**MuSE (this work)**", "Dominant AUROC": f"**{muse_dom:.3f}**",
                               "Dominant dAUROC [95% CI]": "-", "Recessive AUROC": f"**{muse_rec:.3f}**",
                               "Recessive dAUROC [95% CI]": "-"}]), T1], ignore_index=True)
md(T1, RES / "TABLE_1_baselines.md",
   f"Table 1. MuSE against {len(T1) - 1} baselines under the pre-registered held-out-unit protocol.",
   "Mean AUROC across held-out units, by disease mechanism. dAUROC is MuSE minus the baseline, averaged over "
   "units, with a paired position-block bootstrap 95% CI (2,000 resamples). * marks intervals excluding zero. "
   "Hierarchical (unit-then-position) intervals are in Supplementary Table S1 and are wider throughout, "
   "reflecting that only four units contribute to each context.")
T1raw.to_csv(RES / "TABLE_1_baselines.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------------------------------------------
# TABLE 2 - ablations
# ---------------------------------------------------------------------------------------------------------------
A = pd.read_csv(RES / "ablation_summary_labels_v5.tsv", sep="\t", index_col=0)
STATE = {"3JCR": "tri-snRNP", "6QW6": "tri-snRNP", "6QX9": "pre-B", "5O9Z": "B", "6FF7": "Bact",
         "5Z56": "Bact-mature", "5XJC": "C*", "6QDV": "P", "6Y5Q": "U2-snRNP", "7DVQ": "minor-Bact",
         "8Y6O": "minor-preB"}
full_d, full_r = A.loc["ALL (11 structures)", "dominant"], A.loc["ALL (11 structures)", "recessive"]
rows = [dict(Ablation="**None (full MuSE, 11 structures)**", Dominant=f"**{full_d:.3f}**",
             Recessive=f"**{full_r:.3f}**", dDominant="-")]
for k, lab in [("ALL, no graph", "Remove the contact-graph smoothing term"),
               ("NESTED single structure", "Single structure, chosen without the test gene"),
               ("ORACLE best single structure", "Single structure, chosen WITH the test gene (oracle)")]:
    rows.append(dict(Ablation=lab, Dominant=f"{A.loc[k, 'dominant']:.3f}", Recessive=f"{A.loc[k, 'recessive']:.3f}",
                     dDominant=f"{A.loc[k, 'dominant'] - full_d:+.3f}"))
for pdb in sorted(STATE, key=lambda p: -A.loc[p, "dominant"]):
    rows.append(dict(Ablation=f"Only {pdb} ({STATE[pdb]})", Dominant=f"{A.loc[pdb, 'dominant']:.3f}",
                     Recessive=f"{A.loc[pdb, 'recessive']:.3f}",
                     dDominant=f"{A.loc[pdb, 'dominant'] - full_d:+.3f}"))
# the two pre-gated additions that FAILED, reported as negatives
TG = pd.read_csv(RES / "trajectory_gate_observed.tsv", sep="\t").set_index("test").observed
PI = pd.read_csv(RES / "partner_identity_gate.tsv", sep="\t").set_index("model").logo_auroc
rows += [dict(Ablation="ADD conformational-trajectory features (scripts/51)", Dominant="n.s.", Recessive="n.s.",
              dDominant=f"{TG['logo_auroc_gain']:+.3f}, but fails its degree-matched null (p = 1.00) and the gain "
                        f"is one gene with 7 positives"),
         dict(Ablation="ADD protein-partner identity (scripts/52)", Dominant="n.s.", Recessive="n.s.",
              dDominant=f"{PI['burial+degree+identity'] - PI['burial+degree']:+.3f} (does not transfer across genes)")]
T2 = pd.DataFrame(rows)
md(T2, RES / "TABLE_2_ablations.md", "Table 2. Ablation of every component and of each conformational state.",
   "Mean AUROC across held-out units. The final two rows are components that were specified with a kill criterion "
   "in advance, tested, and rejected; they are reported because a negative result about an obvious extension is "
   "evidence about the problem, not an omission.")
T2.to_csv(RES / "TABLE_2_ablations.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------------------------------------------
# TABLE 3 - evaluation regimes
# ---------------------------------------------------------------------------------------------------------------
FS = pd.read_csv(RES / f"final_system_tests_{LAB}.tsv", sep="\t")
pooled = FS[FS.comparison.str.contains("rank-pooled AUROC")].set_index("ctx")
rows = [dict(Regime="Held-out gene, dominant mechanism", What="Every variant of one gene withheld entirely",
             n=int(pooled.loc["dominant", "n"]), Positives=int(pooled.loc["dominant", "n_pos"]),
             AUROC=f"{pooled.loc['dominant', 'delta']:.3f} {pooled.loc['dominant', 'ci']}"),
        dict(Regime="Held-out gene, recessive mechanism", What="Same, recessive positives",
             n=int(pooled.loc["recessive", "n"]), Positives=int(pooled.loc["recessive", "n_pos"]),
             AUROC=f"{pooled.loc['recessive', 'delta']:.3f} {pooled.loc['recessive', 'ci']}"),
        dict(Regime="Held-out gene, pooled", What="Both mechanisms, rank-pooled",
             n=int(pooled.loc["both", "n"]), Positives=int(pooled.loc["both", "n_pos"]),
             AUROC=f"{pooled.loc['both', 'delta']:.3f} {pooled.loc['both', 'ci']}")]
_cf = pd.read_csv(RES / f"constraint_vs_function_rho_{LAB}.tsv", sep="\t")
_mu = _cf[_cf.score == "snrnavep_v2"].iloc[0]
rows.append(dict(Regime="Zero-shot functional transfer",
                 What="RNU4-2 saturation genome editing; no functional measurement used in training",
                 n=int(_mu.n), Positives="continuous",
                 AUROC=f"Spearman {_mu.rho:.3f} [{_mu.lo:.3f}, {_mu.hi:.3f}] (Table 4)"))
dep = sorted(Path(RES).glob("prospective_deposit_*.tsv"))
if dep:
    D = pd.read_csv(dep[-1], sep="\t")
    rows.append(dict(Regime="Prospective, timestamped", What=f"Predictions deposited {dep[-1].stem.split('_')[-1]} "
                                                             f"with SHA-256, evaluation protocol fixed in advance",
                     n=len(D), Positives="pending", AUROC="to be evaluated on future reports"))
T3 = pd.DataFrame(rows)
md(T3, RES / "TABLE_3_regimes.md", "Table 3. Four evaluation regimes of increasing difficulty.",
   "Held-out-gene is strictly harder than the cluster split used by comparable methods papers: no variant, "
   "position or paralogue of the test gene is seen in training. The prospective set is deposited and hashed, "
   "so it cannot be revised after the fact.")
T3.to_csv(RES / "TABLE_3_regimes.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------------------------------------------
# TABLE 4 - zero-shot transfer to measured function (the headline)
# ---------------------------------------------------------------------------------------------------------------
CF = pd.read_csv(RES / f"constraint_vs_function_rho_{LAB}.tsv", sep="\t")
PR = pd.read_csv(RES / f"constraint_vs_function_paired_{LAB}.tsv", sep="\t")
md(CF.round(3), RES / "TABLE_4_function.md",
   "Table 4. Zero-shot Spearman correlation with measured function (RNU4-2 saturation genome editing).",
   "No model saw any functional measurement. This readout is free of the clinical-ascertainment circularity that "
   "affects population-derived scores, because the assay does not consult allele frequency.\n\n"
   + PR.round(3).to_markdown(index=False))
print("\nPaired differences:")
print(PR.round(3).to_string(index=False))

print(f"\n\nwrote TABLE_1..4 (.md and .tsv) to {RES}")
