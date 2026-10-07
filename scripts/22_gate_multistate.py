"""Gate G1 (ANALYSIS_PLAN.md section 2): multi-state vs single-state contacts, leave-one-state-out, shuffled-contact null.

Reuses the feature recomputation of scripts/20_ablation.py (same definitions as the main pipeline) with structures grouped
by splicing STATE (tri-snRNP = 3JCR + 6QW6). Same v2 model, held-out-unit protocol and labels (config.LABELS).
  ALL                 all 11 structures (must reproduce the main features)
  single state S      structures of state S only
  NESTED single state state chosen per held-out unit by the mean AUROC over the OTHER units (primary comparator)
  ORACLE single state best state on the test units (optimistic; descriptive)
  LOSO -S             all structures except state S
  NULL (n=200)        per structure and family, per-position contact profiles permuted among resolved positions
                      (keeps each structure's degree sequence); RNA-RNA contact graph rewired by permuting edge targets
                      within each structure (configuration model: every node keeps its in- and out-degree)
Uncertainty for ALL - NESTED: position-block bootstrap (positions resampled within units) and hierarchical bootstrap
(units, then positions) of the mean per-unit ΔAUROC; 2,000 replicates.
Usage: SNRNA_LABELS=labels_v5c.tsv python scripts/22_gate_multistate.py [--null N]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
src = open(HERE / "20_ablation.py").read().split("set_features(PDBS, True)\nfrom scipy")[0]
G = {"__file__": str(HERE / "20_ablation.py"), "__name__": "gate"}
exec(compile(src, "20_ablation", "exec"), G)
config, d, rows_for_training, design, model = G["config"], G["d"], G["rows_for_training"], G["design"], G["model"]
set_features, ctx_rows, UNITS, PDBS = G["set_features"], G["ctx_rows"], G["UNITS"], G["PDBS"]
PER0, EDGES0, orig = G["per"].copy(), G["edges"].copy(), G["orig"]
N_NULL = int(sys.argv[sys.argv.index("--null") + 1]) if "--null" in sys.argv else 200
rng = np.random.default_rng(0)

state_of = PER0.drop_duplicates("pdb_id").set_index("pdb_id").state.to_dict()
STATES = sorted(set(state_of.values()))
by_state = {s: [p for p in PDBS if state_of[p] == s] for s in STATES}


def predict():
    """Held-out predictions {(unit, ctx): DataFrame(y, s, pos)} for the current feature set."""
    out = {}
    for unit in UNITS:
        tr = rows_for_training(d[d.unit != unit])
        m = model().fit(design(d, tr.i, tr.ctx), tr.y)
        for ctx in ["dominant", "recessive"]:
            t = ctx_rows(unit, ctx)
            if t is not None:
                out[(unit, ctx)] = pd.DataFrame({"y": t.y_path.values, "pos": t.node_name.values,
                                                 "s": m.predict_proba(design(d, t.index, [ctx] * len(t)))[:, 1]})
    return out


def aucs(P):
    return {k: roc_auc_score(v.y, v.s) for k, v in P.items()}


def run(pdbs, per=None, edges=None):
    G["per"] = PER0 if per is None else per
    G["edges"] = EDGES0 if edges is None else edges
    set_features(pdbs, True)
    return predict()


# ---- configurations
preds = {"ALL": run(PDBS)}
from scipy.stats import spearmanr
print("sanity, ALL vs main features (Spearman):",
      {c: round(spearmanr(orig[c], d[c + "_r"], nan_policy="omit")[0], 3) for c in orig}, flush=True)
for s in STATES:
    preds["state:" + s] = run(by_state[s])
    preds["LOSO:-" + s] = run([p for p in PDBS if state_of[p] != s])
A = pd.DataFrame({k: aucs(v) for k, v in preds.items()}).T
keys = list(A.columns)


def mean_ctx(row, ctx):
    return np.mean([row[k] for k in keys if k[1] == ctx])


singles = [c for c in A.index if c.startswith("state:")]
nested_choice, nested_auc = {}, {}
for unit, ctx in keys:
    inner = {c: np.mean([A.loc[c, k] for k in keys if k[1] == ctx and k[0] != unit]) for c in singles}
    nested_choice[(unit, ctx)] = max(inner, key=inner.get)
    nested_auc[(unit, ctx)] = A.loc[nested_choice[(unit, ctx)], (unit, ctx)]
summ = pd.DataFrame({ctx: A.apply(lambda r: mean_ctx(r, ctx), axis=1) for ctx in ["dominant", "recessive"]})
for ctx in ["dominant", "recessive"]:
    summ.loc["NESTED single state", ctx] = np.mean([v for k, v in nested_auc.items() if k[1] == ctx])
    summ.loc["ORACLE single state", ctx] = summ.loc[singles, ctx].max()
# EXPLORATORY (added after the plan; stricter comparator): only states that resolve the test unit's paralog family
FAM = {"RNU4-2": "U4", "RNU2-2P": "U2", "RNU5B-1": "U5", "RNU6": "U6", "RNU4ATAC": "U4", "RNU12": "U2"}
fams_of_state = PER0.groupby("state").paralog_family.agg(set).to_dict()
cov_choice, cov_auc = {}, {}
for unit, ctx in keys:
    ok = [c for c in singles if FAM[unit] in fams_of_state[c.split(":", 1)[1]]]
    inner = {c: np.mean([A.loc[c, k] for k in keys if k[1] == ctx and k[0] != unit]) for c in ok}
    cov_choice[(unit, ctx)] = max(inner, key=inner.get)
    cov_auc[(unit, ctx)] = A.loc[cov_choice[(unit, ctx)], (unit, ctx)]
for ctx in ["dominant", "recessive"]:
    summ.loc["NESTED single state, covering family (exploratory)", ctx] = np.mean([v for k, v in cov_auc.items() if k[1] == ctx])
print("covering-family choices:", {f"{k[0]}/{k[1][:3]}": v for k, v in cov_choice.items()})
A.loc["NESTED single state"] = pd.Series(nested_auc)
A.loc["NESTED covering (exploratory)"] = pd.Series(cov_auc)
A.to_csv(config.result("gate_g1_per_unit.tsv"), sep="\t")
summ.to_csv(config.result("gate_g1_summary.tsv"), sep="\t")
print("\nnested choices:", {f"{k[0]}/{k[1][:3]}": v for k, v in nested_choice.items()})
print("\n=== mean held-out AUROC by configuration ===")
print(summ.round(3).sort_values("dominant", ascending=False).to_string())
print("\nper unit:")
print(A.loc[["ALL", "NESTED single state", "NESTED covering (exploratory)"] + singles].round(3).to_string())


# ---- bootstraps of ALL - NESTED single state (dominant and recessive)
def auc_safe(y, s):
    return roc_auc_score(y, s) if 0 < y.sum() < len(y) else np.nan


def boot_delta(ctx, hierarchical, choice, n=2000):
    ks = [k for k in keys if k[1] == ctx]
    data = {k: (preds["ALL"][k], preds[choice[k]][k]) for k in ks}
    pos_groups = {k: {p: np.where(a.pos.values == p)[0] for p in np.unique(a.pos.values)} for k, (a, _) in data.items()}
    out = []
    for _ in range(n):
        units = [ks[i] for i in rng.integers(0, len(ks), len(ks))] if hierarchical else ks
        ds = []
        for k in units:
            a, b = data[k]
            pg = pos_groups[k]
            names = list(pg)
            ii = np.concatenate([pg[names[j]] for j in rng.integers(0, len(names), len(names))])
            ds.append(auc_safe(a.y.values[ii], a.s.values[ii]) - auc_safe(b.y.values[ii], b.s.values[ii]))
        out.append(np.nanmean(ds))
    return np.array(out)


print("\n=== ALL minus single-state comparator: mean per-unit ΔAUROC ===")
brows = []
for comp, choice, cauc in [("NESTED single state", nested_choice, nested_auc),
                           ("NESTED covering (exploratory)", cov_choice, cov_auc),
                           ("ORACLE single state", {k: max(singles, key=lambda c: summ.loc[c, k[1]]) for k in keys},
                            {k: A.loc[max(singles, key=lambda c: summ.loc[c, k[1]]), k] for k in keys})]:
    for ctx in ["dominant", "recessive"]:
        ks = [k for k in keys if k[1] == ctx]
        obs = np.mean([A.loc["ALL", k] - cauc[k] for k in ks])
        better = sum(A.loc["ALL", k] > cauc[k] for k in ks)
        for h in [False, True]:
            bs = boot_delta(ctx, h, choice)
            brows.append(dict(comparator=comp, ctx=ctx, bootstrap="hierarchical" if h else "position-block",
                              delta=obs, lo=np.quantile(bs, .025), hi=np.quantile(bs, .975), units_better=f"{better}/{len(ks)}"))
B = pd.DataFrame(brows)
print(B.round(3).to_string(index=False), flush=True)


# ---- degree-preserving shuffled-contact null
def shuffle_per(per):
    per = per.copy()
    cols = ["n_protein_res", "n_protein_chains", "min_protein_dist", "n_snrna_contacts", "n_other_rna_contacts"]
    for _, ix in per.groupby(["pdb_id", "paralog_family"]).groups.items():
        perm = rng.permutation(len(ix))
        per.loc[ix, cols] = per.loc[ix, cols].to_numpy()[perm]
    return per


def rewire(edges):
    e = edges.copy()
    for _, ix in e.groupby("pdb_id").groups.items():
        perm = rng.permutation(len(ix))
        e.loc[ix, ["family_b", "pos_b"]] = e.loc[ix, ["family_b", "pos_b"]].to_numpy()[perm]
    return e


obs_dom = summ.loc["ALL", "dominant"]
null = []
for r in range(N_NULL):
    a = aucs(run(PDBS, shuffle_per(PER0), rewire(EDGES0)))
    null.append({ctx: np.mean([v for k, v in a.items() if k[1] == ctx]) for ctx in ["dominant", "recessive"]})
    if (r + 1) % 25 == 0:
        print(f"  null {r + 1}/{N_NULL}: dominant mean so far {np.mean([x['dominant'] for x in null]):.3f}", flush=True)
Nl = pd.DataFrame(null)
Nl.to_csv(config.result("gate_g1_null.tsv"), sep="\t", index=False)
p_null = (1 + (Nl.dominant >= obs_dom).sum()) / (1 + len(Nl))
print(f"\nNULL dominant: mean {Nl.dominant.mean():.3f}, 95th pct {Nl.dominant.quantile(.95):.3f}, max {Nl.dominant.max():.3f};"
      f" observed {obs_dom:.3f}; empirical p {p_null:.4f}")
print(f"NULL recessive: mean {Nl.recessive.mean():.3f}, 95th pct {Nl.recessive.quantile(.95):.3f}; observed {summ.loc['ALL', 'recessive']:.3f}")

# ---- gate decision (pre-specified criteria)
bd = B[(B.comparator == "NESTED single state") & (B.ctx == "dominant") & (B.bootstrap == "position-block")].iloc[0]
c1 = bd.delta > 0 and bd.lo > 0
c2 = int(bd.units_better.split("/")[0]) >= 3
c3 = obs_dom > Nl.dominant.quantile(.95)
print(f"\nG1 criteria: (1) ALL > NESTED with CI excluding 0: {c1}; (2) >=3/4 units: {c2}; (3) beats null 95th pct: {c3}")
print("G1 PASS" if (c1 and c2 and c3) else "G1 FAIL")
B.assign(null_mean=Nl.dominant.mean(), null_p95=Nl.dominant.quantile(.95), p_null=p_null, gate=c1 and c2 and c3).to_csv(
    config.result("gate_g1_tests.tsv"), sep="\t", index=False)
