"""Gate G2 (ANALYSIS_PLAN.md section 2): stronger baselines under the identical held-out-unit protocol.

Final system = v2-nested held-out predictions (results/v2_heldout_predictions*, v2_nested*). Baselines (higher = more damaging):
  vienna_*                     ViennaRNA folding ddG, pairing, duplex ddG (unsupervised); trained_secondary_structure
                               (v2 design and regularisation, contact features replaced by the ViennaRNA features)
  rfam_element_enrichment      annotated element type (Rfam WUSS projection + Sm site; snrna_vep/elements.py); score =
                               log enrichment of the element type among pathogenic variants of the same context in the
                               TRAINING units (pathogenic share / share of all possible variants in that unit), averaged
  pop_density_pos / _w2        minus the number of OTHER gnomAD v4.1 PASS alleles at the position / within +/-2 nt
                               (the variant's own allele excluded)
  paralog_density_mitotip      training-unit pathogenic variants (any mode) at the aligned family position +/-1 nt
  dist_critical                minus the distance (family coordinates) to the nearest training-unit pathogenic position
                               (any mode); units without a same-family training pathogenic variant are tied (AUROC 0.5)
Uncertainty: mean per-unit ΔAUROC (final - baseline) with position-block and hierarchical (units, then positions) bootstrap.
G2 decision rule (pre-specified): fails if pop_density reaches mean dominant AUROC >= final - 0.02, or the position-block CI
of (final - pop_density) includes 0.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
src = open(HERE / "21_robustness.py").read().split("Pv = pd.read_csv")[0]
G = {"__file__": str(HERE / "21_robustness.py"), "__name__": "g2"}
exec(compile(src, "21_robustness", "exec"), G)
config, d, fit_scores, SS, metrics, UNITS = G["config"], G["d"], G["fit_scores"], G["SS"], G["metrics"], G["UNITS"]
from snrna_vep.elements import annotate

rng = np.random.default_rng(0)
d["n_int"] = np.floor(d.n_pos).astype(int)

# ---- element annotation
feat = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t", low_memory=False)
E = annotate(feat, sorted(d.gene_name.unique()))
E.to_csv(config.PROCESSED / "rfam_elements.tsv", sep="\t", index=False)
d = d.merge(E[["gene_name", "n_int", "element"]], on=["gene_name", "n_int"], how="left")
d["element"] = d.element.fillna("unaligned")

# ---- population density (gnomAD v4.1 PASS, transcribed region), own allele excluded
gn = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")
gn = gn[(gn["filter"] == "PASS") & (gn.AC > 0) & (gn.n_pos >= 1)].copy()
gn["key"] = gn.chrom + ":" + gn.pos.astype(str) + ":" + gn.ref + ":" + gn.alt
gn["n_int"] = np.floor(gn.n_pos).astype(int)
cnt = gn.groupby(["gene_name", "n_int"]).size()
in_gn = set(gn.key)
self_ = d.key.isin(in_gn).astype(int)


def count_at(g, p):
    return cnt.get((g, p), 0)


d["pop_density_pos"] = -(np.array([count_at(g, p) for g, p in zip(d.gene_name, d.n_int)]) - self_)
d["pop_density_w2"] = -(np.array([sum(count_at(g, p + k) for k in range(-2, 3)) for g, p in zip(d.gene_name, d.n_int)]) - self_)


# ---- leave-unit-out positional transfer baselines
def ctx_mask(df, ctx):
    pref = "AD-" if ctx == "dominant" else "AR-"
    return (df.y_path == 1) & df["mode"].str.contains(pref)


def transfer_scores(unit, ctx, t):
    tr = d[d.unit != unit]
    # element enrichment (same context), averaged over training units with >= 3 pathogenic in that context
    enr = []
    for u, g in tr.groupby("unit"):
        pth = g[ctx_mask(g, ctx)]
        if len(pth) < 3:
            continue
        allv = d[(d.unit == u)].element.value_counts(normalize=True)
        pv = pth.element.value_counts(normalize=True)
        enr.append(np.log((pv.reindex(allv.index).fillna(0) + 0.01) / (allv + 0.01)))
    e = pd.concat(enr, axis=1).mean(1) if enr else pd.Series(dtype=float)
    out = pd.DataFrame(index=t.index)
    out["rfam_element_enrichment"] = t.element.map(e).fillna(0.0)
    # MitoTIP-style paralog density and distance to critical region (any mode, same paralog family)
    P = tr[tr.y_path == 1][["paralog_family", "family_ref_pos"]].dropna()
    dens, dist = [], []
    for fam, pos in zip(t.paralog_family, t.family_ref_pos):
        q = P[P.paralog_family == fam].family_ref_pos.to_numpy()
        dens.append(((q >= pos - 1) & (q <= pos + 1)).sum())
        dist.append(-np.abs(q - pos).min() if len(q) else 0.0)
    out["paralog_density_mitotip"] = dens
    out["dist_critical"] = dist
    return out


Pv = pd.read_csv(config.result("v2_heldout_predictions.tsv.gz"), sep="\t", low_memory=False)
N = pd.read_csv(config.result("v2_nested.tsv"), sep="\t")
choice = {(u, c): s for u, c, s in zip(N.unit, N.ctx, N.chosen)}
BASE = ["trained_secondary_structure", "vienna_destab", "vienna_paired", "vienna_duplex", "rfam_element_enrichment",
        "pop_density_pos", "pop_density_w2", "paralog_density_mitotip", "dist_critical", "constraint_1_minus_oe",
        "cadd_phred", "phylop447"]
rows, pooled = [], []
for unit in UNITS:
    for ctx, pref in [("dominant", "AD-"), ("recessive", "AR-")]:
        t = d[(d.unit == unit) & ((d.y_path == 0) | ((d.y_path == 1) & d["mode"].str.contains(pref)))].copy()
        if (t.y_path == 1).sum() < 3 or (t.y_path == 0).sum() < 3:
            continue
        pv = Pv[(Pv.unit == unit) & (Pv.ctx == ctx)].set_index("key")
        t["final"] = t.key.map(pv[choice[(unit, ctx)]])
        t["trained_secondary_structure"] = fit_scores([c + "_r" for c in SS], d[d.unit != unit], t, ctx)
        t = t.join(transfer_scores(unit, ctx, t))
        for s in ["final"] + BASE:
            rows.append(dict(unit=unit, ctx=ctx, score=s, auroc=metrics(t.y_path, t[s])["auroc"], n_pos=int(t.y_path.sum()), n=len(t)))
        pooled.append(t.assign(ctx=ctx))
Bt = pd.DataFrame(rows)
Bt.to_csv(config.result("gate_g2_baselines.tsv"), sep="\t", index=False)
print("=== held-out AUROC per unit ===")
print(Bt.pivot_table(index="score", columns=["ctx", "unit"], values="auroc").round(3).to_string())
summ = Bt.groupby(["ctx", "score"]).auroc.mean().unstack(0).round(3).sort_values("dominant", ascending=False)
print("\n=== mean over units ===")
print(summ.to_string(), flush=True)

# ---- clustered bootstrap of final - baseline
P = pd.concat(pooled)


def auc_safe(y, s):
    ok = ~np.isnan(s)
    y, s = y[ok], s[ok]
    return roc_auc_score(y, s) if 0 < y.sum() < len(y) else np.nan


def boot(Q, b, hierarchical, n=2000):
    units = list(Q.unit.unique())
    cache = {}
    for u in units:
        g = Q[Q.unit == u]
        pos = g.node_name.to_numpy()
        names = np.unique(pos)
        cache[u] = (g.y_path.to_numpy(), g.final.to_numpy(float), g[b].to_numpy(float), [np.where(pos == p)[0] for p in names])
    out = []
    for _ in range(n):
        us = [units[i] for i in rng.integers(0, len(units), len(units))] if hierarchical else units
        ds = []
        for u in us:
            y, f, s, groups = cache[u]
            ii = np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])
            ds.append(auc_safe(y[ii], f[ii]) - auc_safe(y[ii], s[ii]))
        out.append(np.nanmean(ds))
    return np.array(out)


res = []
for ctx in ["dominant", "recessive"]:
    Q = P[P.ctx == ctx]
    for b in BASE:
        per = Bt[(Bt.ctx == ctx)].pivot(index="unit", columns="score", values="auroc")
        dl = per["final"] - per[b]
        r = dict(ctx=ctx, baseline=b, final=per["final"].mean(), base=per[b].mean(), delta=dl.mean(), units_better=f"{(dl > 0).sum()}/{len(dl)}")
        for h, lab in [(False, "block"), (True, "hier")]:
            bs = boot(Q, b, h)
            r[lab + "_lo"], r[lab + "_hi"] = np.nanquantile(bs, .025), np.nanquantile(bs, .975)
        res.append(r)
R2 = pd.DataFrame(res)
R2.to_csv(config.result("gate_g2_tests.tsv"), sep="\t", index=False)
print("\n=== final minus baseline (mean per-unit ΔAUROC; 95% CIs: position-block / hierarchical) ===")
print(R2.round(3).to_string(index=False))

pdn = R2[(R2.ctx == "dominant") & (R2.baseline == "pop_density_pos")].iloc[0]
pdw = R2[(R2.ctx == "dominant") & (R2.baseline == "pop_density_w2")].iloc[0]
fail = any(r.base >= r.final - 0.02 or r.block_lo <= 0 for r in [pdn, pdw])
print(f"\nG2 (population density): {'FAIL' if fail else 'PASS'}")
close = R2[(R2.ctx == "dominant") & (R2.delta <= 0.02)]
print("dominant baselines within 0.02 of the final system:", list(close.baseline) or "none")
