"""Rigor analyses (ANALYSIS_PLAN.md section 3).

(a) Per-unit AUROC with position-block 95% CIs (final system and baselines), and a logistic GLMM with a unit random
    intercept: y ~ within-unit rank(score) + (1 | unit), variational Bayes (statsmodels BinomialBayesMixedGLM).
(b) Dominant vs recessive interaction within RNU4-2 and RNU2-2: stacked rows (dominant cases + the gene's controls,
    recessive cases + the gene's controls), y ~ (contact + graph-contact + phyloP) x recessive + gene; position-clustered
    SEs; joint Wald test of the 3 interaction terms. Also the case-only contrast (dominant vs recessive pathogenic).
(c) Family-projection confound: v2 retrained WITHOUT the minor-spliceosome genes (RNU4ATAC, RNU12, RNU6ATAC) in training;
    held-out AUROC of the major units. (Minor-spliceosome STRUCTURES removed: see gate G1 LOSO rows.)
(d) Gradient: Spearman with the continuous RNU4-2 SGE function score (zero-shot; no RNU4-2 labels in training), with
    position-block bootstrap CIs; and with SGE score among variants in the dominant critical region only.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
src = open(HERE / "11_v2_model.py").read().split("results, pooled = [], []")[0]
G = {"__file__": str(HERE / "11_v2_model.py"), "__name__": "rig"}
exec(compile(src, "11_v2_model", "exec"), G)
config, d, rows_for_training, design, model, R = G["config"], G["d"], G["rows_for_training"], G["design"], G["model"], G["R"]
rng = np.random.default_rng(0)
out_lines = []


def log(s=""):
    print(s, flush=True)
    out_lines.append(str(s))


def block_ci(y, s, pos, n=2000):
    ok = ~np.isnan(s)
    y, s, pos = y[ok], s[ok], pos[ok]
    if not 0 < y.sum() < len(y):
        return np.nan, np.nan, np.nan
    groups = [np.where(pos == p)[0] for p in np.unique(pos)]
    v = []
    for _ in range(n):
        ii = np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])
        if 0 < y[ii].sum() < len(ii):
            v.append(roc_auc_score(y[ii], s[ii]))
    return roc_auc_score(y, s), np.quantile(v, .025), np.quantile(v, .975)


# ---- (a) per-unit CIs + GLMM
B = pd.read_csv(config.result("gate_g2_baselines.tsv"), sep="\t") if config.result("gate_g2_baselines.tsv").exists() else None
Pv = pd.read_csv(config.result("v2_heldout_predictions.tsv.gz"), sep="\t", low_memory=False)
N = pd.read_csv(config.result("v2_nested.tsv"), sep="\t")
Pv["final"] = [r[c] for r, c in zip(Pv.to_dict("records"), [dict(zip(zip(N.unit, N.ctx), N.chosen))[(u, c)] for u, c in zip(Pv.unit, Pv.ctx)])]
rows = []
for (unit, ctx), g in Pv.groupby(["unit", "ctx"]):
    for s in ["final", "frac_states_protein_contact", "phylop447", "cadd_phred"]:
        a, lo, hi = block_ci(g.y_path.to_numpy(), g[s].to_numpy(float), g.node_name.to_numpy())
        rows.append(dict(unit=unit, ctx=ctx, score=s, n_pos=int(g.y_path.sum()), n=len(g), auroc=a, lo=lo, hi=hi))
U = pd.DataFrame(rows)
U.to_csv(config.result("rigor_per_unit_ci.tsv"), sep="\t", index=False)
log("=== (a) per-unit AUROC [position-block 95% CI] ===")
U["txt"] = U.apply(lambda r: f"{r.auroc:.2f} [{r.lo:.2f}, {r.hi:.2f}]", axis=1)
log(U.pivot_table(index=["ctx", "unit"], columns="score", values="txt", aggfunc="first").to_string())

from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
log("\n=== (a) GLMM: y ~ z(within-unit rank of score) + (1 | unit); posterior mean (SD) of the score effect, log-odds per SD ===")
gl = []
for ctx in ["dominant", "recessive"]:
    Q = Pv[Pv.ctx == ctx].copy()
    for s in ["final", "frac_states_protein_contact", "phylop447", "cadd_phred"]:
        q = Q[Q[s].notna()].copy()
        q["z"] = q.groupby("unit")[s].rank(pct=True)
        q["z"] = (q.z - q.z.mean()) / q.z.std()
        q["u"] = q.unit.astype("category").cat.codes
        m = BinomialBayesMixedGLM.from_formula("y_path ~ z", {"unit": "0 + C(unit)"}, q).fit_vb()
        i = list(m.model.exog_names).index("z")
        gl.append(dict(ctx=ctx, score=s, beta=m.fe_mean[i], sd=m.fe_sd[i], n=len(q), units=q.unit.nunique()))
GL = pd.DataFrame(gl)
GL["OR_per_SD"] = np.exp(GL.beta)
GL["OR_95"] = GL.apply(lambda r: f"[{np.exp(r.beta - 1.96 * r.sd):.2f}, {np.exp(r.beta + 1.96 * r.sd):.2f}]", axis=1)
GL.to_csv(config.result("rigor_glmm.tsv"), sep="\t", index=False)
log(GL.round(3).to_string(index=False))

# ---- (b) dominant vs recessive interaction within genes with both modes
log("\n=== (b) dominant vs recessive within RNU4-2 and RNU2-2 ===")
F3 = ["frac_states_protein_contact_r", "nb_frac_prot_r", "phylop447_r"]
st = []
for gene in ["RNU4-2", "RNU2-2P"]:
    g = d[d.gene_name == gene]
    ctl = g[g.y_path == 0]
    dom = g[(g.y_path == 1) & g["mode"].str.contains("AD-") & ~g["mode"].str.contains("AR-")]
    rec = g[(g.y_path == 1) & g["mode"].str.contains("AR-") & ~g["mode"].str.contains("AD-")]
    st.append(pd.concat([dom.assign(rec=0), ctl.assign(rec=0), rec.assign(rec=1), ctl.assign(rec=1)]).assign(gene=gene))
    log(f"  {gene}: dominant {len(dom)}, recessive {len(rec)}, controls {len(ctl)}")
S = pd.concat(st)
X = S[F3].copy()
for c in F3:
    X[c + ":rec"] = S[c] * S.rec
X["rec"] = S.rec
X["gene_RNU2"] = (S.gene == "RNU2-2P").astype(float)
X["gene_RNU2:rec"] = X.gene_RNU2 * S.rec
X = sm.add_constant(X.astype(float))
fit = sm.GLM(S.y_path.astype(float), X, family=sm.families.Binomial()).fit(
    cov_type="cluster", cov_kwds={"groups": (S.gene + S.node_name).astype("category").cat.codes})
log(pd.DataFrame({"coef": fit.params, "se": fit.bse, "p": fit.pvalues}).round(4).to_string())
inter = [c + ":rec" for c in F3]
Rm = np.zeros((3, len(fit.params)))
for k, c in enumerate(inter):
    Rm[k, list(fit.params.index).index(c)] = 1
w = fit.wald_test(Rm, scalar=True)
log(f"  joint Wald test of the 3 feature x recessive interactions: chi2 {float(w.statistic):.2f}, p {float(w.pvalue):.2e}")
crit = (fit.params["frac_states_protein_contact_r:rec"] < 0) and (fit.params["phylop447_r:rec"] > 0) and float(w.pvalue) < 0.05
log(f"  pre-specified mechanism criterion (joint p<0.05, contact x rec < 0, phyloP x rec > 0): {'MET' if crit else 'NOT MET'}")
# per gene
for gene in ["RNU4-2", "RNU2-2P"]:
    s_ = S[S.gene == gene]
    Xg = s_[F3].copy()
    for c in F3:
        Xg[c + ":rec"] = s_[c] * s_.rec
    Xg["rec"] = s_.rec
    Xg = sm.add_constant(Xg.astype(float))
    fg = sm.GLM(s_.y_path.astype(float), Xg, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": s_.node_name.astype("category").cat.codes})
    Rg = np.zeros((3, len(fg.params)))
    for k, c in enumerate(inter):
        Rg[k, list(fg.params.index).index(c)] = 1
    wg = fg.wald_test(Rg, scalar=True)
    log(f"  {gene}: contact x rec {fg.params[inter[0]]:+.2f} (p {fg.pvalues[inter[0]]:.3f}), graph x rec {fg.params[inter[1]]:+.2f} "
        f"(p {fg.pvalues[inter[1]]:.3f}), phyloP x rec {fg.params[inter[2]]:+.2f} (p {fg.pvalues[inter[2]]:.3f}); joint p {float(wg.pvalue):.2e}")
pd.DataFrame({"coef": fit.params, "se": fit.bse, "p": fit.pvalues}).to_csv(config.result("rigor_dom_rec_interaction.tsv"), sep="\t")

# ---- (c) family-projection confound: no minor-spliceosome genes in training
log("\n=== (c) v2 trained without minor-spliceosome genes (RNU4ATAC, RNU12, RNU6ATAC) ===")
MINOR = {"RNU4ATAC", "RNU12", "RNU6ATAC"}
rc = []
for unit in ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6"]:
    for excl in [False, True]:
        tr = d[(d.unit != unit) & (~d.unit.isin(MINOR) if excl else True)]
        tr = rows_for_training(tr)
        m = model().fit(design(d, tr.i, tr.ctx), tr.y)
        for ctx, pref in [("dominant", "AD-"), ("recessive", "AR-")]:
            t = d[(d.unit == unit) & ((d.y_path == 0) | ((d.y_path == 1) & d["mode"].str.contains(pref)))]
            if (t.y_path == 1).sum() < 3 or (t.y_path == 0).sum() < 3:
                continue
            p = m.predict_proba(design(d, t.index, [ctx] * len(t)))[:, 1]
            rc.append(dict(unit=unit, ctx=ctx, minor_excluded=excl, auroc=roc_auc_score(t.y_path, p)))
RC = pd.DataFrame(rc)
RC.to_csv(config.result("rigor_minor_excluded.tsv"), sep="\t", index=False)
log(RC.pivot_table(index=["ctx", "unit"], columns="minor_excluded", values="auroc").round(3).to_string())
log(RC.groupby(["ctx", "minor_excluded"]).auroc.mean().round(3).unstack().to_string())

# ---- (d) gradient vs SGE (zero-shot: RNU4-2 and the whole U4 family out of training)
log("\n=== (d) Spearman with -SGE function score, RNU4-2 (no U4-family labels in training) ===")
tr = rows_for_training(d[~d.unit.isin(["RNU4-2", "RNU4ATAC"])])
m = model().fit(design(d, tr.i, tr.ctx), tr.y)
s4 = d[(d.gene_name == "RNU4-2") & d.sge_score.notna()].copy()
s4["v2_dom"] = m.predict_proba(design(d, s4.index, ["dominant"] * len(s4)))[:, 1]
s4["neg_sge"] = -s4.sge_score
groups = [np.where(s4.node_name.to_numpy() == p)[0] for p in s4.node_name.unique()]
sg = []
for sub_name, sub in [("all SGE variants", s4), ("SNVs", s4[s4.vtype == "snv"]), ("n.56-90 (dominant region)", s4[(s4.n_pos >= 56) & (s4.n_pos <= 90)])]:
    gr = [np.where(sub.node_name.to_numpy() == p)[0] for p in sub.node_name.unique()]
    for s in ["v2_dom", "nb_frac_prot", "frac_states_protein_contact", "cadd_phred", "phylop447"]:
        q = sub[[s, "neg_sge"]].to_numpy(float)
        rho = spearmanr(q[:, 0], q[:, 1], nan_policy="omit")[0]
        bs = []
        for _ in range(1000):
            ii = np.concatenate([gr[j] for j in rng.integers(0, len(gr), len(gr))])
            bs.append(spearmanr(q[ii, 0], q[ii, 1], nan_policy="omit")[0])
        sg.append(dict(subset=sub_name, score=s, n=int((~np.isnan(q).any(1)).sum()), rho=rho, lo=np.nanquantile(bs, .025), hi=np.nanquantile(bs, .975)))
SG = pd.DataFrame(sg)
SG.to_csv(config.result("rigor_sge_gradient.tsv"), sep="\t", index=False)
log(SG.round(3).to_string(index=False))
s4[["key", "hgvs", "vtype", "n_pos", "sge_score", "v2_dom", "nb_frac_prot", "frac_states_protein_contact", "cadd_phred", "phylop447"]].to_csv(
    config.result("rigor_sge_zeroshot_predictions.tsv"), sep="\t", index=False)
open(config.result("log_24.txt"), "w", encoding="utf8").write("\n".join(out_lines))
