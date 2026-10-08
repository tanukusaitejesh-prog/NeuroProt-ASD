"""The decisive test for the snRNA paper's framing (EXPLORATORY; added after gate G2 failed, labelled as post-hoc).

Gate G2 showed that gnomAD population constraint and variant density are not separable from the structural model on
CLINICAL labels. This asks the question that decides whether the paper is a method paper or a resource paper:
do those population features also track MEASURED function?

RNU4-2 saturation genome editing (De Jonghe et al. 2026), 485 variants with a structure node. The structural model is
scored zero-shot (no RNU4-2 and no U4-family labels in training). Population features cannot be circular here: the SGE
score is a wet-lab readout, not a curated label derived from population data.

Reported: Spearman with -SGE (higher = more damaging) with position-block bootstrap CIs; paired bootstrap of the
difference between the structural model and each population feature; and the partial Spearman of structure given
constraint (rank residuals), which asks whether structure adds information beyond depletion.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

rng = np.random.default_rng(0)
src = open(Path(__file__).resolve().parent / "11_v2_model.py").read().split("results, pooled = [], []")[0]
G = {"__file__": str(Path(__file__).resolve().parent / "11_v2_model.py"), "__name__": "cvf"}
exec(compile(src, "11_v2_model", "exec"), G)
d, rows_for_training, design, model = G["d"], G["rows_for_training"], G["design"], G["model"]

# population features: per-position gnomAD density (own allele excluded) and local observed/expected constraint
d["n_int"] = np.floor(d.n_pos).astype(int)
gn = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")
gn = gn[(gn["filter"] == "PASS") & (gn.AC > 0) & (gn.n_pos >= 1)].copy()
gn["key"] = gn.chrom + ":" + gn.pos.astype(str) + ":" + gn.ref + ":" + gn.alt
gn["n_int"] = np.floor(gn.n_pos).astype(int)
cnt = gn.groupby(["gene_name", "n_int"]).size()
self_ = d.key.isin(set(gn.key)).astype(int)
d["pop_density_pos"] = -(np.array([cnt.get((g, p), 0) for g, p in zip(d.gene_name, d.n_int)]) - self_)
d["pop_density_w2"] = -(np.array([sum(cnt.get((g, p + k), 0) for k in range(-2, 3)) for g, p in zip(d.gene_name, d.n_int)]) - self_)
d["constraint_1_minus_oe"] = 1 - d.rel_oe_w10

# structural model, zero-shot (no U4 family in training)
tr = rows_for_training(d[~d.unit.isin(["RNU4-2", "RNU4ATAC"])])
m = model().fit(design(d, tr.i, tr.ctx), tr.y)
S = d[(d.gene_name == "RNU4-2") & d.sge_score.notna()].copy()
S["snrnavep_v2"] = m.predict_proba(design(d, S.index, ["dominant"] * len(S)))[:, 1]
S["neg_sge"] = -S.sge_score
SCORES = ["snrnavep_v2", "nb_frac_prot", "frac_states_protein_contact",
          "constraint_1_minus_oe", "pop_density_pos", "pop_density_w2", "cadd_phred", "phylop447"]
groups = [np.where(S.node_name.to_numpy() == p)[0] for p in S.node_name.unique()]


def boot_idx(n=2000):
    for _ in range(n):
        yield np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])


print(f"RNU4-2 SGE variants with a structure node: {len(S)}\n")
print("=== Spearman with -SGE function score [position-block 95% CI] ===")
rows = []
y = S.neg_sge.to_numpy(float)
for s in SCORES:
    v = S[s].to_numpy(float)
    rho = spearmanr(v, y, nan_policy="omit")[0]
    b = [spearmanr(v[i], y[i], nan_policy="omit")[0] for i in boot_idx(1000)]
    rows.append(dict(score=s, n=int((~np.isnan(v)).sum()), rho=rho,
                     lo=np.nanquantile(b, .025), hi=np.nanquantile(b, .975)))
R = pd.DataFrame(rows)
print(R.round(3).to_string(index=False))

print("\n=== paired difference: structural model minus population feature ===")
pr = []
for base in ["constraint_1_minus_oe", "pop_density_pos", "pop_density_w2", "cadd_phred", "phylop447"]:
    for test in ["snrnavep_v2", "nb_frac_prot"]:
        a, bb = S[test].to_numpy(float), S[base].to_numpy(float)
        ok = ~(np.isnan(a) | np.isnan(bb) | np.isnan(y))
        obs = spearmanr(a[ok], y[ok])[0] - spearmanr(bb[ok], y[ok])[0]
        dd = []
        for i in boot_idx(2000):
            i = i[ok[i]]
            if len(i) > 20:
                dd.append(spearmanr(a[i], y[i])[0] - spearmanr(bb[i], y[i])[0])
        dd = np.array(dd)
        p = 2 * min((dd <= 0).mean(), (dd >= 0).mean())
        pr.append(dict(structure=test, population=base, delta_rho=obs,
                       lo=np.quantile(dd, .025), hi=np.quantile(dd, .975), p=min(p, 1.0)))
P = pd.DataFrame(pr)
print(P.round(3).to_string(index=False))


def partial(x, z, yy):
    """Spearman of x with yy after removing the rank-linear effect of z from both."""
    ok = ~(np.isnan(x) | np.isnan(z) | np.isnan(yy))
    rx, rz, ry = rankdata(x[ok]), rankdata(z[ok]), rankdata(yy[ok])
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    return np.corrcoef(ex, ey)[0, 1]


print("\n=== partial correlation: does each add information beyond the other? ===")
for st in ["snrnavep_v2", "nb_frac_prot"]:
    for pop in ["constraint_1_minus_oe", "pop_density_w2"]:
        a = partial(S[st].to_numpy(float), S[pop].to_numpy(float), y)
        b = partial(S[pop].to_numpy(float), S[st].to_numpy(float), y)
        print(f"  {st:20s} | {pop:22s}  structure given population {a:+.3f}   population given structure {b:+.3f}")

R.to_csv(config.result("constraint_vs_function_rho.tsv"), sep="\t", index=False)
P.to_csv(config.result("constraint_vs_function_paired.tsv"), sep="\t", index=False)
S[["key", "hgvs", "n_pos", "sge_score"] + SCORES].to_csv(config.result("constraint_vs_function_variants.tsv"), sep="\t", index=False)
