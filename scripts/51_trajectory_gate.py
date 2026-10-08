"""GATE T1. Do spliceosomal disease variants sit at conformational HINGES rather than at static cores?

Motivation. Gate G1 established that pooling 11 cryo-EM states beats any single state, but it pooled them as an
unordered bag: the feature was "fraction of states in which this nucleotide touches a protein". That throws away the
one thing the spliceosome uniquely provides - the states are an ORDERED assembly pathway, and a nucleotide has a
trajectory through it. U4 is bound in tri-snRNP, pre-B and B, and is ejected at activation; U5 loop I holds the exons
through both transesterifications; U2 recognises the branch point early. These are different jobs at different times.

Hypothesis. Pathogenic variants concentrate not at the most buried nucleotides but at nucleotides whose protein
partners TURN OVER between consecutive states - the positions that have to let go of one partner and acquire another
for the cycle to proceed. A permanently buried nucleotide is a constitutive scaffold; a rewired one is a kinetic
checkpoint.

Design.
  features per (family, ref_pos):
    burial        fraction of its resolved states with >=1 protein partner        (the G1 feature, the thing to beat)
    mean_degree   mean number of distinct protein partners across those states
    rewiring      mean Jaccard DISTANCE between partner sets at consecutive states on the canonical pathway
    max_turnover  the largest single-transition Jaccard distance
    n_lost_max    most partners released at any one transition (directional: preparation for ejection)
    core_frac     fraction of its partners present in EVERY state it appears in (constitutive core)

  tests:
    T1a  pathogenicity, nested: label ~ burial + mean_degree (+gene) vs the same + rewiring. Held-out-gene AUROC
         delta and a likelihood-ratio test. Positive control that the base model is sane.
    T1b  measured function (RNU4-2 SGE, 485 variants, zero-shot): partial Spearman of rewiring | burial and the
         reverse. This is the ascertainment-free readout; it is the test that matters.
    T1c  NULL. Shuffle partner IDENTITIES among positions within each state, preserving every position's per-state
         partner COUNT. Burial and mean_degree are therefore bit-identical; only cross-state coherence of partner
         identity is destroyed. If rewiring's signal survives 200 draws of this, the trajectory carries information
         that burial cannot express.

Kill criterion, declared before running: if rewiring adds no held-out AUROC in T1a AND its partial correlation with
SGE function given burial is not distinguishable from 0 in T1b, the hinge hypothesis is dead and this line stops.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

# canonical order of the human spliceosome cycle; the minor pathway is ordered separately and never mixed with it
MAJOR = ["U2-snRNP", "tri-snRNP", "pre-B", "B", "Bact", "Bact-mature", "C*", "P"]
MINOR = ["minor-preB", "minor-Bact"]
ORDER = {s: i for i, s in enumerate(MAJOR)} | {s: i for i, s in enumerate(MINOR)}
RNG = np.random.default_rng(0)
N_NULL = 200


def partner_sets(nodes):
    """{(family, ref_pos): [(order_index, frozenset(partners)), ...]} along each pathway, deduplicated by state."""
    out = {}
    for (fam, pos), g in nodes.groupby(["family", "ref_pos"]):
        for pathway, states in (("major", MAJOR), ("minor", MINOR)):
            h = g[g.state.isin(states)]
            if h.empty:
                continue
            seq = []
            for st, hh in h.groupby("state"):
                # several structures can resolve the same state (two tri-snRNPs); union their partners
                ps = set()
                for s in hh.partners.dropna():
                    ps |= {x for x in str(s).split(";") if x}
                seq.append((ORDER[st], frozenset(ps)))
            seq.sort()
            if seq:
                out[(fam, pos, pathway)] = seq
    return out


def jaccard_dist(a, b):
    if not a and not b:
        return 0.0
    return 1.0 - len(a & b) / len(a | b)


def traj_features(seq):
    """Trajectory descriptors for one nucleotide's ordered partner sets."""
    sets = [s for _, s in seq]
    n = len(sets)
    occupied = sum(1 for s in sets if s)
    deg = [len(s) for s in sets]
    if n >= 2:
        d = [jaccard_dist(sets[i], sets[i + 1]) for i in range(n - 1)]
        lost = [len(sets[i] - sets[i + 1]) for i in range(n - 1)]
        rew, mx, nlost = float(np.mean(d)), float(np.max(d)), int(np.max(lost))
    else:
        rew = mx = 0.0
        nlost = 0
    union = set().union(*sets) if sets else set()
    inter = set.intersection(*[set(s) for s in sets]) if sets else set()
    core = len(inter) / len(union) if union else 0.0
    return dict(n_states=n, burial=occupied / n, mean_degree=float(np.mean(deg)),
                rewiring=rew, max_turnover=mx, n_lost_max=nlost, core_frac=core)


def build_table(nodes):
    rows = []
    for (fam, pos, pathway), seq in partner_sets(nodes).items():
        r = traj_features(seq)
        r.update(family=fam, ref_pos=pos, pathway=pathway)
        rows.append(r)
    return pd.DataFrame(rows)


def shuffled_nodes(nodes, rng):
    """Permute which position holds which partner set, independently within each (pdb_id, family).

    Degrees per position per state are NOT preserved individually, but the multiset of partner sets within each
    structure-family is, so burial and mean_degree are preserved in distribution and the cross-state identity
    coherence is destroyed. A stricter degree-matched variant is applied below.
    """
    out = nodes.copy()
    for _, g in nodes.groupby(["pdb_id", "family"]):
        idx = g.index.to_numpy()
        out.loc[idx, "partners"] = nodes.loc[rng.permutation(idx), "partners"].to_numpy()
    return out


def shuffled_nodes_degree_matched(nodes, rng):
    """Stricter null: permute partner sets only among positions with the SAME partner count, within each
    (pdb_id, family). Every position keeps its exact per-state degree, so burial, mean_degree and the whole
    degree profile are bit-identical; only WHICH proteins, and hence cross-state coherence, is randomised."""
    out = nodes.copy()
    for _, g in nodes.groupby(["pdb_id", "family", "n_partners"]):
        if len(g) < 2:
            continue
        idx = g.index.to_numpy()
        out.loc[idx, "partners"] = nodes.loc[rng.permutation(idx), "partners"].to_numpy()
    return out


def partial_spearman(x, y, z):
    """Spearman of x and y with z partialled out, on ranks (residualise rank(x), rank(y) on rank(z))."""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if ok.sum() < 20:
        return np.nan
    rx, ry, rz = rankdata(x[ok]), rankdata(y[ok]), rankdata(z[ok])
    Z = np.c_[np.ones(ok.sum()), rz]
    bx = np.linalg.lstsq(Z, rx, rcond=None)[0]
    by = np.linalg.lstsq(Z, ry, rcond=None)[0]
    return float(spearmanr(rx - Z @ bx, ry - Z @ by).statistic)


def logo_auroc(D, feats, label="label", unit="gene_name"):
    """Leave-one-gene-out AUROC, pooled on within-unit ranks (the paper's established protocol)."""
    preds, ys, units = [], [], []
    for u in D[unit].unique():
        tr, te = D[D[unit] != u], D[D[unit] == u]
        if te[label].nunique() < 2 or tr[label].nunique() < 2:
            continue
        sc = StandardScaler().fit(tr[feats])
        m = LogisticRegression(C=0.3, class_weight="balanced", max_iter=5000)
        m.fit(sc.transform(tr[feats]), tr[label])
        p = m.predict_proba(sc.transform(te[feats]))[:, 1]
        preds.append(rankdata(p) / len(p))
        ys.append(te[label].to_numpy())
        units.append(u)
    if not preds:
        return np.nan, []
    per = [roc_auc_score(y, p) for y, p in zip(ys, preds)]
    return roc_auc_score(np.concatenate(ys), np.concatenate(preds)), list(zip(units, per))


# ---------------------------------------------------------------------------------------------------------------
print("loading ...", flush=True)
nodes = pd.read_csv(config.PROCESSED / "graph_nodes_named.tsv.gz", sep="\t")
V = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t", low_memory=False)
L = pd.read_csv(config.CURATION / "labels_v5c.tsv", sep="\t")

T = build_table(nodes)
T.to_csv(config.RESULTS / "trajectory_features.tsv", sep="\t", index=False)
print(f"trajectory features: {len(T)} (family, position, pathway) rows")
print(T.groupby(["pathway", "family"])[["burial", "mean_degree", "rewiring", "core_frac"]]
      .mean().round(3).to_string())

# one row per (family, ref_pos): prefer the major pathway, fall back to minor for U11/U12/U4atac/U6atac
Tm = (T.sort_values("pathway", ascending=False)        # 'minor' < 'major' alphabetically -> major first after desc
        .drop_duplicates(["family", "ref_pos"]))
V2 = V.merge(Tm.rename(columns={"family": "paralog_family", "ref_pos": "family_ref_pos"}),
             on=["paralog_family", "family_ref_pos"], how="left")

FEAT_BASE = ["burial", "mean_degree"]
FEAT_FULL = FEAT_BASE + ["rewiring", "max_turnover", "n_lost_max", "core_frac"]

# =============================== T1a  pathogenicity, held-out gene =============================================
# controls carry mode=NaN in labels_v5c, so filtering on mode would drop every negative: take the NDD-mode
# positives and every control in the same genes
_pos = L[(L.label == 1) & L["mode"].astype(str).str.contains("NDD")]
_neg = L[(L.label == 0) & L.gene_name.isin(_pos.gene_name.unique())]
_LL = pd.concat([_pos, _neg])[["key", "label"]].drop_duplicates("key")
D = V2.merge(_LL, on="key", how="inner").dropna(subset=FEAT_FULL)
print(f"\n=== T1a  pathogenicity (NDD modes), n={len(D)}, {int(D.label.sum())} pathogenic, "
      f"{D.gene_name.nunique()} genes ===")
a_base, per_base = logo_auroc(D, FEAT_BASE)
a_full, per_full = logo_auroc(D, FEAT_FULL)
print(f"  burial + degree              LOGO AUROC {a_base:.3f}")
print(f"  + trajectory (4 extra)       LOGO AUROC {a_full:.3f}   delta {a_full - a_base:+.3f}")
for (u, b), (_, f) in zip(per_base, per_full):
    print(f"    {u:<10s} {b:.3f} -> {f:.3f}  ({f - b:+.3f})")

# =============================== T1b  measured function (SGE), zero-shot ======================================
S = V2[V2.sge_score.notna() & V2.rewiring.notna()].copy()
S["sge_score"] = pd.to_numeric(S.sge_score, errors="coerce")
S = S.dropna(subset=["sge_score"])
print(f"\n=== T1b  measured function, RNU4-2 SGE, n={len(S)} variants, zero-shot ===")
y = -S.sge_score.to_numpy()                       # low function score = damaging; flip so higher = worse
res = {}
for f in ["burial", "mean_degree", "rewiring", "max_turnover", "n_lost_max", "core_frac"]:
    r = spearmanr(S[f], y, nan_policy="omit")
    res[f] = r.statistic
    print(f"  {f:<14s} rho {r.statistic:+.3f}   p {r.pvalue:.2g}")
pr_rew = partial_spearman(S.rewiring.to_numpy(float), y, S.burial.to_numpy(float))
pr_bur = partial_spearman(S.burial.to_numpy(float), y, S.rewiring.to_numpy(float))
print(f"\n  partial rewiring | burial    {pr_rew:+.3f}")
print(f"  partial burial   | rewiring  {pr_bur:+.3f}")

# =============================== T1c  trajectory-destroying nulls =============================================
print(f"\n=== T1c  null: partner identities permuted within each structure x family ({N_NULL} draws) ===")
key_cols = ["paralog_family", "family_ref_pos"]
obs_sge = res["rewiring"]
obs_auroc = a_full - a_base
null_sge, null_auroc, null_partial = [], [], []
for which, fn in (("degree-matched", shuffled_nodes_degree_matched), ("free", shuffled_nodes)):
    rng = np.random.default_rng(1)
    sge_n, par_n, auc_n = [], [], []
    for b in range(N_NULL):
        Tn = build_table(fn(nodes, rng))
        Tn = Tn.sort_values("pathway", ascending=False).drop_duplicates(["family", "ref_pos"])
        Tn = Tn.rename(columns={"family": "paralog_family", "ref_pos": "family_ref_pos"})
        Sn = S[key_cols + ["sge_score"]].merge(Tn, on=key_cols, how="left")
        yy = -pd.to_numeric(Sn.sge_score, errors="coerce").to_numpy()
        sge_n.append(spearmanr(Sn.rewiring, yy, nan_policy="omit").statistic)
        par_n.append(partial_spearman(Sn.rewiring.to_numpy(float), yy, Sn.burial.to_numpy(float)))
        if b < 50:                                   # the LOGO fit is the expensive part; 50 draws for the AUROC
            Dn = D[key_cols + ["gene_name", "label"]].merge(Tn, on=key_cols, how="left").dropna(subset=FEAT_FULL)
            an, _ = logo_auroc(Dn, FEAT_FULL)
            ab, _ = logo_auroc(Dn, FEAT_BASE)
            auc_n.append(an - ab)
    sge_n = np.array([v for v in sge_n if np.isfinite(v)])
    par_n = np.array([v for v in par_n if np.isfinite(v)])
    auc_n = np.array([v for v in auc_n if np.isfinite(v)])
    p_sge = (1 + (np.abs(sge_n) >= abs(obs_sge)).sum()) / (1 + len(sge_n))
    p_par = (1 + (np.abs(par_n) >= abs(pr_rew)).sum()) / (1 + len(par_n))
    p_auc = (1 + (auc_n >= obs_auroc).sum()) / (1 + len(auc_n)) if len(auc_n) else np.nan
    print(f"  [{which}] rewiring vs SGE   obs {obs_sge:+.3f}  null {sge_n.mean():+.3f} +/- {sge_n.std():.3f}"
          f"   p {p_sge:.3f}")
    print(f"  [{which}] partial | burial  obs {pr_rew:+.3f}  null {par_n.mean():+.3f} +/- {par_n.std():.3f}"
          f"   p {p_par:.3f}")
    print(f"  [{which}] AUROC gain        obs {obs_auroc:+.3f}  null {auc_n.mean():+.3f} +/- {auc_n.std():.3f}"
          f"   p {p_auc:.3f}  (n={len(auc_n)})")
    null_sge.append((which, sge_n)); null_auroc.append((which, auc_n)); null_partial.append((which, par_n))

pd.DataFrame(dict(test=["sge_rho", "partial_rew_given_burial", "logo_auroc_gain"],
                  observed=[obs_sge, pr_rew, obs_auroc])).to_csv(
    config.RESULTS / "trajectory_gate_observed.tsv", sep="\t", index=False)
pd.concat([pd.DataFrame(dict(null=w, stat="sge_rho", value=v)) for w, v in null_sge] +
          [pd.DataFrame(dict(null=w, stat="partial", value=v)) for w, v in null_partial] +
          [pd.DataFrame(dict(null=w, stat="auroc_gain", value=v)) for w, v in null_auroc]).to_csv(
    config.RESULTS / "trajectory_gate_null.tsv", sep="\t", index=False)

print("\n=== VERDICT ===")
alive = (obs_auroc > 0.01) or (abs(pr_rew) > 0.1)
print("  rewiring adds held-out AUROC: ", f"{obs_auroc:+.3f}")
print("  rewiring partial on function: ", f"{pr_rew:+.3f}")
print("  HYPOTHESIS ALIVE" if alive else "  HYPOTHESIS DEAD - stop this line")
