"""Train and evaluate the snRNA-VEP spliceosome network (snrna_vep/net.py).

E1  SGE zero-shot: no RNU4-2 data in training (gene-out, and stricter family-out); Spearman with -SGE.
E2  Held-out pathogenicity: eval units RNU4-2, RNU2-2P, RNU5B-1 and the RNU6 RP genes (as one unit);
    gene-out (other paralog-family members may train) and family-out (whole family removed).
All models and baselines are scored on the identical test variants (those with a structure node).
5-seed ensembles; hyperparameters fixed in advance (hidden 32, 2 layers, 300 epochs, AdamW 3e-3, wd 1e-3).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import bootstrap_auroc, metrics, paired_bootstrap_delta
from snrna_vep.model import FEATURE_GROUPS, fit_predict
from snrna_vep.net import MODES, build_graph, train_predict, variant_tensor

FEATURE_GROUPS["conservation"] = ["phylop447"]
SEEDS = [0, 1, 2, 3, 4]
torch.set_num_threads(4)

feat = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t", low_memory=False)
pat = pd.read_csv(config.CURATION / "patient_variants_curated.tsv", sep="\t")
ctl = pd.read_csv(config.CURATION / "controls_curated.tsv", sep="\t")
gn = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")
NODE_FEATS = "graph_node_features_agnostic.tsv.gz" if "--per-state" not in sys.argv else "graph_node_features.tsv.gz"
nodes_f = pd.read_csv(config.PROCESSED / NODE_FEATS, sep="\t")
edges = pd.read_csv(config.PROCESSED / "graph_edges.tsv.gz", sep="\t")
idx, X, A = build_graph(nodes_f, edges)

# ---- labels v2 (scripts/05e_labels_v2.py): curated supplements + ClinVar 2026-10, mode-aware controls
L2 = pd.read_csv(config.CURATION / config.LABELS, sep="\t")
LAB = L2.groupby("key").label.max()                  # a key labelled pathogenic anywhere is pathogenic
mode_of = L2[L2.label == 1].drop_duplicates("key").set_index("key")["mode"]
pat = L2[L2.label == 1]


def label(g, k):
    return float(LAB[k]) if k in LAB.index else np.nan


d = feat[feat.paralog_family.notna() & feat.family_ref_pos.notna()].copy()
d["node_name"] = d.paralog_family + ":" + d.family_ref_pos.astype(float).astype(int).astype(str)
d = d[d.node_name.isin(idx)].drop_duplicates("key").reset_index(drop=True)
# Design correction (before any tuning): population variants in non-disease paralog copies (RNU4-1, RNU5D/E/F-1,
# RNU11, ...) sit at the same structural positions as pathogenic variants in the expressed disease copy and are
# tolerated there, so labelling them benign teaches "structurally critical -> benign". Training and evaluation use
# only genes with patient variants (within-gene contrasts); applied identically to every model.
DISEASE_GENES = sorted(set(pat.gene_name))
d = d[d.gene_name.isin(DISEASE_GENES)].reset_index(drop=True)
d["y_path"] = [label(g, k) for g, k in zip(d.gene_name, d.key)]
d["y_mode"] = [MODES.index(mode_of[k]) if k in mode_of.index and mode_of[k] in MODES else -1 for k in d.key]
sge = d.sge_score if "sge_score" in d else pd.Series(np.nan, index=d.index)
d["y_sge"] = -(sge - sge.mean()) / sge.std()
d["phylop447"] = d.phylop447 / 10.0
V = variant_tensor(d)
V = (V - V.mean(0)) / (V.std(0) + 1e-6)
FAMILY = {"RNU4-2": "U4", "RNU2-2P": "U2", "RNU5B-1": "U5", "RNU6": "U6", "RNU4ATAC": "U4", "RNU12": "U2", "RNU6ATAC": "U6"}

print(f"graph: {len(idx)} nodes; variants with a structure node: {len(d)}; labelled: {int(d.y_path.notna().sum())} "
      f"(pathogenic {int((d.y_path == 1).sum())}); SGE: {int(d.y_sge.notna().sum())}")


def pack(sub, use_sge=True):
    t = lambda a, dt=torch.float32: torch.tensor(np.asarray(a), dtype=dt)
    return dict(node=t([idx[n] for n in sub.node_name], torch.long), vfeat=t(V[sub.index]),
                y_path=t(sub.y_path.values), y_mode=t(sub.y_mode.values, torch.long),
                y_sge=t(sub.y_sge.values if use_sge else np.full(len(sub), np.nan)))


ARGS = dict(n_node_feat=X.shape[1], n_var_feat=V.shape[1], hidden=32, layers=2, dropout=0.2)


def ensemble(train, test, use_sge=True):
    ps, ss = [], []
    for s in SEEDS:
        p, _, sg = train_predict(ARGS, X, A, pack(train, use_sge), pack(test), seed=s)
        ps.append(p); ss.append(sg)
    return np.mean(ps, 0), np.mean(ss, 0)


def unit_mask(df, unit):
    return df.gene_name.str.startswith("RNU6-") if unit == "RNU6" else (df.gene_name == unit)


rows, pooled = [], []
for mode in ["gene-out", "family-out"]:
    for unit in ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6", "RNU4ATAC", "RNU12", "RNU6ATAC"]:
        test = d[unit_mask(d, unit) & d.y_path.notna()]
        if (test.y_path == 1).sum() < 3:
            continue
        drop = unit_mask(d, unit) if mode == "gene-out" else (d.paralog_family == FAMILY[unit])
        train = d[~drop]
        use_sge = unit != "RNU4-2" and not (mode == "family-out" and FAMILY[unit] == "U4")
        p, _ = ensemble(train[train.y_path.notna() | train.y_sge.notna()], test, use_sge=use_sge)
        tr_l = train[train.y_path.notna()]
        lc, _, _ = fit_predict(tr_l.assign(label=tr_l.y_path.astype(int)), test.assign(label=0), ["contacts"], "logistic")
        t = test.assign(net=p, logistic_contacts=lc)
        pooled.append(t.assign(eval_mode=mode, unit=unit))
        for s in ["net", "logistic_contacts", "frac_states_protein_contact", "max_protein_res", "cadd_phred", "phylop447"]:
            m = metrics(t.y_path, t[s])
            rows.append(dict(eval=mode, unit=unit, score=s, **m))

res = pd.DataFrame(rows)
res.to_csv(config.result("net_heldout.tsv"), sep="\t", index=False)
print("\n=== Held-out pathogenicity, AUROC (identical test variants) ===")
print(res.pivot_table(index="score", columns=["eval", "unit"], values="auroc").round(3).to_string())
P = pd.concat(pooled)
for mode in ["gene-out", "family-out"]:
    t = P[P.eval_mode == mode]
    lo, hi = bootstrap_auroc(t.y_path, t.net, n=1000)
    print(f"\n[{mode}] pooled net AUROC {metrics(t.y_path, t.net)['auroc']:.3f} [{lo:.3f}, {hi:.3f}]  n {len(t)} pathogenic {int(t.y_path.sum())}")
    for b in ["frac_states_protein_contact", "logistic_contacts", "phylop447"]:
        dl, l, h, pv = paired_bootstrap_delta(t.y_path.values, t.net.values, t[b].values, n=1000)
        print(f"   net - {b:28s}: {dl:+.3f} [{l:+.3f}, {h:+.3f}] p {pv:.3f}")
    s = t[t.cadd_phred.notna()]
    dl, l, h, pv = paired_bootstrap_delta(s.y_path.values, s.net.values, s.cadd_phred.values, n=1000)
    print(f"   net - CADD (SNVs, n {len(s)}, pathogenic {int(s.y_path.sum())}): {dl:+.3f} [{l:+.3f}, {h:+.3f}] p {pv:.3f}")
P.to_csv(config.result("net_heldout_predictions.tsv.gz"), sep="\t", index=False)

# ---- E1: SGE zero-shot
print("\n=== RNU4-2 SGE zero-shot (Spearman with -SGE) ===")
test = d[(d.gene_name == "RNU4-2") & d.y_sge.notna()]
out = []
for mode, drop in [("gene-out", d.gene_name == "RNU4-2"), ("family-out", d.paralog_family == "U4")]:
    train = d[~drop & (d.y_path.notna())]
    p, _ = ensemble(train, test, use_sge=False)
    rho = spearmanr(p, test.y_sge)[0]
    out.append(dict(eval=mode, score="net", spearman=rho, n=len(test)))
for s in ["frac_states_protein_contact", "max_protein_res", "cadd_phred", "phylop447"]:
    x = test[[s, "y_sge"]].dropna()
    out.append(dict(eval="-", score=s, spearman=spearmanr(x[s], x.y_sge)[0], n=len(x)))
z = pd.DataFrame(out)
z.to_csv(config.result("net_sge_zeroshot.tsv"), sep="\t", index=False)
print(z.round(3).to_string(index=False))
