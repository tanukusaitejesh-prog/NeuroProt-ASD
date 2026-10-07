"""Ablation: does combining 11 cryo-EM structures beat the best single structure? Is the contact graph needed?

Structural features are recomputed with the SAME definitions as the main pipeline (snrna_vep.contacts.aggregate_states,
from data/processed/contacts_per_structure.tsv) for a chosen set of structures: fraction of splicing STATES with protein
contact, max protein residues, fraction of states with snRNA contact; graph-smoothed contact fraction uses the RNA-RNA
contact graph of the selected structures only. Positions not resolved in the selection count as 'no contact'.
'ALL' must reproduce the main model (sanity check printed).
Same model (v2 design, phyloP447 included), same held-out-unit protocol, labels as in config.LABELS.
Configurations: ALL (multi-state, recomputed), ALL-nograph (graph feature replaced by raw contact fraction),
each single structure, ORACLE single (best structure chosen on the test units: optimistic for single-state),
NESTED single (structure chosen by inner leave-one-unit-out over the other units' held-out AUROCs).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import metrics
from snrna_vep.net import build_graph

src = open(Path(__file__).resolve().parent / "11_v2_model.py").read().split("results, pooled = [], []")[0]
G = {"__file__": str(Path(__file__).resolve().parent / "11_v2_model.py"), "__name__": "abl"}
exec(compile(src, "11_v2_model", "exec"), G)
d, rows_for_training, design, model, nodes_f = G["d"], G["rows_for_training"], G["design"], G["model"], G["nodes_f"]
orig = {c: d[c + "_r"].copy() for c in ["frac_states_protein_contact", "max_protein_res", "frac_states_snrna_contact", "nb_frac_prot"]}
per = pd.read_csv(config.PROCESSED / "contacts_per_structure.tsv", sep="\t").dropna(subset=["paralog_family", "family_ref_pos"])
per["node"] = per.paralog_family + ":" + per.family_ref_pos.astype(int).astype(str)
nodes = pd.read_csv(config.PROCESSED / "graph_nodes.tsv.gz", sep="\t")
edges = pd.read_csv(config.PROCESSED / "graph_edges.tsv.gz", sep="\t")
nodes["node"] = nodes.family + ":" + nodes.ref_pos.astype(int).astype(str)
PDBS = sorted(nodes.pdb_id.unique())
UNITS = ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6", "RNU4ATAC", "RNU12"]


def set_features(pdbs, graph=True):
    n = per[per.pdb_id.isin(pdbs)]
    st = n.assign(p=n.n_protein_res > 0, r=n.n_snrna_contacts > 0).groupby(["node", "state"]).agg(p=("p", "any"), r=("r", "any"))
    g = pd.DataFrame({"frac": st.p.groupby(level=0).mean(), "rna": st.r.groupby(level=0).mean(),
                      "mx": n.groupby("node").n_protein_res.max()})
    e = edges[edges.pdb_id.isin(pdbs)]
    nf = nodes_f[["family", "ref_pos"]].copy()
    nf["node"] = nf.family + ":" + nf.ref_pos.astype(int).astype(str)
    nf["frac_prot"] = nf.node.map(g.frac).fillna(0.0)
    idx, _, A = build_graph(nf[["family", "ref_pos", "frac_prot"]], e)
    fp = torch.tensor(nf.set_index("node").reindex(list(idx)).frac_prot.to_numpy(), dtype=torch.float32).unsqueeze(1)
    nbm = dict(zip(idx, torch.sparse.mm(A, fp).squeeze(1).numpy()))
    raw = {"frac_states_protein_contact": d.node_name.map(g.frac).fillna(0.0),
           "max_protein_res": d.node_name.map(g.mx).fillna(0.0),
           "frac_states_snrna_contact": d.node_name.map(g.rna).fillna(0.0)}
    raw["nb_frac_prot"] = d.node_name.map(nbm).fillna(0.0) if graph else raw["frac_states_protein_contact"]
    for c, v in raw.items():
        d[c + "_r"] = v.groupby(d.gene_name).rank(pct=True)


def ctx_rows(unit, ctx):
    pref = "AD-" if ctx == "dominant" else "AR-"
    t = d[(d.unit == unit) & ((d.y_path == 0) | ((d.y_path == 1) & d["mode"].str.contains(pref)))]
    return t if (t.y_path == 1).sum() >= 3 and (t.y_path == 0).sum() >= 3 else None


def evaluate(test_units=UNITS):
    out = {}
    for unit in test_units:
        tr = rows_for_training(d[d.unit != unit])
        m = model().fit(design(d, tr.i, tr.ctx), tr.y)
        for ctx in ["dominant", "recessive"]:
            t = ctx_rows(unit, ctx)
            if t is not None:
                out[(unit, ctx)] = metrics(t.y_path, m.predict_proba(design(d, t.index, [ctx] * len(t)))[:, 1])["auroc"]
    return out


set_features(PDBS, True)
from scipy.stats import spearmanr
print("sanity, recomputed ALL vs main-pipeline features (Spearman):",
      {c: round(spearmanr(orig[c], d[c + "_r"], nan_policy="omit")[0], 3) for c in orig})
configs = {"ALL (11 structures)": (PDBS, True), "ALL, no graph": (PDBS, False)}
configs.update({p: ([p], True) for p in PDBS})
res = {}
for name, (pdbs, graph) in configs.items():
    set_features(pdbs, graph)
    res[name] = evaluate()
    print(f"  {name:22s} " + " ".join(f"{k[0]}/{k[1][:3]} {v:.3f}" for k, v in res[name].items()), flush=True)
T = pd.DataFrame(res).T
keys = list(T.columns)
state = nodes.drop_duplicates("pdb_id").set_index("pdb_id").state
summ = pd.DataFrame({ctx: T[[k for k in keys if k[1] == ctx]].mean(1) for ctx in ["dominant", "recessive"]})
summ["state"] = [state.get(i, "") for i in summ.index]
singles = summ.loc[PDBS]
nested = {}
for unit, ctx in keys:
    inner = {p: np.mean([T.loc[p, k] for k in keys if k[1] == ctx and k[0] != unit]) for p in PDBS}
    nested[(unit, ctx)] = T.loc[max(inner, key=inner.get), (unit, ctx)]
for ctx in ["dominant", "recessive"]:
    summ.loc["NESTED single structure", ctx] = np.mean([v for k, v in nested.items() if k[1] == ctx])
    summ.loc["ORACLE best single structure", ctx] = singles[ctx].max()
T.to_csv(config.result("ablation_per_unit.tsv"), sep="\t")
summ.to_csv(config.result("ablation_summary.tsv"), sep="\t")
print("\n=== mean held-out AUROC ===")
print(summ.round(3).sort_values("dominant", ascending=False).to_string())
