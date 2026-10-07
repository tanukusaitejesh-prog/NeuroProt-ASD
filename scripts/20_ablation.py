"""Ablation: does combining 11 cryo-EM structures beat the best single structure? Is the contact graph needed?

All structural features are recomputed from data/processed/graph_nodes.tsv.gz / graph_edges.tsv.gz for a chosen set of
structures (positions not resolved in a structure count as 'no contact'):
  contact fraction (structures with protein contact / structures resolving the position), max protein atoms,
  snRNA-contact fraction, graph-smoothed contact fraction (graph from RNA-RNA contacts in the selected structures).
Same model (v2 design, phyloP447 included), same held-out-unit protocol, labels as in config.LABELS.
Configurations: ALL (multi-state, recomputed), ALL-nograph (graph feature replaced by raw contact fraction),
each single structure, ORACLE single (best structure chosen on the test units: optimistic for single-state),
NESTED single (structure chosen by inner leave-one-unit-out on training units only).
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
d, rows_for_training, design, model, R, nodes_f = G["d"], G["rows_for_training"], G["design"], G["model"], G["R"], G["nodes_f"]
nodes = pd.read_csv(config.PROCESSED / "graph_nodes.tsv.gz", sep="\t")
edges = pd.read_csv(config.PROCESSED / "graph_edges.tsv.gz", sep="\t")
nodes["node"] = nodes.family + ":" + nodes.ref_pos.astype(int).astype(str)
PDBS = sorted(nodes.pdb_id.unique())
UNITS = ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6", "RNU4ATAC", "RNU12"]


def set_features(pdbs, graph=True):
    n = nodes[nodes.pdb_id.isin(pdbs)]
    g = n.groupby("node").agg(frac=("n_protein_atoms", lambda x: (x > 0).mean()), mx=("n_protein_atoms", "max"),
                              rna=("n_rna_nb", lambda x: (x > 0).mean()))
    e = edges[edges.pdb_id.isin(pdbs)]
    nf = nodes_f[["family", "ref_pos"]].copy()
    nf["node"] = nf.family + ":" + nf.ref_pos.astype(int).astype(str)
    nf["frac_prot"] = nf.node.map(g.frac).fillna(0.0)
    idx, _, A = build_graph(nf[["family", "ref_pos", "frac_prot"]], e)
    fp = torch.tensor(nf.set_index("node").reindex(list(idx)).frac_prot.to_numpy(), dtype=torch.float32).unsqueeze(1)
    nb = torch.sparse.mm(A, fp).squeeze(1).numpy()
    nbm = dict(zip(idx, nb))
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


def evaluate(test_units=UNITS, train_exclude=()):
    """held-out AUROC per (unit, ctx) for the current feature columns"""
    out = {}
    for unit in test_units:
        tr = rows_for_training(d[(d.unit != unit) & ~d.unit.isin(train_exclude)])
        m = model().fit(design(d, tr.i, tr.ctx), tr.y)
        for ctx in ["dominant", "recessive"]:
            t = ctx_rows(unit, ctx)
            if t is not None:
                out[(unit, ctx)] = metrics(t.y_path, m.predict_proba(design(d, t.index, [ctx] * len(t)))[:, 1])["auroc"]
    return out


configs = {"ALL (11 structures)": (PDBS, True), "ALL, no graph": (PDBS, False)}
configs.update({p: ([p], True) for p in PDBS})
res = {}
for name, (pdbs, graph) in configs.items():
    set_features(pdbs, graph)
    res[name] = evaluate()
    print(f"  {name:22s} " + " ".join(f"{k[0]}/{k[1][:3]} {v:.3f}" for k, v in res[name].items()), flush=True)
R_ = pd.DataFrame(res).T
state = nodes.drop_duplicates("pdb_id").set_index("pdb_id").state
keys = list(R_.columns)
summ = pd.DataFrame({ctx: R_[[k for k in keys if k[1] == ctx]].mean(1) for ctx in ["dominant", "recessive"]})
summ["state"] = [state.get(i, "") for i in summ.index]
singles = summ.loc[PDBS]
oracle = {ctx: singles[ctx].max() for ctx in ["dominant", "recessive"]}

# nested single-structure choice: for each test unit, pick the structure with best mean inner AUROC on the other units
nested = {}
for unit in UNITS:
    for ctx in ["dominant", "recessive"]:
        if (unit, ctx) not in keys:
            continue
        inner = {}
        for p in PDBS:
            vals = [R_.loc[p, k] for k in keys if k[1] == ctx and k[0] != unit]
            inner[p] = np.mean(vals)
        best = max(inner, key=inner.get)
        set_features([best], True)
        nested[(unit, ctx)] = evaluate(test_units=[unit])[(unit, ctx)]
nest = {ctx: np.mean([v for k, v in nested.items() if k[1] == ctx]) for ctx in ["dominant", "recessive"]}
summ.loc["NESTED single structure"] = [nest["dominant"], nest["recessive"], ""]
summ.loc["ORACLE best single structure"] = [oracle["dominant"], oracle["recessive"], ""]
R_.to_csv(config.result("ablation_per_unit.tsv"), sep="\t")
summ.to_csv(config.result("ablation_summary.tsv"), sep="\t")
print("\n=== mean held-out AUROC ===")
print(summ.round(3).sort_values("dominant", ascending=False).to_string())
