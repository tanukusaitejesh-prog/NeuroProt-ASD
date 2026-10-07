"""snRNA-VEP network: message passing over the multi-state spliceosome contact graph, multi-task heads.

Nodes   : snRNA nucleotide positions in paralog-family reference coordinates (U1, U2, U4, U5, U6; minor
          spliceosome snRNAs are projected onto their major analogs).
Node x  : per-state contact features + partner identity (snrna_vep.graph.node_feature_matrix).
Edges   : RNA-RNA contacts in any cryo-EM state (weighted by number of states) + backbone neighbours.
Variant : node embedding of the affected position + the change itself (ref/alt base, SNV/ins/del,
          transition, conservation) -> heads for pathogenicity, mechanism and SGE function score.
"""
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

MODES = ["AD-NDD", "AR-NDD", "AD-RP"]
VAR_COLS = ["is_snv", "is_ins", "is_del", "transition", "phylop447", "rel_pos", "local_gc"]


def build_graph(node_feat, edges):
    """Return node index, standardized feature tensor and a normalized sparse adjacency."""
    nf = node_feat.copy()
    nf["node"] = nf.family + ":" + nf.ref_pos.astype(int).astype(str)
    idx = {n: i for i, n in enumerate(nf.node)}
    X = nf.drop(columns=["family", "ref_pos", "node"]).astype(float).to_numpy()
    X = (X - X.mean(0)) / (X.std(0) + 1e-6)
    e = edges.assign(a=edges.family_a + ":" + edges.pos_a.astype(int).astype(str),
                     b=edges.family_b + ":" + edges.pos_b.astype(int).astype(str))
    w = e.groupby(["a", "b"]).state.nunique().reset_index()
    src, dst, wt = [], [], []
    for a, b, k in zip(w.a, w.b, w.state):
        if a in idx and b in idx:
            src.append(idx[a]); dst.append(idx[b]); wt.append(float(k))
    # backbone neighbours within a family
    for n, i in idx.items():
        fam, pos = n.split(":")
        nb = f"{fam}:{int(pos) + 1}"
        if nb in idx:
            for s, d in ((i, idx[nb]), (idx[nb], i)):
                src.append(s); dst.append(d); wt.append(1.0)
    n = len(idx)
    src += list(range(n)); dst += list(range(n)); wt += [1.0] * n          # self loops
    A = torch.sparse_coo_tensor(torch.tensor([dst, src]), torch.tensor(wt, dtype=torch.float32), (n, n)).coalesce()
    deg = torch.sparse.sum(A, 1).to_dense()
    vals = A.values() / deg[A.indices()[0]]
    A = torch.sparse_coo_tensor(A.indices(), vals, (n, n)).coalesce()
    return idx, torch.tensor(X, dtype=torch.float32), A


class SpliceosomeNet(nn.Module):
    def __init__(self, n_node_feat, n_var_feat, hidden=32, layers=2, dropout=0.2):
        super().__init__()
        self.inp = nn.Sequential(nn.Linear(n_node_feat, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.mp = nn.ModuleList([nn.Linear(2 * hidden, hidden) for _ in range(layers)])
        self.drop = nn.Dropout(dropout)
        self.var = nn.Sequential(nn.Linear(hidden + n_var_feat, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.path = nn.Linear(hidden, 1)
        self.mode = nn.Linear(hidden, len(MODES))
        self.sge = nn.Linear(hidden, 1)

    def node_embed(self, X, A):
        h = self.inp(X)
        for lin in self.mp:
            h = self.drop(torch.relu(lin(torch.cat([h, torch.sparse.mm(A, h)], 1)))) + h
        return h

    def forward(self, X, A, node_idx, vfeat):
        h = self.node_embed(X, A)[node_idx]
        z = self.var(torch.cat([h, vfeat], 1))
        return self.path(z).squeeze(1), self.mode(z), self.sge(z).squeeze(1)


def variant_tensor(df, base_map=("A", "C", "G", "T")):
    v = df[VAR_COLS].astype(float).fillna(0.0).to_numpy()
    ref = np.zeros((len(df), 4)); alt = np.zeros((len(df), 4))
    for i, h in enumerate(df.get("hgvs", pd.Series([""] * len(df))).fillna("").astype(str)):
        if ">" in h:
            r, a = h[-3], h[-1]
            if r in base_map: ref[i, base_map.index(r)] = 1
            if a in base_map: alt[i, base_map.index(a)] = 1
    return np.concatenate([v, ref, alt], 1)


def train_predict(model_args, X, A, train, test, epochs=300, lr=3e-3, wd=1e-3, sge_weight=0.5, seed=0):
    """train/test: dicts with node (LongTensor), vfeat, y_path (float, nan=unlabelled), y_mode (long, -1 = none),
    y_sge (float, nan = none). Returns test pathogenicity probabilities, mode probabilities and SGE predictions."""
    torch.manual_seed(seed); np.random.seed(seed)
    m = SpliceosomeNet(**model_args)
    opt = torch.optim.AdamW(m.parameters(), lr=lr, weight_decay=wd)
    yp, ym, ys = train["y_path"], train["y_mode"], train["y_sge"]
    has_p, has_m, has_s = ~torch.isnan(yp), ym >= 0, ~torch.isnan(ys)
    pos_w = torch.tensor(float((yp[has_p] == 0).sum() / max(1, (yp[has_p] == 1).sum())))
    bce = nn.BCEWithLogitsLoss(pos_weight=pos_w)
    for _ in range(epochs):
        m.train(); opt.zero_grad()
        lp, lm, ls = m(X, A, train["node"], train["vfeat"])
        loss = torch.tensor(0.0)
        if has_p.any():
            loss = loss + bce(lp[has_p], yp[has_p])
        if has_m.any():
            loss = loss + 0.5 * nn.functional.cross_entropy(lm[has_m], ym[has_m])
        if has_s.any():
            loss = loss + sge_weight * nn.functional.mse_loss(ls[has_s], ys[has_s])
        loss.backward(); opt.step()
    m.eval()
    with torch.no_grad():
        lp, lm, ls = m(X, A, test["node"], test["vfeat"])
    return torch.sigmoid(lp).numpy(), torch.softmax(lm, 1).numpy(), ls.numpy()
