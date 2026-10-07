"""Precomputed snRNA-VEP scores for every enumerated variant in every spliceosomal snRNA gene with a structure node
(disease genes and candidate/paralog genes: RNU5A-1, RNU4-1, RNU6-7, RNU11, RNU5D/E/F-1, ...).

Final model = the system selected by nested cross-validation: v2 (trained on all labelled variants) for the dominant
mechanism; contacts+phyloP (untrained) for the recessive mechanism. Features are within-gene percentile ranks over all
possible variants of the gene, exactly as in training. Scores are also given as within-gene percentiles.
Usage: SNRNA_LABELS=labels_v5.tsv python scripts/19_score_resource.py  -> results/snrnavep_scores.tsv.gz
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

src = open(Path(__file__).resolve().parent / "11_v2_model.py").read().split("results, pooled = [], []")[0]
G = {"__file__": str(Path(__file__).resolve().parent / "11_v2_model.py"), "__name__": "res"}
exec(compile(src, "11_v2_model", "exec"), G)
d, idx, A, nodes_f = G["d"], G["idx"], G["A"], G["nodes_f"]
rows_for_training, design, model, FEATS, R = G["rows_for_training"], G["design"], G["model"], G["FEATS"], G["R"]
tr = rows_for_training(d)
m = model().fit(design(d, tr.i, tr.ctx), tr.y)

# all genes with structure nodes (same preprocessing as scripts/09 + 11, without the disease-gene restriction)
f = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t", low_memory=False)
f = f[f.paralog_family.notna() & f.family_ref_pos.notna()].copy()
f["node_name"] = f.paralog_family + ":" + f.family_ref_pos.astype(float).astype(int).astype(str)
f = f[f.node_name.isin(idx)].drop_duplicates("key").reset_index(drop=True)
f["phylop447"] = f.phylop447 / 10.0
fp = torch.tensor(nodes_f.set_index(nodes_f.family + ":" + nodes_f.ref_pos.astype(int).astype(str))
                  .reindex(list(idx)).frac_prot.fillna(0).to_numpy(), dtype=torch.float32).unsqueeze(1)
nb = torch.sparse.mm(A, fp).squeeze(1).numpy()
f["nb_frac_prot"] = nb[[idx[n] for n in f.node_name]]
for c in FEATS:
    f[c + "_r"] = f.groupby("gene_name")[c].rank(pct=True)
f["score_dominant"] = m.predict_proba(design(f, f.index, ["dominant"] * len(f)))[:, 1]
f["score_recessive"] = f[["frac_states_protein_contact_r", "phylop447_r"]].mean(1)
for c in ["score_dominant", "score_recessive"]:
    f[c + "_gene_pct"] = f.groupby("gene_name")[c].rank(pct=True)
L = pd.read_csv(config.CURATION / config.LABELS, sep="\t")
lab = L.groupby("key").label.max()
f["training_label"] = f.key.map(lab).map({1: "pathogenic", 0: "control"}).fillna("")
f.loc[~f.gene_name.isin(set(d.gene_name)), "training_label"] = ""     # only disease genes were used for training
cols = ["gene_name", "chrom", "pos", "ref", "alt", "key", "hgvs", "vtype", "paralog_family", "family_ref_pos",
        "score_dominant", "score_dominant_gene_pct", "score_recessive", "score_recessive_gene_pct",
        "frac_states_protein_contact", "nb_frac_prot", "frac_states_snrna_contact", "max_protein_res", "phylop447",
        "cadd_phred", "training_label"]
out = f[cols].sort_values(["gene_name", "pos"])
out.to_csv(config.RESULTS / "snrnavep_scores.tsv.gz", sep="\t", index=False, float_format="%.4g")
print(f"{len(out)} variants in {out.gene_name.nunique()} genes")
print(out.groupby("gene_name").agg(n=("key", "size"), labelled=("training_label", lambda s: (s != "").sum())).to_string())
