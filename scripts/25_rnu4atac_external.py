"""RNU4ATAC external comparison (Benoit-Pilven et al. 2020, PLoS One 15:e0235655, S5 Table).

Added after the analysis plan was registered (the published structure score was pre-specified as a G2 baseline; the
cellular-assay validation comes from block 4 of the reviewer plan). Source numbering is NR_023343.1, 1 nt upstream of the
GENCODE model: source n.k = our n.(k-1) (reference bases match at 234/241 SNVs).
(A) Published RNAstructure bimolecule (U4atac:U6atac) score vs our held-out scores, recessive pathogenic vs controls
    (labels config.LABELS), RNU4ATAC held out of training.
(B) Spearman with the U12-type splicing reporter (single variants transfected into TALS patient fibroblasts; lower =
    less rescue = more damaging), n = 23 SNVs + 1 duplication (duplication excluded: no single node).
Scores: v2 (recessive context; trained without RNU4ATAC), contacts+phyloP mean rank (the final system's choice for
RNU4ATAC), contact fraction, phyloP447, CADD v1.7 (ours) and CADD from the source table, published structure score.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
src = open(HERE / "11_v2_model.py").read().split("results, pooled = [], []")[0]
G = {"__file__": str(HERE / "11_v2_model.py"), "__name__": "ext"}
exec(compile(src, "11_v2_model", "exec"), G)
config, d, rows_for_training, design, model = G["config"], G["d"], G["rows_for_training"], G["design"], G["model"]
SRC = config.DATA / "external" / "benoitpilven2020_S5.xlsx"
rng = np.random.default_rng(0)

tr = rows_for_training(d[d.unit != "RNU4ATAC"])
m = model().fit(design(d, tr.i, tr.ctx), tr.y)
a = d[d.gene_name == "RNU4ATAC"].copy()
a["v2_rec"] = m.predict_proba(design(d, a.index, ["recessive"] * len(a)))[:, 1]
a["combo"] = a[["frac_states_protein_contact_r", "phylop447_r"]].mean(1)
a["n_int"] = np.floor(a.n_pos).astype(int)

S = pd.read_excel(SRC, sheet_name="Scores for all variants")
S = S[S.Consequence.str.match(r"^n\.\d+[ACGT]>[ACGT]$", na=False)].copy()
S["our_n"] = S.Consequence.str.extract(r"n\.(\d+)")[0].astype(int) - 1
S["key"] = "chr2:" + (S.Position - 757576).astype(str) + ":" + S.Reference + ":" + S.Alternate
S = S.rename(columns={"Score modif structure bimolecule": "published_structure_score", "CADD phredScore": "cadd_source"})
a = a.merge(S[["key", "published_structure_score", "cadd_source"]], on="key", how="left")
print(f"source SNVs {len(S)}; matched to our variants {a.published_structure_score.notna().sum()}")
SC = ["v2_rec", "combo", "frac_states_protein_contact", "phylop447", "cadd_phred", "cadd_source", "published_structure_score"]

# (A) labels
lab = a[(a.y_path == 0) | ((a.y_path == 1) & a["mode"].str.contains("AR-"))]
rows = []
for s in SC:
    q = lab[lab[s].notna()]
    rows.append(dict(analysis="A: recessive P/LP vs controls", score=s, n=len(q), n_pos=int(q.y_path.sum()),
                     value=roc_auc_score(q.y_path, q[s]) if 0 < q.y_path.sum() < len(q) else np.nan))
qa = lab.dropna(subset=["published_structure_score"])
print(f"(A) variants with a published score: {len(qa)} ({int(qa.y_path.sum())} pathogenic)")
for s in SC:
    rows.append(dict(analysis="A': same variants as published score", score=s, n=len(qa), n_pos=int(qa.y_path.sum()),
                     value=roc_auc_score(qa.y_path, qa[s]) if qa[s].notna().all() else np.nan))

# (B) cellular assay
C = pd.read_excel(SRC, sheet_name="Results cell.assay", header=None).iloc[1:37, :4]
C.columns = ["cells", "v1", "v2", "rescue"]
C = C[C.v2.astype(str).str.strip().eq("none") & C.v1.astype(str).str.match(r"\s*n\.\d+[ACGT]>[ACGT]")].copy()
C["src_n"] = C.v1.str.extract(r"n\.(\d+)")[0].astype(int)
C["ref"], C["alt"] = C.v1.str.extract(r"([ACGT])>")[0], C.v1.str.extract(r">([ACGT])")[0]
Sx = S.set_index("Consequence")
C["key"] = [Sx.loc[v.strip(), "key"] if v.strip() in Sx.index else None for v in C.v1]
C = C.merge(a[["key"] + [s for s in SC if s not in ("published_structure_score", "cadd_source")]], on="key", how="left").merge(
    S[["key", "published_structure_score", "cadd_source"]], on="key", how="left")
C["rescue"] = C.rescue.astype(float)
print(f"(B) assay SNVs {len(C)}, with our scores {C.v2_rec.notna().sum()}")
for s in SC:
    q = C[[s, "rescue"]].dropna().to_numpy(float)
    rho = spearmanr(q[:, 0], -q[:, 1])[0]
    bs = [spearmanr(*(q[ii][:, [0, 1]] * [1, -1]).T)[0] for ii in (rng.integers(0, len(q), len(q)) for _ in range(2000))]
    rows.append(dict(analysis="B: Spearman with -reporter rescue (24-variant assay)", score=s, n=len(q), n_pos=np.nan,
                     value=rho, lo=np.nanquantile(bs, .025), hi=np.nanquantile(bs, .975)))
R = pd.DataFrame(rows)
R.to_csv(config.result("rnu4atac_external.tsv"), sep="\t", index=False)
C.to_csv(config.result("rnu4atac_assay_scores.tsv"), sep="\t", index=False)
print(R.round(3).to_string(index=False))
