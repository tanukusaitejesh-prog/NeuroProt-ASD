"""GATE T2. Does it matter WHICH protein a nucleotide touches, beyond HOW MANY?

Motivation. scripts/18b produced large per-protein Mantel-Haenszel odds ratios for pathogenic contact (PPIL2 31.8,
SF3B1 22.2, SRRM2 17.5, PRPF8 7.7, SART1 6.7, SNRNP200 4.7), stratified by gene. That is suggestive, but it is
confounded: pathogenic variants cluster spatially within a gene, so ANY protein whose footprint covers the cluster
gets a large OR without carrying independent information. Stratifying by gene does not remove within-gene clustering.

SE(3)-PROTACs (Briefings in Bioinformatics 2026) makes partner identity central by embedding each partner protein
with ESM-2. Before building anything like that here, the premise has to be tested: does partner identity transfer?

Design. The confound is removed by requiring the information to generalise to a gene the model has never seen.
  features  base      = [burial, mean_degree]                  (count only; no identity)
            identity  = per-position occupancy of each protein partner, restricted to proteins that appear at
                        labelled positions in >=3 of the 5 labelled genes, so identity CAN in principle transfer
  test      leave-one-gene-out AUROC. If identity encodes real chemistry it should help on a held-out gene; if the
            per-protein ORs were positional clustering it will not, because the test gene's cluster was never seen.
  second    on RNU4-2 SGE function, partial Spearman of each partner's occupancy given mean_degree.

Kill criterion, declared before running: if identity does not improve held-out-gene AUROC over base, partner
identity carries no transferable information and no ESM-2 partner-embedding layer is worth building.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

_g = {"__file__": str(Path(__file__).resolve().parent / "51_trajectory_gate.py"), "__name__": "gate"}
_src = open(Path(__file__).resolve().parent / "51_trajectory_gate.py").read().split("# =============================== T1a")[0]
exec(compile(_src, "51_trajectory_gate", "exec"), _g)
nodes, V2, L = _g["nodes"], _g["V2"], _g["L"]
logo_auroc, partial_spearman = _g["logo_auroc"], _g["partial_spearman"]

# ---- partner occupancy per (family, ref_pos): fraction of that position's resolved states in which P is a partner
rows = []
for (fam, pos), g in nodes.groupby(["family", "ref_pos"]):
    n, c = len(g), {}
    for s in g.partners.dropna():
        for x in str(s).split(";"):
            if x:
                c[x] = c.get(x, 0) + 1
    rows.append(dict(paralog_family=fam, family_ref_pos=pos, **{("P_" + k): v / n for k, v in c.items()}))
P = pd.DataFrame(rows).fillna(0.0)
pcols = [c for c in P.columns if c.startswith("P_")]
print(f"distinct protein partners across all structures: {len(pcols)}")

V3 = V2.merge(P, on=["paralog_family", "family_ref_pos"], how="left")
pos = L[(L.label == 1) & L["mode"].astype(str).str.contains("NDD")]
neg = L[(L.label == 0) & L.gene_name.isin(pos.gene_name.unique())]
LL = pd.concat([pos, neg])[["key", "label"]].drop_duplicates("key")
D = V3.merge(LL, on="key", how="inner").dropna(subset=["burial", "mean_degree"])
D[pcols] = D[pcols].fillna(0.0)

seen = [c for c in pcols if (D.groupby("gene_name")[c].max() > 0).sum() >= 3]
print(f"transferable vocabulary (present in >=3 of {D.gene_name.nunique()} labelled genes): {len(seen)}")
for c in sorted(seen):
    print("   ", c[2:])

res = []
for name, feats in [("burial+degree", ["burial", "mean_degree"]),
                    ("identity ONLY", seen),
                    ("burial+degree+identity", ["burial", "mean_degree"] + seen)]:
    a, per = logo_auroc(D, feats)
    res.append(dict(model=name, logo_auroc=a, **{u: v for u, v in per}))
    print(f"\n  {name:<24s} LOGO AUROC {a:.3f}")
    for u, v in per:
        print(f"      {u:<10s} {v:.3f}")
R = pd.DataFrame(res)
R.to_csv(config.RESULTS / "partner_identity_gate.tsv", sep="\t", index=False)

print("\n=== RNU4-2 SGE function: partner occupancy given partner count ===")
S = V3[V3.sge_score.notna() & V3.burial.notna()].copy()
S["sge_score"] = pd.to_numeric(S.sge_score, errors="coerce")
S = S.dropna(subset=["sge_score"])
S[pcols] = S[pcols].fillna(0.0)
y = -S.sge_score.to_numpy()
out = []
for c in [c for c in seen if S[c].std() > 0]:
    r = spearmanr(S[c], y)
    out.append(dict(partner=c[2:], rho=r.statistic, p=r.pvalue,
                    partial_given_degree=partial_spearman(S[c].to_numpy(float), y, S.mean_degree.to_numpy(float))))
O = pd.DataFrame(out).sort_values("rho", key=lambda s: -s.abs())
print(O.round(3).to_string(index=False))
O.to_csv(config.RESULTS / "partner_identity_sge.tsv", sep="\t", index=False)

base = R.loc[R.model == "burial+degree", "logo_auroc"].iloc[0]
full = R.loc[R.model == "burial+degree+identity", "logo_auroc"].iloc[0]
print(f"\n=== VERDICT ===  base {base:.3f} -> +identity {full:.3f}  ({full - base:+.3f})")
print("  IDENTITY TRANSFERS" if full > base + 0.01 else
      "  IDENTITY DOES NOT TRANSFER - the per-protein ORs are positional clustering; no ESM-2 partner layer")
