"""Label-free signal check: do structure features track within-gene gnomAD depletion?

Per gene and position, correlate the mean feature over that position's SNVs with the 10-nt
relative O/E. A negative Spearman (more damaging feature -> lower O/E) means the feature carries
constraint signal without using any patient labels.
"""
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

feat = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t")
snv = feat[feat.vtype == "snv"].copy()
snv["abs_ddG_fold"] = snv.ddG_fold.abs()
snv["paired_wt"] = 1 - snv.p_unpaired_wt
FEATS = ["abs_ddG_fold", "ddG_duplex", "paired_wt"] + [c for c in ["phylop", "max_protein_res", "frac_states_protein_contact", "frac_states_snrna_contact", "max_protein_chains"] if c in snv]

per_pos = snv.groupby(["gene_name", "n_pos"]).agg({**{f: "mean" for f in FEATS}, "rel_oe_w10": "first"}).reset_index()
rows = []
for g, d in [("all", per_pos)] + list(per_pos.groupby("gene_name")):
    for f in FEATS:
        x = d[[f, "rel_oe_w10"]].dropna()
        if len(x) > 20:
            rho, p = spearmanr(x[f], x.rel_oe_w10)
            rows.append(dict(gene=g, feature=f, n_pos=len(x), spearman=rho, p=p))
out = pd.DataFrame(rows)
out.to_csv(config.RESULTS / "label_free_feature_vs_depletion.tsv", sep="\t", index=False)
print(out[out.gene.isin(["all", "RNU4-2", "RNU2-2P", "RNU5B-1", "RNU5A-1", "RNU4ATAC"])].round(4).to_string(index=False))

# Smoothed O/E is autocorrelated along the RNA, so naive Spearman p-values are optimistic.
# Circular-shift null: rotate the depletion track against the feature track (keeps both
# autocorrelation structures intact) and recompute rho for every non-trivial shift.
import numpy as np

perm_rows = []
for g, d in per_pos.groupby("gene_name"):
    d = d.sort_values("n_pos")
    for f in FEATS:
        x = d[[f, "rel_oe_w10"]].dropna()
        n = len(x)
        if n < 40:
            continue
        a, b = x[f].to_numpy(), x.rel_oe_w10.to_numpy()
        obs = spearmanr(a, b)[0]
        null = np.array([spearmanr(a, np.roll(b, k))[0] for k in range(10, n - 10)])
        perm_rows.append(dict(gene=g, feature=f, n_pos=n, spearman=obs,
                              p_shift=(1 + (null <= obs).sum()) / (1 + len(null))))
perm = pd.DataFrame(perm_rows)
perm.to_csv(config.RESULTS / "label_free_feature_vs_depletion_shiftperm.tsv", sep="\t", index=False)
print("\n== circular-shift permutation (one-sided: feature predicts depletion) ==")
print(perm[perm.gene.isin(["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU5A-1", "RNU4ATAC"])].round(4).to_string(index=False))
