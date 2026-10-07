"""Interim position-level conservation: phyloP (Zoonomia 241 mammals) as annotated on gnomAD records.

gnomAD only annotates positions where it observed a variant, and the unobserved positions are the most
constrained ones. A missingness flag would therefore leak gnomAD observation into the features, so gaps
are filled by linear interpolation along the transcript and no flag is kept. The full phyloP track from
scripts/04b_external_scores.py (column `phylop`) supersedes this column (`phylop_gnomadpos`).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.coords import g_to_n

gn = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")
genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t")
path = config.PROCESSED / "variant_features.tsv.gz"
feat = pd.read_csv(path, sep="\t").drop(columns=["phylop_gnomadpos"], errors="ignore")

tracks = []
for g in genes.itertuples():
    x = gn[(gn.gene_name == g.gene_name) & (gn.pos >= g.start) & (gn.pos <= g.end) & gn.phylop.notna()]
    if x.empty:
        continue
    by_pos = x.groupby("pos").phylop.first()
    n = pd.Series({g_to_n(p, g.start, g.end, g.strand): v for p, v in by_pos.items()}).sort_index()
    full = n.reindex(range(1, g.length + 1)).interpolate(limit_direction="both")
    tracks.append(pd.DataFrame({"gene_name": g.gene_name, "n_int": full.index, "phylop_gnomadpos": full.values}))
track = pd.concat(tracks)
feat["n_int"] = np.floor(feat.n_pos).astype(int).clip(lower=1)
feat = feat.merge(track, on=["gene_name", "n_int"], how="left").drop(columns="n_int")
feat.to_csv(path, sep="\t", index=False)
print(f"phylop_gnomadpos on {feat.phylop_gnomadpos.notna().sum()} / {len(feat)} variants")
