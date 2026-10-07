"""Merge data/external/*.tsv (CADD, phyloP, RNA-LM, Evo 2, AlphaGenome) into the feature table
without recomputing RNA folding. Re-running replaces previously merged columns."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

path = config.PROCESSED / "variant_features.tsv.gz"
feat = pd.read_csv(path, sep="\t")
for f in sorted((config.DATA / "external").glob("*.tsv")):
    e = pd.read_csv(f, sep="\t").drop_duplicates("key")
    feat = feat.drop(columns=[c for c in e.columns if c != "key" and c in feat])
    feat = feat.merge(e, on="key", how="left")
    print(f"{f.name}: {e.columns[1:].tolist()} on {feat[e.columns[1]].notna().sum()} / {len(feat)} variants")
feat.to_csv(path, sep="\t", index=False)
