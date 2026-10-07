"""Compute external baseline scores for every enumerated variant (run locally; needs internet/GPU).

Usage: python scripts/04b_external_scores.py cadd phylop [rnalm] [evo2] [alphagenome=/path/to/atlas.tsv]
Writes data/external/<name>.tsv, which scripts/04_build_features.py merges on the next run.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config, external
from snrna_vep.coords import revcomp
from snrna_vep.reference import fetch_region

feat = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t",
                   usecols=["gene_name", "key", "hgvs", "chrom", "pos", "ref", "alt", "vtype"])
out = config.DATA / "external"
out.mkdir(parents=True, exist_ok=True)

for arg in sys.argv[1:]:
    name, _, path = arg.partition("=")
    if name == "cadd":
        df = external.cadd_scores(feat[feat.vtype == "snv"])
    elif name == "phylop":
        df = external.phylop_scores(feat)
    elif name == "alphagenome":
        df = external.alphagenome_atlas(path)
    elif name == "rnalm":
        genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t")
        rna = {g.gene_name: (lambda s: s if g.strand == "+" else revcomp(s))(fetch_region(g.chrom, g.start, g.end))
               for g in genes.itertuples()}
        df = external.rna_lm_delta_ll(feat[feat.gene_name.isin(rna)], rna)
    elif name == "evo2":
        df = external.evo2_delta_ll(feat, fetch_region)
    else:
        raise SystemExit(f"unknown score {name}")
    df.to_csv(out / f"{name}.tsv", sep="\t", index=False)
    print(f"{name}: {df.iloc[:, 1].notna().sum()} scored -> {out / (name + '.tsv')}")
