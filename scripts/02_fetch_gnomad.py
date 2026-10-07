"""Fetch gnomAD v4.1 genome variants (+/- flank) and the reference for the core snRNA genes."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.gnomad import fetch_locus
from snrna_vep.reference import fetch_region

FLANK = 2000
EXTRA = ["RN7SK", "RN7SL1", "RN7SL2", "RPPH1"]

cat = pd.read_csv(config.PROCESSED / "small_rna_catalogue.tsv", sep="\t")
core = cat[cat.paralog_group.notna() | cat.is_disease_gene | cat.gene_name.isin(EXTRA)]
core = core.drop_duplicates("gene_name")
core.to_csv(config.PROCESSED / "core_genes.tsv", sep="\t", index=False)

frames, refs = [], []
for g in core.itertuples():
    df = fetch_locus(g, flank=FLANK)
    frames.append(df)
    refs.append(dict(gene_name=g.gene_name, chrom=g.chrom, start=g.start - FLANK, end=g.end + FLANK,
                     seq=fetch_region(g.chrom, g.start - FLANK, g.end + FLANK)))
    inside = df[(df.n_pos >= 1) & (df.n_pos <= g.length)]
    print(f"{g.gene_name:9s} {len(df):5d} records  ({len(inside)} in gene, "
          f"{(inside['filter'] == 'PASS').sum()} PASS, segdup={bool(inside.segdup.any())})")

pd.concat(frames).to_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t", index=False)
pd.DataFrame(refs).to_csv(config.PROCESSED / "reference_core_loci.tsv.gz", sep="\t", index=False)
print("saved")
