"""Write wild-type / mutant RNA sequences of every enumerated variant in the disease genes, as input for the RNA
language-model scorer (scripts/04h_rnalm_score.py, run in a separate environment)."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.coords import revcomp
from snrna_vep.features_rna import genomic_edit_to_rna
from snrna_vep.reference import fetch_region

GENES = ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6-1", "RNU6-2", "RNU6-8", "RNU6-9", "RNU4ATAC", "RNU12"]
genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t").drop_duplicates("gene_name").set_index("gene_name")
feat = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t", usecols=["key", "gene_name"])
rows = []
for g in GENES:
    gene = genes.loc[g]
    gene = gene.copy(); gene["gene_name"] = g
    seq_plus = fetch_region(gene.chrom, gene.start, gene.end)
    wt = seq_plus if gene.strand == "+" else revcomp(seq_plus)
    for k in feat[feat.gene_name == g].key.unique():
        r = genomic_edit_to_rna(k, gene, seq_plus)
        if r is not None:
            rows.append(dict(key=k, gene_name=g, n_pos=r[1], vtype=r[2], wt=wt, mut=r[0]))
out = pd.DataFrame(rows)
out.to_csv(config.ROOT / "data/external/rnalm_inputs.tsv.gz", sep="\t", index=False)
print(out.groupby("gene_name").size().to_string())
