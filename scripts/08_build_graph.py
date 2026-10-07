"""Build the state-resolved spliceosome contact graph (nodes, RNA-RNA edges, per-position feature matrix)."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.align import FAMILY_REFERENCE
from snrna_vep.coords import revcomp
from snrna_vep.genes import PARALOG_GROUPS
from snrna_vep.graph import build, node_feature_matrix
from snrna_vep.reference import fetch_region

FAMILY_OF_GROUP = {"U4": "U4", "U4atac": "U4", "U6": "U6", "U6atac": "U6", "U5": "U5",
                   "U2": "U2", "U12": "U2", "U1": "U1", "U11": "U1", "U7": "U7"}
cat = pd.read_csv(config.PROCESSED / "small_rna_catalogue.tsv", sep="\t")
members = {g for gs in PARALOG_GROUPS.values() for g in gs}
cat = cat[cat.gene_name.isin(members)].drop_duplicates("gene_name")
rna = {g.gene_name: (lambda s: s if g.strand == "+" else revcomp(s))(fetch_region(g.chrom, int(g.start), int(g.end)))
       for g in cat.itertuples()}
cifs = sorted((config.RAW / "pdb").glob("*.cif*"))
nodes, edges = build(cifs, rna, FAMILY_REFERENCE, FAMILY_OF_GROUP, PARALOG_GROUPS)
nodes.to_csv(config.PROCESSED / "graph_nodes.tsv.gz", sep="\t", index=False)
edges.to_csv(config.PROCESSED / "graph_edges.tsv.gz", sep="\t", index=False)
fm, states, partners = node_feature_matrix(nodes)
fm.to_csv(config.PROCESSED / "graph_node_features.tsv.gz", sep="\t", index=False)
print(f"{len(nodes)} node-state rows, {len(edges)} RNA-RNA edges, {len(fm)} family positions x {fm.shape[1]-2} features")
print("states:", states)
print("partners:", partners)
print(edges.groupby(["family_a", "family_b"]).size().sort_values(ascending=False).head(12).to_string())
