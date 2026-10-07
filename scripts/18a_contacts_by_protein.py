"""Per-state contacts with named proteins (raw mmCIF entity descriptions), for the mechanistic analysis.
Same extraction as scripts/08_build_graph.py, but partners keep their entity name instead of the coarse group.
Output: data/processed/graph_nodes_named.tsv.gz"""
import re
import sys
from pathlib import Path

src = open(Path(__file__).resolve().parent / "08_build_graph.py").read().split("nodes, edges = build(")[0]
G = {"__file__": str(Path(__file__).resolve().parent / "08_build_graph.py"), "__name__": "named"}
exec(compile(src, "08_build_graph", "exec"), G)
import snrna_vep.graph as graph

graph._protein_key = lambda name: re.sub(r"[;|]", ",", name.upper().strip()) or "UNNAMED"
nodes, _ = graph.build(G["cifs"], G["rna"], G["FAMILY_REFERENCE"], G["FAMILY_OF_GROUP"], G["PARALOG_GROUPS"])
nodes.to_csv(G["config"].PROCESSED / "graph_nodes_named.tsv.gz", sep="\t", index=False)
p = nodes.partners.dropna().str.split(";").explode()
print(len(nodes), "node-state rows;", p.nunique(), "distinct protein entities")
print(p.value_counts().head(60).to_string())
