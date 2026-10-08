"""Per-residue NMDA receptor contact profiles across gating states (data/curation/nmdar_states.tsv).

For every GluN subunit residue, per structure: contacts with the OTHER subunits, and minimum distance to the
agonists (glutamate / glycine), to channel blockers and allosteric modulators, and to any non-glycan ligand.
Chains are anchored to UniProt by the deposition's own _struct_ref blocks (snrna_vep.protein_contacts.uniprot_map).
A tetramer carries two copies of each subunit, so values are aggregated over copies within a structure (max).

Aggregation across states yields the features a single-structure method cannot have:
  frac_states_subunit_contact   fraction of states where the residue touches another subunit
  sd_subunit_contact            SD of the inter-subunit contact count across states (state-dependence)
  delta_open_nonactive          mean contacts in open/closed(active) states minus non-active/pre-active states
  min_dist_agonist / _ligand    static ligand proximity (the published baseline feature)
Output: data/processed/nmdar_residue_features.tsv.gz
"""
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.protein_contacts import uniprot_map

CUTOFF = 4.5
SUBUNITS = {"Q05586": "GluN1", "Q12879": "GluN2A", "Q13224": "GluN2B", "O15399": "GluN2D", "Q14957": "GluN2C"}
AGONISTS = {"GLU", "GLY"}
SKIP_HET = {"HOH", "NAG", "BMA", "MAN", "FUC", "GOL", "EDO", "SO4", "ACT", "CL", "NA", "ZN", "POV", "PLM", "CLR"}
ACTIVE_STATES = {"open", "closed"}
INACTIVE_STATES = {"non-active", "pre-active"}

states = pd.read_csv(config.CURATION / "nmdar_states.tsv", sep="\t")
rows = []
for r in states.itertuples():
    p = config.RAW / "nmdar" / f"{r.pdb_id}.cif.gz"
    if not p.exists():
        print(f"  {r.pdb_id}: missing", flush=True)
        continue
    st = gemmi.read_structure(str(p))
    st.setup_entities()
    st.remove_waters()
    model = st[0]
    umap = uniprot_map(p)
    sub_chains = {c: umap[c][0] for c in umap if umap[c][0] in SUBUNITS}
    if not sub_chains:
        print(f"  {r.pdb_id}: no GluN chains mapped", flush=True)
        continue
    ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
    het = {}                                   # (chain, seqid) -> comp_id for ligands
    for ch in model:
        for res in ch:
            tab = gemmi.find_tabulated_residue(res.name)
            if (tab is None or not tab.is_amino_acid()) and res.name not in SKIP_HET:
                het[(ch.name, res.seqid.num)] = res.name
    for chain in model:
        if chain.name not in sub_chains:
            continue
        acc, off = umap[chain.name]
        for res in chain:
            tab = gemmi.find_tabulated_residue(res.name)
            if tab is None or not tab.is_amino_acid():
                continue
            other, dag, dlig = set(), np.inf, np.inf
            for atom in res:
                if atom.element == gemmi.Element("H"):
                    continue
                for mark in ns.find_atoms(atom.pos, "\0", radius=6.0):
                    cra = mark.to_cra(model)
                    if cra.chain.name == chain.name:
                        continue
                    d = cra.atom.pos.dist(atom.pos)
                    key = (cra.chain.name, cra.residue.seqid.num)
                    if key in het:
                        dlig = min(dlig, d)
                        if het[key] in AGONISTS:
                            dag = min(dag, d)
                    elif cra.chain.name in sub_chains and d <= CUTOFF:
                        other.add(key)
            rows.append(dict(pdb_id=r.pdb_id, state=r.state, subunit=SUBUNITS[acc], uniprot=acc,
                             pos=res.seqid.num + off, n_subunit_contacts=len(other),
                             min_dist_agonist=min(dag, 20.0), min_dist_ligand=min(dlig, 20.0)))
    print(f"  {r.pdb_id} ({r.state}): {len(sub_chains)} GluN chains", flush=True)

per = pd.DataFrame(rows)
# collapse the two copies of each subunit within a structure
per = per.groupby(["pdb_id", "state", "subunit", "uniprot", "pos"], as_index=False).agg(
    n_subunit_contacts=("n_subunit_contacts", "max"), min_dist_agonist=("min_dist_agonist", "min"),
    min_dist_ligand=("min_dist_ligand", "min"))
per.to_csv(config.PROCESSED / "nmdar_per_structure.tsv.gz", sep="\t", index=False)

g = per.groupby(["uniprot", "subunit", "pos"])
F = pd.DataFrame({
    "n_states_seen": g.state.nunique(),
    "frac_states_subunit_contact": g.apply(lambda d: (d.n_subunit_contacts > 0).groupby(d.state).any().mean(), include_groups=False),
    "mean_subunit_contact": g.n_subunit_contacts.mean(),
    "max_subunit_contact": g.n_subunit_contacts.max(),
    "sd_subunit_contact": g.apply(lambda d: d.groupby("state").n_subunit_contacts.max().std(), include_groups=False),
    "min_dist_agonist": g.min_dist_agonist.min(),
    "min_dist_ligand": g.min_dist_ligand.min(),
}).reset_index()
act = per[per.state.isin(ACTIVE_STATES)].groupby(["uniprot", "pos"]).n_subunit_contacts.mean()
ina = per[per.state.isin(INACTIVE_STATES)].groupby(["uniprot", "pos"]).n_subunit_contacts.mean()
F["delta_open_nonactive"] = [act.get((u, p), np.nan) - ina.get((u, p), np.nan) for u, p in zip(F.uniprot, F.pos)]
F["sd_subunit_contact"] = F.sd_subunit_contact.fillna(0.0)
F.to_csv(config.PROCESSED / "nmdar_residue_features.tsv.gz", sep="\t", index=False)
print(f"\nresidues with features: {len(F)}; per subunit:\n{F.subunit.value_counts().to_string()}")
print(f"states covered per subunit:\n{per.groupby('subunit').state.nunique().to_string()}")
