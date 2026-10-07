"""Per-nucleotide contact profiles across spliceosome cryo-EM states (mmCIF via gemmi).

snRNA chains are identified by sequence alignment against the family references (no hand-entered
chain IDs), then each nucleotide gets, per structure: number of contacting protein residues,
minimum protein distance, contacts with other snRNAs, and contacts with non-snRNA RNA (pre-mRNA).
Aggregating over states gives the "contact profile across the splicing cycle".
"""
import gemmi
import numpy as np
import pandas as pd

from .align import identity, map_positions

# Human spliceosome structures by assembly state. IDs are from the literature; verify each locally
# (the cloud session could not reach the PDB). Any extra mmCIF dropped in data/raw/pdb/ is used too.
SPLICEOSOME_STATES = {
    "3JCR": "tri-snRNP", "6QW6": "tri-snRNP", "6QX9": "pre-B", "5O9Z": "B",
    "6FF7": "Bact", "5Z56": "Bact-mature", "5XJC": "C*", "6QDV": "P",
    "6Y5Q": "U2-snRNP", "7DVQ": "minor-Bact", "8Y6O": "minor-preB",
}

CUTOFF = 4.5
RNA_RES = {"A", "C", "G", "U"}


def chain_sequence(chain):
    seq, residues = [], []
    for r in chain:
        if r.name in RNA_RES:
            seq.append(r.name)
            residues.append(r)
    return "".join(seq).replace("U", "T"), residues


def assign_snrna_chains(structure, references, min_identity=0.85, min_len=30):
    """{chain_name: (family_gene, {residue_index: ref n_pos})} for chains that match an snRNA.

    Identity is measured over the modelled residues, so partially built snRNAs still match.
    """
    out = {}
    for chain in structure[0]:
        seq, _ = chain_sequence(chain)
        if len(seq) < min_len:
            continue
        best = max(references.items(), key=lambda kv: identity(seq, kv[1], over="query"))
        if identity(seq, best[1], over="query") >= min_identity:
            out[chain.name] = (best[0], map_positions(seq, best[1]))
    return out


def nucleotide_contacts(structure, references, pdb_id="?", state="?"):
    model = structure[0]
    ns = gemmi.NeighborSearch(model, structure.cell, 6).populate()
    snrna = assign_snrna_chains(structure, references)
    rows = []
    for chain in model:
        if chain.name not in snrna:
            continue
        gene, mapping = snrna[chain.name]
        _, residues = chain_sequence(chain)
        for idx, res in enumerate(residues, start=1):
            ref_pos = mapping.get(idx)
            if ref_pos is None:
                continue
            prot, rna_sn, rna_other, dmin = set(), set(), set(), np.inf
            for atom in res:
                for mark in ns.find_atoms(atom.pos, "\0", radius=CUTOFF):
                    cra = mark.to_cra(model)
                    if cra.chain.name == chain.name:
                        continue
                    d = cra.atom.pos.dist(atom.pos)
                    key = (cra.chain.name, cra.residue.seqid.num)
                    if cra.residue.name in RNA_RES:
                        (rna_sn if cra.chain.name in snrna else rna_other).add(key)
                    elif gemmi.find_tabulated_residue(cra.residue.name) is not None and \
                            gemmi.find_tabulated_residue(cra.residue.name).is_amino_acid():
                        prot.add(key)
                        dmin = min(dmin, d)
            rows.append(dict(pdb_id=pdb_id, state=state, chain=chain.name, family_gene=gene,
                             ref_n_pos=ref_pos, n_protein_res=len(prot),
                             n_protein_chains=len({c for c, _ in prot}),
                             min_protein_dist=dmin if np.isfinite(dmin) else CUTOFF + 1,
                             n_snrna_contacts=len(rna_sn), n_other_rna_contacts=len(rna_other)))
    return pd.DataFrame(rows)


def aggregate_states(per_structure):
    """Collapse structure-level contacts into one profile per (family_gene, ref_n_pos)."""
    g = per_structure.groupby(["family_gene", "ref_n_pos"])
    return pd.DataFrame({
        "n_states_seen": g.state.nunique(),
        "frac_states_protein_contact": g.apply(lambda d: (d.n_protein_res > 0).groupby(d.state).any().mean()),
        "max_protein_res": g.n_protein_res.max(),
        "mean_protein_res": g.n_protein_res.mean(),
        "max_protein_chains": g.n_protein_chains.max(),
        "min_protein_dist": g.min_protein_dist.min(),
        "frac_states_snrna_contact": g.apply(lambda d: (d.n_snrna_contacts > 0).groupby(d.state).any().mean()),
        "max_other_rna_contacts": g.n_other_rna_contacts.max(),
    }).reset_index()


def load_structures(paths, references):
    frames = []
    for p in paths:
        pdb_id = p.stem.split(".")[0].upper()
        st = gemmi.read_structure(str(p))
        st.setup_entities()
        frames.append(nucleotide_contacts(st, references, pdb_id, SPLICEOSOME_STATES.get(pdb_id, "other")))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def family_contact_features(cif_paths, rna, family_reference, family_of_group, paralog_groups):
    """Contact profiles keyed on (paralog_family, family_ref_pos), ready to merge onto variants.

    `rna` maps gene name -> transcript sequence and should include every expressed paralog so that
    each structure chain is matched to its true gene before being projected onto family coordinates.
    """
    from .align import map_positions
    per = load_structures(cif_paths, rna)
    group_of = {g: grp for grp, members in paralog_groups.items() for g in members}
    maps = {}
    for gname in per.family_gene.unique():
        grp = group_of.get(gname)
        if grp:
            maps[gname] = (family_of_group[grp], map_positions(rna[gname], rna[family_reference[grp]]))
    per["paralog_family"] = per.family_gene.map(lambda x: maps.get(x, (None,))[0])
    per["family_ref_pos"] = [maps[g][1].get(p) if g in maps else None
                             for g, p in zip(per.family_gene, per.ref_n_pos)]
    per = per.dropna(subset=["paralog_family", "family_ref_pos"])
    agg = aggregate_states(per.assign(family_gene=per.paralog_family, ref_n_pos=per.family_ref_pos))
    agg = agg.rename(columns={"family_gene": "paralog_family", "ref_n_pos": "family_ref_pos"})
    return per, agg
