"""State-resolved spliceosome contact graph for the snRNA-VEP network.

For every cryo-EM structure: snRNA chains are matched by sequence, and each nucleotide (projected to its
paralog-family reference coordinate) gets per-state node features, the protein partners it touches
(named from the mmCIF entity descriptions), and RNA-RNA contact edges to other snRNA nucleotides.
"""
import re
from collections import defaultdict

import gemmi
import numpy as np
import pandas as pd

from .align import map_positions
from .contacts import CUTOFF, RNA_RES, SPLICEOSOME_STATES, assign_snrna_chains, chain_sequence

STATE_ORDER = ["tri-snRNP", "pre-B", "B", "Bact", "Bact-mature", "C*", "P", "U2-snRNP", "minor-preB", "minor-Bact"]


def _entity_names(path):
    """auth chain id -> protein/RNA entity description, from the mmCIF atom_site and entity tables."""
    doc = gemmi.cif.read(str(path))
    b = doc.sole_block()
    desc = dict(zip(b.find_loop("_entity.id"), b.find_loop("_entity.pdbx_description")))
    chain_ent = {}
    for auth, ent in zip(b.find_loop("_atom_site.auth_asym_id"), b.find_loop("_atom_site.label_entity_id")):
        chain_ent.setdefault(auth, ent)
    return {c: gemmi.cif.as_string(desc.get(e, "")).upper() for c, e in chain_ent.items()}


PARTNER_RULES = [  # ordered: first match wins
    (r"PRP8|SPLICING FACTOR 8\b|PRPF8", "PRPF8"), (r"BRR2|200 KDA|SNRNP200", "SNRNP200"),
    (r"SNU114|116 KDA|EFTUD2", "EFTUD2"), (r"PRP31|PRPF31", "PRPF31"), (r"PRP3\b|PRPF3\b", "PRPF3"),
    (r"PRP4\b|PRPF4\b", "PRPF4"), (r"PRP6|PRPF6|PROCESSING FACTOR 6", "PRPF6"),
    (r"NHP2-LIKE|SNU13|15\.5", "SNU13"), (r"TRI-SNRNP-ASSOCIATED PROTEIN 1|SART1|SNU66", "SART1"),
    (r"TRI-SNRNP-ASSOCIATED PROTEIN 2|SAD1|USP39", "USP39"), (r"DDX23|PRP28", "DDX23"),
    (r"40 KDA|40K", "SNRNP40"), (r"THIOREDOXIN-LIKE PROTEIN 4A|DIM1", "TXNL4A"),
    (r"RIBONUCLEOPROTEIN [DEFG]\d?\b|ASSOCIATED PROTEINS B AND|^SM[BDEFG]\d?$", "Sm"), (r"\bLSM", "LSm"),
    (r"SPLICING FACTOR 3B|PHD FINGER-LIKE DOMAIN-CONTAINING PROTEIN 5A", "SF3B"),
    (r"SPLICING FACTOR 3A", "SF3A"), (r"U2 SMALL NUCLEAR RIBONUCLEOPROTEIN A|B''", "U2-A'/B''"),
    (r"PROCESSING FACTOR 19|CELL DIVISION CYCLE 5|PLEIOTROPIC REGULATOR|CROOKED NECK|SYF1|SPF27|BUD31|"
     r"RBM22|AQUARIUS|SNW DOMAIN|ISOMERASE-LIKE 1|XAB2|CWC", "PRP19C/NTC"),
    (r"PROCESSING FACTOR 17|SLU7|PRP18|PRKR-INTERACTING|CACTIN|FAM32|DHX8|PRP22|DHX16|DHX15|PRP2\b|PRP43", "stepII/helicase"),
    (r"MAGO NASHI|RNA-BINDING PROTEIN 8A|EIF4A3|CASC3", "EJC"), (r"U1|SNRNP70|SNRPA\b|SNRPC", "U1-snRNP"),
]


def _protein_key(name):
    """Collapse an mmCIF entity description to a spliceosome partner label."""
    n = name.upper().strip()
    for pat, lab in PARTNER_RULES:
        if re.search(pat, n):
            return lab
    return "other"


def build(cif_paths, rna, family_reference, family_of_group, paralog_groups):
    group_of = {g: grp for grp, members in paralog_groups.items() for g in members}
    maps = {}
    node_rows, edge_rows = [], []
    for path in cif_paths:
        pdb_id = path.name.split(".")[0].upper()
        state = SPLICEOSOME_STATES.get(pdb_id, "other")
        st = gemmi.read_structure(str(path))
        st.setup_entities()
        names = _entity_names(path)
        model = st[0]
        ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
        snrna = assign_snrna_chains(st, rna)
        # residue (chain, seqnum) -> (family, family_ref_pos)
        res2fam = {}
        for chain in model:
            if chain.name not in snrna:
                continue
            gene, mapping = snrna[chain.name]
            grp = group_of.get(gene)
            if not grp:
                continue
            if gene not in maps:
                maps[gene] = map_positions(rna[gene], rna[family_reference[grp]])
            _, residues = chain_sequence(chain)
            for idx, res in enumerate(residues, start=1):
                rp = mapping.get(idx)
                fp = maps[gene].get(rp) if rp else None
                if fp:
                    res2fam[(chain.name, res.seqid.num)] = (family_of_group[grp], fp, gene)
        for chain in model:
            for res in chain:
                key = (chain.name, res.seqid.num)
                if key not in res2fam:
                    continue
                fam, fp, gene = res2fam[key]
                partners, rna_nb, other_rna, dmin = defaultdict(int), set(), 0, np.inf
                for atom in res:
                    for mark in ns.find_atoms(atom.pos, "\0", radius=CUTOFF):
                        cra = mark.to_cra(model)
                        if cra.chain.name == chain.name and abs(cra.residue.seqid.num - res.seqid.num) <= 1:
                            continue
                        k2 = (cra.chain.name, cra.residue.seqid.num)
                        if cra.residue.name in RNA_RES:
                            if k2 in res2fam:
                                rna_nb.add(k2)
                            elif cra.chain.name != chain.name:
                                other_rna += 1
                        elif gemmi.find_tabulated_residue(cra.residue.name) is not None and \
                                gemmi.find_tabulated_residue(cra.residue.name).is_amino_acid():
                            partners[_protein_key(names.get(cra.chain.name, ""))] += 1
                            dmin = min(dmin, cra.atom.pos.dist(atom.pos))
                node_rows.append(dict(pdb_id=pdb_id, state=state, family=fam, ref_pos=fp, gene=gene,
                                      n_protein_atoms=sum(partners.values()), n_partners=len(partners),
                                      partners=";".join(sorted(partners)), min_protein_dist=min(dmin, CUTOFF + 1),
                                      n_rna_nb=len(rna_nb), other_rna=int(other_rna > 0)))
                for k2 in rna_nb:
                    f2, p2, g2 = res2fam[k2]
                    edge_rows.append(dict(state=state, pdb_id=pdb_id, family_a=fam, pos_a=fp, family_b=f2, pos_b=p2,
                                          intermolecular=k2[0] != chain.name))
    nodes = pd.DataFrame(node_rows)
    edges = pd.DataFrame(edge_rows).drop_duplicates(["state", "family_a", "pos_a", "family_b", "pos_b"])
    return nodes, edges


def node_feature_matrix(nodes, top_partners=None):
    """One row per (family, ref_pos): per-state contact features + fraction of states touching each partner."""
    states = [s for s in STATE_ORDER if s in set(nodes.state)]
    g = nodes.groupby(["family", "ref_pos", "state"]).agg(
        prot=("n_protein_atoms", "mean"), npart=("n_partners", "mean"), dmin=("min_protein_dist", "min"),
        rnb=("n_rna_nb", "mean"), orna=("other_rna", "max")).reset_index()
    g["prot"] = np.log1p(g.prot)
    wide = g.pivot_table(index=["family", "ref_pos"], columns="state",
                         values=["prot", "npart", "dmin", "rnb", "orna"])
    wide.columns = [f"st_{v}_{s}" for v, s in wide.columns]
    seen = g.pivot_table(index=["family", "ref_pos"], columns="state", values="prot", aggfunc="size").notna()
    seen.columns = [f"st_seen_{s}" for s in seen.columns]
    # partner identity across states
    p = nodes.assign(partner=nodes.partners.str.split(";")).explode("partner")
    p = p[p.partner.notna() & (p.partner != "")]
    if top_partners is None:
        top_partners = [x for x in p.partner.value_counts().index if x != "other"][:20]
    pp = p[p.partner.isin(top_partners)].groupby(["family", "ref_pos", "partner"]).state.nunique().unstack(fill_value=0)
    nstate = nodes.groupby(["family", "ref_pos"]).state.nunique()
    pp = pp.div(nstate.reindex(pp.index), axis=0)
    pp.columns = [f"partner_{c}" for c in pp.columns]
    out = wide.join(seen, how="left").join(pp, how="left").reset_index()
    fill = {c: (5.5 if c.startswith("st_dmin") else 0.0) for c in out.columns if c not in ("family", "ref_pos")}
    return out.fillna(fill), states, top_partners


def node_feature_matrix_agnostic(nodes, edges):
    """Family-agnostic node features: summaries over the states in which a nucleotide is resolved, plus their
    percentile ranks within the nucleotide's own snRNA family. No state- or partner-identity columns, so features
    mean the same thing in U1, U2, U4, U5 and U6 (needed to generalise to a held-out family)."""
    n = nodes.copy()
    n["prot"] = np.log1p(n.n_protein_atoms)
    n["has_prot"] = (n.n_protein_atoms > 0).astype(float)
    n["has_rna"] = (n.n_rna_nb > 0).astype(float)
    inter = edges[edges.intermolecular].groupby(["family_a", "pos_a"]).state.nunique().rename("n_states_inter_rna")
    intra = edges[~edges.intermolecular].groupby(["family_a", "pos_a"]).state.nunique().rename("n_states_intra_rna")
    g = n.groupby(["family", "ref_pos"]).agg(
        frac_prot=("has_prot", "mean"), mean_prot=("prot", "mean"), max_prot=("prot", "max"),
        mean_partners=("n_partners", "mean"), max_partners=("n_partners", "max"),
        min_dist=("min_protein_dist", "min"), frac_rna=("has_rna", "mean"), mean_rna_nb=("n_rna_nb", "mean"),
        other_rna=("other_rna", "max"), n_states=("state", "nunique"))
    fam_states = n.groupby("family").state.nunique()
    g = g.reset_index()
    g["frac_resolved"] = g.n_states / g.family.map(fam_states)
    g = g.merge(inter.reset_index().rename(columns={"family_a": "family", "pos_a": "ref_pos"}), how="left",
                on=["family", "ref_pos"]).merge(
        intra.reset_index().rename(columns={"family_a": "family", "pos_a": "ref_pos"}), how="left", on=["family", "ref_pos"])
    g[["n_states_inter_rna", "n_states_intra_rna"]] = g[["n_states_inter_rna", "n_states_intra_rna"]].fillna(0)
    g["frac_inter_rna"] = g.n_states_inter_rna / g.n_states
    g["frac_intra_rna"] = g.n_states_intra_rna / g.n_states
    base = ["frac_prot", "mean_prot", "max_prot", "mean_partners", "max_partners", "min_dist", "frac_rna",
            "mean_rna_nb", "frac_inter_rna", "frac_intra_rna", "other_rna", "frac_resolved"]
    for c in base:
        g[f"rank_{c}"] = g.groupby("family")[c].rank(pct=True)
    return g[["family", "ref_pos"] + base + [f"rank_{c}" for c in base]].drop(columns=["n_states"], errors="ignore")
