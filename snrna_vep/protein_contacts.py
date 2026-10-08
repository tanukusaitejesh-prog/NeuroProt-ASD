"""Per-residue contact profiles for spliceosome PROTEIN subunits across cryo-EM states.

Protein-side analogue of contacts.py. Chains are identified by the mmCIF _struct_ref / _struct_ref_seq blocks
(UniProt accession plus the author-numbering offset), so residue numbering comes from the deposition itself rather
than from sequence alignment. For each protein residue that maps to a UniProt position, per structure:
  n_protein_res / n_protein_chains : contacting residues and chains of OTHER proteins
  partners                         : UniProt accessions of contacted proteins (and snRNA families)
  n_snrna_contacts                 : contacting nucleotides of snRNA chains
  n_premrna_contacts               : contacting nucleotides of non-snRNA RNA (pre-mRNA / intron)
  min_protein_dist / min_rna_dist  : closest heavy-atom distances
snRNA chains are labelled by family using the same references as contacts.py, so "which snRNA does this residue hold"
is answerable. Aggregation over states mirrors aggregate_states().
"""
import gemmi
import numpy as np
import pandas as pd

from .align import identity
from .contacts import CUTOFF, RNA_RES, SPLICEOSOME_STATES, chain_sequence

AA_MIN_LEN = 20


def uniprot_map(cif_path):
    """{chain: (uniprot_acc, offset)} where uniprot_pos = auth_seq_id + offset, from _struct_ref_seq."""
    block = gemmi.cif.read(str(cif_path)).sole_block()
    acc = {}
    for r in block.find("_struct_ref.", ["id", "db_name", "pdbx_db_accession"]):
        if r[1].strip('"\'') == "UNP":
            acc[r[0]] = r[2].strip('"\'')
    out = {}
    cols = ["ref_id", "pdbx_strand_id", "db_align_beg", "pdbx_auth_seq_align_beg"]
    for r in block.find("_struct_ref_seq.", cols):
        if r[0] not in acc:
            continue
        try:
            db_beg, auth_beg = int(r[2]), int(r[3])
        except ValueError:
            continue
        for ch in r[1].strip('"\'').split(","):
            out.setdefault(ch.strip(), (acc[r[0]], db_beg - auth_beg))
    return out


def classify_chains(structure, rna_references):
    """{chain: ('protein', acc) | ('snrna', family_gene) | ('rna_other', None)}."""
    kinds = {}
    for chain in structure[0]:
        seq, res = chain_sequence(chain)
        n_rna = sum(1 for r in chain if r.name in RNA_RES)
        n_aa = sum(1 for r in chain if gemmi.find_tabulated_residue(r.name) is not None
                   and gemmi.find_tabulated_residue(r.name).is_amino_acid())
        if n_rna >= 10 and n_rna > n_aa:
            best = max(rna_references.items(), key=lambda kv: identity(seq, kv[1], over="query")) if seq else (None, "")
            ok = best[0] and identity(seq, rna_references[best[0]], over="query") >= 0.85
            kinds[chain.name] = ("snrna", best[0]) if ok else ("rna_other", None)
        elif n_aa >= AA_MIN_LEN:
            kinds[chain.name] = ("protein", None)
    return kinds


def residue_contacts(cif_path, rna_references, pdb_id, state):
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    st.remove_ligands_and_waters()
    model = st[0]
    kinds = classify_chains(st, rna_references)
    umap = uniprot_map(cif_path)
    ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
    label = {}
    for ch, (kind, fam) in kinds.items():
        if kind == "protein":
            label[ch] = umap.get(ch, (f"PDB:{pdb_id}:{ch}", 0))[0]
        else:
            label[ch] = f"snRNA:{fam}" if kind == "snrna" else "RNA:other"
    rows = []
    for chain in model:
        if kinds.get(chain.name, (None,))[0] != "protein" or chain.name not in umap:
            continue
        acc, off = umap[chain.name]
        for res in chain:
            tab = gemmi.find_tabulated_residue(res.name)
            if tab is None or not tab.is_amino_acid():
                continue
            prot, snr, pre, partners = set(), set(), set(), set()
            dmin_p, dmin_r = np.inf, np.inf
            for atom in res:
                for mark in ns.find_atoms(atom.pos, "\0", radius=CUTOFF):
                    cra = mark.to_cra(model)
                    if cra.chain.name == chain.name or cra.chain.name not in kinds:
                        continue
                    d = cra.atom.pos.dist(atom.pos)
                    key = (cra.chain.name, cra.residue.seqid.num)
                    k = kinds[cra.chain.name][0]
                    if k == "protein":
                        prot.add(key)
                        dmin_p = min(dmin_p, d)
                        partners.add(label[cra.chain.name])
                    elif k == "snrna":
                        snr.add(key)
                        dmin_r = min(dmin_r, d)
                        partners.add(label[cra.chain.name])
                    else:
                        pre.add(key)
                        dmin_r = min(dmin_r, d)
                        partners.add("RNA:other")
            rows.append(dict(pdb_id=pdb_id, state=state, chain=chain.name, uniprot=acc,
                             uniprot_pos=res.seqid.num + off, aa=gemmi.find_tabulated_residue(res.name).one_letter_code.upper(),
                             n_protein_res=len(prot), n_protein_chains=len({c for c, _ in prot}),
                             n_snrna_contacts=len(snr), n_premrna_contacts=len(pre),
                             min_protein_dist=dmin_p if np.isfinite(dmin_p) else CUTOFF + 1,
                             min_rna_dist=dmin_r if np.isfinite(dmin_r) else CUTOFF + 1,
                             partners="|".join(sorted(partners))))
    return pd.DataFrame(rows)


def load_all(cif_paths, rna_references):
    frames = []
    for p in cif_paths:
        pdb_id = p.name.split(".")[0].upper()
        frames.append(residue_contacts(p, rna_references, pdb_id, SPLICEOSOME_STATES.get(pdb_id, "other")))
        print(f"  {pdb_id}: {len(frames[-1])} protein residues mapped to UniProt", flush=True)
    return pd.concat(frames, ignore_index=True)


def aggregate(per):
    """One profile per (uniprot, uniprot_pos) across states."""
    g = per.groupby(["uniprot", "uniprot_pos"])
    out = pd.DataFrame({
        "n_states_seen": g.state.nunique(),
        "frac_states_protein_contact": g.apply(lambda d: (d.n_protein_res > 0).groupby(d.state).any().mean(), include_groups=False),
        "frac_states_snrna_contact": g.apply(lambda d: (d.n_snrna_contacts > 0).groupby(d.state).any().mean(), include_groups=False),
        "frac_states_premrna_contact": g.apply(lambda d: (d.n_premrna_contacts > 0).groupby(d.state).any().mean(), include_groups=False),
        "max_protein_res": g.n_protein_res.max(),
        "max_snrna_contacts": g.n_snrna_contacts.max(),
        "max_premrna_contacts": g.n_premrna_contacts.max(),
        "min_protein_dist": g.min_protein_dist.min(),
        "min_rna_dist": g.min_rna_dist.min(),
    }).reset_index()
    return out
