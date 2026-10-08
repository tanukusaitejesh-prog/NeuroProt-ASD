"""Allele-specific structural features for snRNA variants (pre-specified in ANALYSIS_PLAN.md section 3).

The multi-state contact model scores POSITIONS: all three substitutions at a nucleotide receive the same score.
These features distinguish alleles, which is what makes the method a variant-level rather than position-level
predictor. Everything is read from the same cryo-EM structures used for the contact profiles.

Per structure, for each snRNA nucleotide:
  base pairing   Watson-Crick edge geometry. A purine N1 within PAIR_CUT of a pyrimidine N3 (or two purines'
                 N1...N1, for non-canonical pairs) is called a pair. The partner's identity is recorded, so the
                 effect of each alternative allele can be evaluated: does it keep a canonical pair (A-U, G-C),
                 a wobble (G-U), or break the pair?
  contact target Protein contacts are split into BASE contacts (to ring/exocyclic atoms, which depend on the
                 identity of the base) and BACKBONE contacts (to phosphate and ribose, which do not). A variant
                 can only perturb a protein interaction through a base contact, so this separates allele-sensitive
                 from allele-insensitive contacts - a distinction invisible to a position-level score.

Aggregated over states, these yield per-(position, alt allele) features:
  frac_states_paired, frac_states_pair_broken_by_alt, frac_states_wobble_by_alt,
  frac_states_base_contact, max_base_contacts, frac_states_backbone_contact
"""
import gemmi
import numpy as np
import pandas as pd

from .align import identity, map_positions
from .contacts import CUTOFF, RNA_RES, SPLICEOSOME_STATES, chain_sequence

PAIR_CUT = 3.5                      # N1...N3 heavy-atom distance for a Watson-Crick-like pair
PURINE, PYRIMIDINE = {"A", "G"}, {"C", "U"}
BASE_ATOMS = {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9", "N6", "O6", "N2", "O2", "N4", "O4"}
CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
WOBBLE = {("G", "U"), ("U", "G")}


def pair_class(a, b):
    if (a, b) in CANONICAL:
        return "canonical"
    if (a, b) in WOBBLE:
        return "wobble"
    return "broken"


def _wc_atom(res):
    """The Watson-Crick edge nitrogen: N1 for purines, N3 for pyrimidines."""
    want = "N1" if res.name in PURINE else "N3"
    for a in res:
        if a.name == want:
            return a
    return None


def residue_allele_features(cif_path, references, pdb_id, state):
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    st.remove_waters()
    model = st[0]
    # identify snRNA chains the same way contacts.py does
    snrna = {}
    for chain in model:
        seq, _ = chain_sequence(chain)
        if len(seq) < 30:
            continue
        best = max(references.items(), key=lambda kv: identity(seq, kv[1], over="query"))
        if identity(seq, references[best[0]], over="query") >= 0.85:
            snrna[chain.name] = (best[0], map_positions(seq, references[best[0]]))
    if not snrna:
        return pd.DataFrame()
    ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
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
            wt = res.name
            # ---- base pairing: nearest Watson-Crick partner
            partner, pdist = None, np.inf
            wc = _wc_atom(res)
            if wc is not None:
                for mark in ns.find_atoms(wc.pos, "\0", radius=PAIR_CUT):
                    cra = mark.to_cra(model)
                    if cra.residue.name not in RNA_RES:
                        continue
                    if cra.chain.name == chain.name and abs(cra.residue.seqid.num - res.seqid.num) <= 1:
                        continue
                    want = "N1" if cra.residue.name in PURINE else "N3"
                    if cra.atom.name != want:
                        continue
                    d = cra.atom.pos.dist(wc.pos)
                    # a genuine pair needs one purine and one pyrimidine on the WC edge
                    if d < pdist and ((wt in PURINE) != (cra.residue.name in PURINE)):
                        partner, pdist = cra.residue.name, d
            # ---- protein contacts split into base vs backbone
            base_prot, bb_prot = set(), set()
            for atom in res:
                is_base = atom.name in BASE_ATOMS
                for mark in ns.find_atoms(atom.pos, "\0", radius=CUTOFF):
                    cra = mark.to_cra(model)
                    if cra.chain.name == chain.name or cra.residue.name in RNA_RES:
                        continue
                    tab = gemmi.find_tabulated_residue(cra.residue.name)
                    if tab is None or not tab.is_amino_acid():
                        continue
                    (base_prot if is_base else bb_prot).add((cra.chain.name, cra.residue.seqid.num))
            rows.append(dict(pdb_id=pdb_id, state=state, family=gene, ref_pos=ref_pos, wt=wt,
                             paired=partner is not None, partner=partner or "",
                             n_base_prot=len(base_prot), n_backbone_prot=len(bb_prot)))
    return pd.DataFrame(rows)


def aggregate_alleles(per, families=("U1", "U2", "U4", "U5", "U6")):
    """Expand per-(structure, position) records into per-(position, alt allele) features across states."""
    out = []
    for (fam, pos), g in per.groupby(["family", "ref_pos"]):
        n_states = g.state.nunique()
        wt = g.wt.mode().iat[0] if len(g.wt.mode()) else None
        if wt is None:
            continue
        by_state = g.groupby("state")
        frac_paired = by_state.paired.any().mean()
        frac_base = by_state.apply(lambda d: (d.n_base_prot > 0).any(), include_groups=False).mean()
        frac_bb = by_state.apply(lambda d: (d.n_backbone_prot > 0).any(), include_groups=False).mean()
        for alt in "ACGU":
            if alt == wt:
                continue
            broken = wob = 0
            for state, d in by_state:
                p = d[d.paired & (d.partner != "")]
                if not len(p):
                    continue
                part = p.partner.iat[0]
                was = pair_class(wt, part)
                now = pair_class(alt, part)
                if was in ("canonical", "wobble") and now == "broken":
                    broken += 1
                elif now == "wobble" and was == "canonical":
                    wob += 1
            out.append(dict(family=fam, ref_pos=pos, wt=wt, alt=alt, n_states=n_states,
                            frac_states_paired=frac_paired,
                            frac_states_pair_broken_by_alt=broken / max(n_states, 1),
                            frac_states_wobble_by_alt=wob / max(n_states, 1),
                            frac_states_base_contact=frac_base,
                            frac_states_backbone_contact=frac_bb,
                            max_base_contacts=g.n_base_prot.max(),
                            max_backbone_contacts=g.n_backbone_prot.max(),
                            transition=int((wt in PURINE) == (alt in PURINE))))
    return pd.DataFrame(out)
