"""Allele-specific scoring (ANALYSIS_PLAN.md section 3, pre-specified and previously unrun).

The multi-state contact model scores positions, so all three substitutions at a nucleotide receive one score.
This script adds allele-resolved structural features (snrna_vep/allele_features.py) and tests them where the
answer is measurable: WITHIN positions of RNU4-2 where all three alternative alleles were assayed by saturation
genome editing (141 such positions). A position-level score has, by construction, zero within-position signal;
any correlation here is attributable to allele identity alone.

Primary metric (pre-specified): mean within-position Spearman between score and -SGE across positions with all
three alleles measured, plus the fraction of positions where the most damaging allele is ranked first.
Comparators: CADD (allele-resolved), phyloP (position-level, so a negative control), and the position-level
contact score (also a negative control by construction).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.allele_features import aggregate_alleles, residue_allele_features

from snrna_vep.coords import revcomp
_ref = pd.read_csv(config.PROCESSED / "reference_core_loci.tsv.gz", sep="\t")
_strand = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t").set_index("gene_name").strand.to_dict()
# reference_core_loci holds plus-strand genomic sequence; the pipeline (04c) reverse-complements minus-strand
# genes to recover the transcript. Without this, minus-strand genes (RNU4-2 among them) never match a chain.
# reference_core_loci stores each gene with 2000 bp flanks on each side; trim to the transcript before aligning,
# otherwise mapped positions come out offset by exactly 2000.
FLANK = 2000
_len = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="	").set_index("gene_name").length.to_dict()
def _transcript(name, seq):
    body = seq[FLANK:FLANK + int(_len[name])] if name in _len else seq
    return body if _strand.get(name, "+") == "+" else revcomp(body)
REF = {r.gene_name: _transcript(r.gene_name, r.seq) for r in _ref.itertuples() if r.gene_name in _len}
cifs = sorted((config.RAW / "pdb").glob("*.cif*"))
print(f"structures: {len(cifs)}", flush=True)
frames = []
for p in cifs:
    pdb_id = p.name.split(".")[0].upper()
    f = residue_allele_features(p, REF, pdb_id, __import__("snrna_vep.contacts", fromlist=["x"]).SPLICEOSOME_STATES.get(pdb_id, "other"))
    frames.append(f)
    print(f"  {pdb_id}: {len(f)} nucleotides", flush=True)
PER = pd.concat(frames, ignore_index=True)
PER.to_csv(config.PROCESSED / "allele_per_structure.tsv.gz", sep="\t", index=False)
print(f"\nper-structure nucleotide records: {len(PER)}")
print(f"  paired in >=1 state: {PER.paired.mean():.2f}; with base-protein contact: {(PER.n_base_prot > 0).mean():.2f}; "
      f"backbone-only contact: {((PER.n_base_prot == 0) & (PER.n_backbone_prot > 0)).mean():.2f}")

A = aggregate_alleles(PER)
# project gene-level coordinates onto the paralog-family reference, as 04c does
from snrna_vep.align import FAMILY_REFERENCE, map_positions
from snrna_vep.genes import PARALOG_GROUPS
_group_of = {g: grp for grp, mem in PARALOG_GROUPS.items() for g in mem}
_FAM = {"U4": "U4", "U4atac": "U4", "U6": "U6", "U6atac": "U6", "U5": "U5", "U2": "U2", "U12": "U2", "U1": "U1", "U11": "U1", "U7": "U7"}
_maps = {}
for _g in A.family.unique():
    _grp = _group_of.get(_g)
    if _grp and _grp in FAMILY_REFERENCE:
        _maps[_g] = (_FAM[_grp], map_positions(REF[_g], REF[FAMILY_REFERENCE[_grp]]))
A["paralog_family"] = A.family.map(lambda g: _maps.get(g, (None,))[0])
A["family_ref_pos"] = [_maps[g][1].get(int(p)) if g in _maps else None for g, p in zip(A.family, A.ref_pos)]
A = A.dropna(subset=["paralog_family", "family_ref_pos"])
A["ref_pos"] = A.family_ref_pos.astype(int)
A = A.sort_values("n_states", ascending=False).drop_duplicates(["paralog_family", "ref_pos", "alt"])
A["family"] = A.paralog_family
A.to_csv(config.PROCESSED / "allele_features.tsv.gz", sep="\t", index=False)
print(f"per-(position, allele) feature rows: {len(A)}")

# ---- join to RNU4-2 SGE, which measures all three alleles at 141 positions
feat = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t", low_memory=False)
sge = pd.read_csv(config.DATA / "external" / "sge_rnu4-2.tsv", sep="\t")
V = feat[(feat.gene_name == "RNU4-2") & (feat.vtype == "snv")].drop(columns=["sge_score", "sge_type"], errors="ignore").merge(sge, on="key")
V["alt_rna"] = V.hgvs.str.extract(r">([ACGT])")[0].str.replace("T", "U")
V["ref_pos"] = V.family_ref_pos
V = V.merge(A[A.family == "U4"], left_on=["ref_pos", "alt_rna"], right_on=["ref_pos", "alt"], how="left")
print(f"\nRNU4-2 SNVs with SGE and allele features: {V.frac_states_pair_broken_by_alt.notna().sum()} of {len(V)}")

n = V.groupby("n_pos").size()
full = V[V.n_pos.isin(n[n == 3].index) & V.frac_states_pair_broken_by_alt.notna()].copy()
pos_ok = [p for p, g in full.groupby("n_pos") if len(g) == 3]
full = full[full.n_pos.isin(pos_ok)]
print(f"positions with all 3 alleles measured AND structural features: {len(pos_ok)}")

full["allele_score"] = (full.frac_states_pair_broken_by_alt * full.frac_states_base_contact
                        + 0.5 * full.frac_states_pair_broken_by_alt
                        + 0.25 * full.frac_states_wobble_by_alt)
SCORES = {"allele_score (new)": "allele_score",
          "pair broken by alt": "frac_states_pair_broken_by_alt",
          "CADD": "cadd_phred",
          "phyloP447 (position-level control)": "phylop447",
          "contact fraction (position-level control)": "frac_states_protein_contact"}
print("\n=== WITHIN-POSITION allele ranking against -SGE (141 positions, 3 alleles each) ===")
rows = []
for lab, col in SCORES.items():
    rhos, first = [], []
    for p, g in full.groupby("n_pos"):
        if g[col].nunique() < 2:
            continue
        r = spearmanr(g[col], -g.sge_score)[0]
        if np.isfinite(r):
            rhos.append(r)
        first.append(int(g.loc[g[col].idxmax(), "sge_score"] == g.sge_score.min()))
    rows.append(dict(score=lab, n_positions_informative=len(rhos), mean_within_rho=np.mean(rhos) if rhos else np.nan,
                     top_allele_correct=np.mean(first) if first else np.nan))
R = pd.DataFrame(rows)
R["expected_by_chance"] = 1 / 3
print(R.round(3).to_string(index=False))

best = full.groupby("n_pos").apply(lambda g: spearmanr(g.allele_score, -g.sge_score)[0]
                                   if g.allele_score.nunique() > 1 else np.nan, include_groups=False).dropna()
if len(best) > 5:
    st, p = wilcoxon(best, alternative="greater")
    print(f"\nwithin-position rho for the allele score > 0: Wilcoxon p = {p:.4f} (n = {len(best)} informative positions)")
R.to_csv(config.RESULTS / "allele_specific_within_position.tsv", sep="\t", index=False)
full.to_csv(config.RESULTS / "allele_specific_variants.tsv", sep="\t", index=False)
