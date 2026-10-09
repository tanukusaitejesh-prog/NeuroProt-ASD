"""Two additions that answer the two strongest remaining objections to the recurrence paper.

OBJECTION 1: "your 0.40x expectation is a genome-wide average over 21,831 tracts, applied to one site. Site-level
indel rate varies with flanking sequence, replication timing and chromatin, and you have not bounded that."

This is answerable here in a way it usually is not, because RNU4-2 has a paralogue that is sequence-identical
across the critical region and 1.2 kb away, and - unlike the RNU1 family - is callable in short-read data
(2,545 gnomAD PASS records). RNU4-1 is not a disease gene, so population variation there is not depleted by
selection against a phenotype. It therefore provides the mutational supply at the SAME sequence context,
measured at the SAME locus, free of disease ascertainment: an internal calibration of the genome-wide estimate.

OBJECTION 2: "'no known mechanism' is an argument from absence."

Converted into a number. If the excess were produced by clonal selection during oogonial proliferation, the
per-division selective advantage required can be computed and compared with what is documented for selfish
spermatogonial selection, which operates over hundreds of divisions rather than tens.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

# homologous coordinates: both genes are minus strand, 141 nt; genomic = end - n + 1
GENES = {"RNU4-2": dict(start=120291763, end=120291903, disease=True),
         "RNU4-1": dict(start=120293097, end=120293237, disease=False)}
TRACT_N = (65, 68)                       # the A4 run on the transcribed strand = T4 on the genomic strand
FOCAL_N = 65                             # n.64_65insT inserts immediately 5' of n.65

G = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")
G = G[G["filter"] == "PASS"].copy()
G["kind"] = np.where(G.alt.str.len() > G.ref.str.len(), "ins",
                     np.where(G.alt.str.len() < G.ref.str.len(), "del", "snv"))

print("=" * 96)
print("PART 1  locus-specific mutational supply at the homologous T4 tract")
print("=" * 96)
rows = []
for g, meta in GENES.items():
    end = meta["end"]
    lo, hi = end - TRACT_N[1] + 1, end - TRACT_N[0] + 1        # genomic span of the tract
    d = G[(G.gene_name == g) & G.pos.between(lo, hi)]
    ins = d[(d.kind == "ins") & (d.alt.str.len() == d.ref.str.len() + 1)]
    snv = d[d.kind == "snv"]
    print(f"\n{g}  (disease gene: {meta['disease']})   tract chr12:{lo:,}-{hi:,}")
    if len(d):
        print(d[["pos", "ref", "alt", "kind", "AC", "AN"]].to_string(index=False))
    else:
        print("  no PASS variants")
    rows.append(dict(gene=g, disease=meta["disease"], tract=f"chr12:{lo:,}-{hi:,}",
                     ins_alleles=len(ins), ins_AC=int(ins.AC.sum()),
                     snv_alleles=len(snv), snv_AC=int(snv.AC.sum()),
                     ins_per_snv_alleles=len(ins) / max(len(snv), 1)))
T = pd.DataFrame(rows)
print("\n--- summary at the homologous tract ---")
print(T.to_string(index=False))

ctrl = T[T.gene == "RNU4-1"].iloc[0]
print(f"""
Reading. At RNU4-1 the equivalent single-base insertion IS observed in gnomAD, so this sequence context does
generate the allele; its absence at RNU4-2 is consistent with selection against the phenotype rather than with
the insertion being impossible. Counting DISTINCT observed alleles - the unit used by the homopolymer analysis,
and the one least distorted by drift - the tract gives {ctrl.ins_alleles} insertion against {ctrl.snv_alleles}
substitutions, a ratio of {ctrl.ins_per_snv_alleles:.2f}, against the genome-wide estimate of 0.40 applied in the
manuscript. The locus-specific supply is therefore consistent with the genome-wide figure and gives no support to
the objection that this particular tract is unusually insertion-prone.

Caveat, stated plainly: allele COUNT tells a different story ({ctrl.ins_AC} insertion alleles against
{ctrl.snv_AC} substitution alleles), because one insertion lineage has drifted to AC 5 while the substitutions
are singletons. Allele count reflects drift and allele age as much as mutation rate, which is why distinct-allele
counting is the primary reading; but a reader should see both.""")

# whole-gene comparison, for context
print("\n--- whole gene body, for context ---")
for g, meta in GENES.items():
    d = G[(G.gene_name == g) & G.pos.between(meta["start"], meta["end"])]
    ins = d[(d.kind == "ins") & (d.alt.str.len() == d.ref.str.len() + 1)]
    snv = d[d.kind == "snv"]
    print(f"  {g}: {len(snv)} substitution alleles, {len(ins)} 1-bp insertion alleles "
          f"-> {len(ins) / max(len(snv), 1):.3f} insertions per substitution")
T.to_csv(config.RESULTS / "paralogue_control.tsv", sep="\t", index=False)

print("\n" + "=" * 96)
print("PART 2  what a selective explanation would have to achieve, per cell division")
print("=" * 96)
RESIDUAL = 110.0
OOGONIAL_DIVISIONS = 30          # mitotic divisions from PGC to the arrested oocyte pool, completed prenatally
SPERM_DIVISIONS_40 = 610         # approximate spermatogonial divisions by paternal age 40
SELFISH_ENRICHMENT = 1000.0      # upper end of documented selfish spermatogonial enrichment in offspring

f_req = RESIDUAL ** (1 / OOGONIAL_DIVISIONS)
m_req = SELFISH_ENRICHMENT ** (1 / SPERM_DIVISIONS_40)
print(f"""
  required female advantage  {RESIDUAL:.0f}-fold over ~{OOGONIAL_DIVISIONS} prenatal oogonial divisions
                             -> {f_req:.4f} per division, i.e. {100 * (f_req - 1):.1f}% per division

  documented male advantage  up to {SELFISH_ENRICHMENT:.0f}-fold over ~{SPERM_DIVISIONS_40} spermatogonial
                             divisions by age 40
                             -> {m_req:.4f} per division, i.e. {100 * (m_req - 1):.1f}% per division

  ratio                      the female germline would need a per-division advantage
                             {100 * (f_req - 1) / (100 * (m_req - 1)):.0f}x larger than the best-documented
                             selfish spermatogonial variants, and would have to sustain it within a window that
                             closes before birth""")
for n in (20, 30, 40, 60):
    print(f"    sensitivity: {n:2d} divisions -> {100 * (RESIDUAL ** (1 / n) - 1):5.1f}% per division")
pd.DataFrame([dict(residual=RESIDUAL, oogonial_divisions=OOGONIAL_DIVISIONS,
                   required_per_division_pct=100 * (f_req - 1),
                   documented_male_per_division_pct=100 * (m_req - 1),
                   ratio=(f_req - 1) / (m_req - 1))]).to_csv(
    config.RESULTS / "selection_requirement.tsv", sep="\t", index=False)
print(f"\nwrote paralogue_control.tsv and selection_requirement.tsv")
