"""Could the recurrent ReNU allele be produced by interlocus gene conversion from the tandem paralogue?

Motivation. n.64_65insT accounts for 70-77% of ReNU syndrome, and scripts/48-50 showed its sequence context is
mutationally DISfavoured (single-base insertions at a T4 tract occur at 0.40x the substitution frequency), leaving
a ~20-fold residual. Non-allelic gene conversion is the standard explanation for a recurrent allele whose
mutation rate is too low to account for it: a near-identical donor elsewhere in the genome repeatedly overwrites
the acceptor, and the "mutation" is actually a copying event. It is the mechanism behind recurrent pathogenic
alleles in CYP21A2 (from CYP21A1P), PMS2 (from PMS2CL), SMN1 (from SMN2) and GBA (from GBAP1).

RNU4-2 has exactly the configuration that makes this plausible: RNU4-1 sits 1.2 kb away on the same strand of
chr12, is the same length, and is a direct paralogue.

The test is simple and decisive. Gene conversion can only transfer sequence the donor actually has. If RNU4-1
carries five thymines where RNU4-2 carries four, conversion produces n.64_65insT for free and the recurrence is
explained. If RNU4-1 carries the same T4 tract, conversion cannot produce the insertion and the hypothesis is
dead.

The same logic is applied to every observed patient variant in the critical region: a variant is
conversion-explicable only if the donor carries the alternate allele at the aligned position.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

COMP = str.maketrans("ACGTN", "TGCAN")
FLANK = 2000                      # reference_core_loci stores each gene with +/-2 kb of flanking sequence
CR = (62, 79)                     # the 18-nt critical region in n. coordinates (chr12:120,291,825-120,291,842)

R = pd.read_csv(config.PROCESSED / "reference_core_loci.tsv.gz", sep="\t").set_index("gene_name")
G = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t").set_index("gene_name")


def mature(gene):
    """Mature (transcribed-strand) sequence: strip the stored flanks, reverse-complement if minus strand."""
    s = R.loc[gene].seq.upper()
    core = s[FLANK:len(s) - FLANK]
    return core.translate(COMP)[::-1] if G.loc[gene].strand == "-" else core


acc, don = mature("RNU4-2"), mature("RNU4-1")
print(f"RNU4-2  chr12:{G.loc['RNU4-2'].start:,}-{G.loc['RNU4-2'].end:,} ({G.loc['RNU4-2'].strand}), {len(acc)} nt")
print(f"RNU4-1  chr12:{G.loc['RNU4-1'].start:,}-{G.loc['RNU4-1'].end:,} ({G.loc['RNU4-1'].strand}), {len(don)} nt")
gap = G.loc["RNU4-1"].start - G.loc["RNU4-2"].end
print(f"separation {gap:,} bp, same strand -> a plausible conversion donor on proximity and homology\n")

n = min(len(acc), len(don))
ident = sum(a == b for a, b in zip(acc, don))
diffs = [(i + 1, a, b) for i, (a, b) in enumerate(zip(acc, don)) if a != b]
print(f"sequence identity {ident}/{n} ({ident / n:.1%})")
print(f"all differences (n.pos, RNU4-2, RNU4-1): {diffs}")
print(f"differences inside the critical region n.{CR[0]}-{CR[1]}: "
      f"{[d for d in diffs if CR[0] <= d[0] <= CR[1]] or 'NONE'}\n")

# ---- the focal question: the homopolymer that n.64_65insT extends
# the insertion is a T added to a T4 run on the GENOMIC (plus) strand; on the transcribed minus strand that run
# reads as A4 at n.65-68
print("the homopolymer that the recurrent allele extends")
print(f"  RNU4-2 n.60-75  {acc[59:75]}")
print(f"  RNU4-1 n.60-75  {don[59:75]}")


def run_at(seq, pos):
    """Length of the homopolymer run containing 1-based position pos."""
    i = pos - 1
    b, lo, hi = seq[i], i, i
    while lo > 0 and seq[lo - 1] == b:
        lo -= 1
    while hi < len(seq) - 1 and seq[hi + 1] == b:
        hi += 1
    return b, hi - lo + 1


ba, la = run_at(acc, 66)
bd, ld = run_at(don, 66)
print(f"  RNU4-2 run at n.66: {ba}x{la}   (= T{la} on the genomic plus strand)")
print(f"  RNU4-1 run at n.66: {bd}x{ld}   (= T{ld} on the genomic plus strand)")
verdict_ins = ld > la
print(f"\n  donor carries a LONGER run than the acceptor? {verdict_ins}")
print("  -> conversion from RNU4-1 CAN produce n.64_65insT" if verdict_ins else
      "  -> conversion from RNU4-1 CANNOT produce n.64_65insT: the donor has the identical tract")

# ---- generalise: is ANY observed patient variant conversion-explicable?
C = pd.read_csv(config.CURATION / "renu_variant_counts.tsv", sep="\t", comment="#")
C["ndd"] = C.gel_ndd + C.nongel_ndd
end = int(G.loc["RNU4-2"].end)
rows = []
for r in C[C.ndd > 0].itertuples():
    pos = int(r.key.split(":")[1])
    npos = end - pos + 1                                   # minus-strand transcript coordinate
    if r.vtype == "snv":
        a = acc[npos - 1] if 1 <= npos <= len(acc) else "?"
        d = don[npos - 1] if 1 <= npos <= len(don) else "?"
        alt = r.key.split(":")[3].translate(COMP)
        expl = (d == alt)
    else:
        a = d = "-"
        expl = verdict_ins
    rows.append(dict(hgvs=r.hgvs, vtype=r.vtype, n_pos=npos, ndd_carriers=int(r.ndd),
                     rnu4_2=a, rnu4_1=d, conversion_explicable=expl))
D = pd.DataFrame(rows).sort_values("ndd_carriers", ascending=False)
print("\n=== can interlocus conversion from RNU4-1 explain each observed patient variant? ===")
print(D.to_string(index=False))
print(f"\nexplicable: {int(D.conversion_explicable.sum())} of {len(D)} alleles, "
      f"{int(D[D.conversion_explicable].ndd_carriers.sum())} of {int(D.ndd_carriers.sum())} carriers")
D.to_csv(config.RESULTS / "gene_conversion_test.tsv", sep="\t", index=False)

print("\n=== VERDICT ===")
print("  Interlocus gene conversion from the tandem paralogue is EXCLUDED as the source of the recurrent allele."
      if not verdict_ins else "  Conversion remains viable - investigate further.")
