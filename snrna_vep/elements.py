"""Annotated secondary-structure elements from Rfam seed alignments (data/external/rfam/RFxxxxx.sto, downloaded
2026-10-07 from rfam.org/family/<acc>/alignment/stockholm).

Each human gene is aligned to its closest seed sequence; the consensus structure (WUSS notation) is projected onto the gene.
Element type per position from WUSS: paired stem (<>()[]{}), hairpin loop (_), bulge / interior loop (-), multiloop
junction (,), exterior single strand (:), pseudoknot (letters), unaligned. A canonical Sm-site motif (A U4-6 G R) in the
transcript is flagged separately because it is a protein-binding element not visible in the RNA-only structure.
"""
import re
from pathlib import Path

import pandas as pd
from Bio import Align

from . import config

RFAM = {"RNU4-2": "RF00015", "RNU4-1": "RF00015", "RNU2-2P": "RF00004", "RNU2-1": "RF00004",
        "RNU5B-1": "RF00020", "RNU5A-1": "RF00020", "RNU5D-1": "RF00020", "RNU5E-1": "RF00020", "RNU5F-1": "RF00020",
        "RNU6-1": "RF00026", "RNU6-2": "RF00026", "RNU6-7": "RF00026", "RNU6-8": "RF00026", "RNU6-9": "RF00026",
        "RNU4ATAC": "RF00618", "RNU12": "RF00007", "RNU6ATAC": "RF00619"}


def read_stockholm(path):
    seqs, ss = {}, ""
    for line in open(path):
        if line.startswith("#=GC SS_cons"):
            ss += line.split()[-1]
        elif line.strip() and not line.startswith(("#", "//")):
            name, s = line.split()
            seqs[name] = seqs.get(name, "") + s
    return seqs, ss


def wuss_type(c):
    if c in "<>()[]{}":
        return "stem"
    return {"_": "hairpin_loop", "-": "interior_bulge", ",": "multiloop", ":": "exterior", "~": "exterior"}.get(
        c, "pseudoknot" if c.isalpha() else "unaligned")


def transcript_seq(feat, gene):
    """Transcript-strand sequence from the enumerated SNV HGVS strings (n.<pos><ref>><alt>)."""
    x = feat[(feat.gene_name == gene) & (feat.vtype == "snv")]
    ref = {}
    for h in x.hgvs:
        m = re.match(r"n\.(\d+)([ACGT])>", h)
        if m:
            ref[int(m.group(1))] = m.group(2)
    return "".join(ref[i] for i in range(1, max(ref) + 1)).replace("T", "U")


def annotate(feat, genes):
    a = Align.PairwiseAligner()
    a.mode, a.match_score, a.mismatch_score, a.open_gap_score, a.extend_gap_score = "local", 2, -1, -3, -1
    rows = []
    for gene in genes:
        if gene not in RFAM:
            continue
        seq = transcript_seq(feat, gene)
        seqs, ss = read_stockholm(config.DATA / "external" / "rfam" / f"{RFAM[gene]}.sto")
        best = max(seqs.items(), key=lambda kv: a.score(seq, kv[1].replace("-", "").replace(".", "").upper().replace("T", "U")))
        gapped = best[1]
        cols = [i for i, c in enumerate(gapped) if c not in "-."]                 # ungapped seed index -> column
        aln = a.align(seq, gapped.replace("-", "").replace(".", "").upper().replace("T", "U"))[0]
        col_of = {}
        for (qs, qe), (rs, re_) in zip(*aln.aligned):
            for k in range(qe - qs):
                col_of[qs + k + 1] = cols[rs + k]
        sm = set()
        for m in re.finditer(r"A[AU]U{4,6}G[AG]", seq):
            sm.update(range(m.start() + 1, m.end() + 1))
        ident = sum(seq[q - 1] == gapped[c].upper().replace("T", "U") for q, c in col_of.items()) / len(seq)
        for p in range(1, len(seq) + 1):
            t = wuss_type(ss[col_of[p]]) if p in col_of else "unaligned"
            rows.append(dict(gene_name=gene, n_int=p, element=("Sm_site" if p in sm else t), wuss=ss[col_of[p]] if p in col_of else "",
                             rfam=RFAM[gene], seed=best[0], seed_identity=round(ident, 3)))
    return pd.DataFrame(rows)
