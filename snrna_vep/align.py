"""Paralog / analog alignment: map each gene's positions onto a family reference coordinate.

This lets evidence transfer between copies (RNU4-1 <-> RNU4-2, RNU5A-1 <-> RNU5B-1, ...) and
between major/minor spliceosome analogs (U4atac -> U4, U6atac -> U6) whose structures are homologous.
"""
from Bio import Align

FAMILY_REFERENCE = {
    "U4": "RNU4-2", "U4atac": "RNU4-2",
    "U6": "RNU6-1", "U6atac": "RNU6-1",
    "U5": "RNU5A-1",
    "U2": "RNU2-2P", "U12": "RNU2-2P",
    "U1": "RNU1-1", "U11": "RNU1-1",
    "U7": "RNU7-1",
}


def _aligner():
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.match_score, a.mismatch_score = 2, -1
    a.open_gap_score, a.extend_gap_score = -3, -1
    a.end_gap_score = 0
    return a


def map_positions(query, reference):
    """Return {query n_pos (1-based): reference n_pos or None} from a global alignment."""
    aln = _aligner().align(query, reference)[0]
    mapping = {i + 1: None for i in range(len(query))}
    for (qs, qe), (rs, re_) in zip(*aln.aligned):
        for k in range(qe - qs):
            mapping[qs + k + 1] = rs + k + 1
    return mapping


def identity(query, reference, over="longer"):
    """Identical aligned positions / length. over="query" scores partial models (e.g. a 42-nt
    fragment of U6 in a cryo-EM map) by how well the modelled part matches."""
    aln = _aligner().align(query, reference)[0]
    matches = sum(query[qs + k] == reference[rs + k]
                  for (qs, qe), (rs, _) in zip(*aln.aligned) for k in range(qe - qs))
    denom = len(query) if over == "query" else max(len(query), len(reference))
    return matches / denom
