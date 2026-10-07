"""Offline tests (no network): coordinates, variant enumeration, context model, model CV, metrics."""
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep.coords import g_to_n, hgvs_n_to_vcf, n_to_g, revcomp
from snrna_vep.evaluate import metrics, paired_bootstrap_delta
from snrna_vep.model import leave_one_gene_out
from snrna_vep.mutation_model import context_key, fit, observed_expected, site_table
from snrna_vep.variants import apply_variant, enumerate_variants

# toy chromosome: positions 1..60, gene on 21..40
CHROM = "ACGTTGCAAGGCTTACGATCGATCGGATCCATTTAAACCCGGGTTTACGATGCATGCATG"


def fetch(chrom, pos):
    return CHROM[pos - 1]


def gene(strand):
    return SimpleNamespace(gene_name="G", chrom="t", start=21, end=40, strand=strand)


def test_coordinate_roundtrip():
    for strand in "+-":
        for n in range(1, 21):
            assert g_to_n(n_to_g(n, 21, 40, strand), 21, 40, strand) == n


def test_snv_minus_strand():
    g = gene("-")
    rna = revcomp(CHROM[20:40])
    c, p, ref, alt = hgvs_n_to_vcf(f"n.1{rna[0]}>A", g.chrom, g.start, g.end, g.strand, fetch)
    assert p == 40 and ref == CHROM[39] and alt == "T"


def test_insertion_left_normalized_in_homopolymer():
    # positions 33-35 are TTT on the plus strand (CATTTAAA): inserting T anywhere in it is one variant
    g = SimpleNamespace(gene_name="G", chrom="t", start=21, end=40, strand="+")
    keys = {hgvs_n_to_vcf(f"n.{a}_{a + 1}insT", "t", 21, 40, "+", fetch) for a in (12, 13, 14)}
    assert len(keys) == 1


def test_enumeration_and_apply():
    for strand in "+-":
        g = gene(strand)
        v = enumerate_variants(g, CHROM[20:40], fetch)
        assert (v.vtype == "snv").sum() == 60
        assert v.key.is_unique
    rna = "ACGU".replace("U", "T")
    assert apply_variant(rna, "n.2C>G") == "AGGT"
    assert apply_variant(rna, "n.2del") == "AGT"
    assert apply_variant(rna, "n.2_3insA") == "ACAGT"


def test_context_model():
    assert context_key("ACG", "T") == context_key(revcomp("ACG"), revcomp("T"))
    s = site_table("t", 1, 60, CHROM, observed_snvs={(10, "A"), (11, "T")})
    m = fit(s)
    per = observed_expected(s, m)
    assert abs(per.obs.sum() - 2) < 1e-9 and per.exp.sum() > 0


def test_logo_learns_signal():
    rng = np.random.default_rng(0)
    rows = []
    for g in ["A", "B", "C", "D"]:
        for i in range(200):
            y = int(rng.random() < 0.15)
            rows.append(dict(gene_name=g, key=f"{g}{i}", label=y,
                             p_unpaired_wt=rng.normal(1.5 * y, 1), ddG_fold=rng.normal(0, 1),
                             ddG_duplex=np.nan, rel_pos=rng.random(), paralog_family="U5",
                             family_ref_pos=rng.integers(1, 100)))
    p = leave_one_gene_out(pd.DataFrame(rows), ["position", "structure2d"], "logistic")
    assert metrics(p.label, p.pred)["auroc"] > 0.75


def test_paired_bootstrap():
    rng = np.random.default_rng(1)
    y = rng.integers(0, 2, 400)
    good, bad = y + rng.normal(0, 0.5, 400), rng.normal(0, 1, 400)
    d, lo, hi, p = paired_bootstrap_delta(y, good, bad, n=300)
    assert d > 0 and lo > 0 and p < 0.05
