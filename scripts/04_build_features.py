"""Enumerate every SNV and 1-nt indel in the callable core snRNA genes and attach features.

Always computed: position, RNA structure (ViennaRNA), paralog-family alignment, within-gene depletion.
Merged when present locally: structure contacts (data/raw/pdb/*.cif[.gz]) and external scores
(data/external/<name>.tsv with a `key` column; produced by scripts/04b_external_scores.py).
"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.align import FAMILY_REFERENCE, map_positions
from snrna_vep.coords import revcomp
from snrna_vep.features_rna import DUPLEX_PARTNER, gene_features
from snrna_vep.genes import PARALOG_GROUPS
from snrna_vep.reference import fetch_base, fetch_region
from snrna_vep.variants import enumerate_variants

FAMILY_OF_GROUP = {"U4": "U4", "U4atac": "U4", "U6": "U6", "U6atac": "U6", "U5": "U5",
                   "U2": "U2", "U12": "U2", "U1": "U1", "U11": "U1", "U7": "U7"}
UNCALLABLE = {"RNU1-1", "RNU1-2", "RNU1-3", "RNU1-4", "RNU2-1"}


def _one_gene(job):
    g, seq, partner_seq = job
    v = enumerate_variants(g, fetch_region(g.chrom, g.start, g.end), fetch_base)
    f = gene_features(v, seq, partner_seq)
    return v.merge(f, on="key")


def main():
    genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t")
    genes = genes[~genes.gene_name.isin(UNCALLABLE)].reset_index(drop=True)
    rna = {g.gene_name: (lambda s: s if g.strand == "+" else revcomp(s))(fetch_region(g.chrom, g.start, g.end))
           for g in genes.itertuples()}
    for name in set(FAMILY_REFERENCE.values()) | set(DUPLEX_PARTNER.values()):
        if name not in rna:   # family references that are uncallable still provide coordinates
            cat = pd.read_csv(config.PROCESSED / "small_rna_catalogue.tsv", sep="\t").set_index("gene_name")
            r = cat.loc[name]
            s = fetch_region(r.chrom, int(r.start), int(r.end))
            rna[name] = s if r.strand == "+" else revcomp(s)

    dep = pd.read_csv(config.RESULTS / "depletion_per_position.tsv", sep="\t")

    # ViennaRNA folding dominates the runtime; one process per gene.
    jobs = [(SimpleNamespace(**r), rna[r["gene_name"]], rna.get(DUPLEX_PARTNER.get(r["gene_name"])))
            for r in genes.to_dict("records")]
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as pool:
        per_gene = dict(zip(genes.gene_name, pool.map(_one_gene, jobs)))

    frames = []
    for g in genes.itertuples():
        v = per_gene[g.gene_name]
        partner = DUPLEX_PARTNER.get(g.gene_name)
        grp = g.paralog_group if isinstance(g.paralog_group, str) else None
        fam = FAMILY_OF_GROUP.get(grp) if grp else None
        v["paralog_family"] = fam
        if grp:
            m = map_positions(rna[g.gene_name], rna[FAMILY_REFERENCE[grp]])
            v["family_ref_pos"] = [m.get(int(np.floor(p))) for p in v.n_pos]
        else:
            v["family_ref_pos"] = np.nan
        d = dep[dep.gene_name == g.gene_name].set_index("n_pos")
        if len(d):
            w = d[["obs", "exp_gene"]].rolling(10, center=True, min_periods=5).sum()
            v["rel_oe_w10"] = v.n_pos.apply(lambda p: (w.obs / w.exp_gene).get(int(np.floor(p)), np.nan))
        frames.append(v)
        print(f"{g.gene_name:9s} {len(v):5d} variants  partner={partner}")

    feat = pd.concat(frames, ignore_index=True)

    # gnomAD observation (for labels / filtering, never as a default feature)
    gn = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")
    gn["key"] = gn.chrom + ":" + gn.pos.astype(str) + ":" + gn.ref + ":" + gn.alt
    gn = gn.drop_duplicates("key").set_index("key")
    feat["gnomad_AC"] = feat.key.map(gn.AC)
    feat["gnomad_nhomalt"] = feat.key.map(gn.nhomalt)
    feat["gnomad_pass"] = feat.key.map(gn["filter"] == "PASS")
    feat["gnomad_phylop"] = feat.key.map(gn.phylop)        # only defined where observed
    feat["gnomad_cadd"] = feat.key.map(gn.cadd_phred)

    # Structure contacts across spliceosome states
    pdb_dir = config.RAW / "pdb"
    cifs = sorted(pdb_dir.glob("*.cif*")) if pdb_dir.exists() else []
    if cifs:
        from snrna_vep.contacts import family_contact_features
        per, agg = family_contact_features(cifs, rna, FAMILY_REFERENCE, FAMILY_OF_GROUP, PARALOG_GROUPS)
        per.to_csv(config.PROCESSED / "contacts_per_structure.tsv", sep="\t", index=False)
        feat = feat.merge(agg, on=["paralog_family", "family_ref_pos"], how="left")
        print(f"contacts from {len(cifs)} structures merged")
    else:
        print("no structures in data/raw/pdb/ -> contact features skipped (run scripts/00_download_local.py)")

    ext_dir = config.DATA / "external"
    for f in sorted(ext_dir.glob("*.tsv")) if ext_dir.exists() else []:
        e = pd.read_csv(f, sep="\t")
        feat = feat.merge(e.drop_duplicates("key"), on="key", how="left")
        print(f"merged external scores: {f.name}")

    out = config.PROCESSED / "variant_features.tsv.gz"
    feat.to_csv(out, sep="\t", index=False)
    print(f"{len(feat)} variants x {feat.shape[1]} columns -> {out}")


if __name__ == "__main__":
    main()
