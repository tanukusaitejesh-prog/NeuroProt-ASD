"""Featurize patient variants missing from the enumerated table (multi-nt indels, complex edits).

Usage: python scripts/04d_add_patient_variants.py data/curation/patient_variants*.tsv
Variants reaching outside the transcribed region (e.g. RMRP promoter duplications) are reported and
skipped: the features describe the RNA, which they do not change directly.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.align import FAMILY_REFERENCE, map_positions
from snrna_vep.coords import revcomp
from snrna_vep.features_rna import DUPLEX_PARTNER, featurize_keys
from snrna_vep.labels import load_patients
from snrna_vep.model import FEATURE_GROUPS
from snrna_vep.reference import fetch_base, fetch_region

FAMILY_OF_GROUP = {"U4": "U4", "U4atac": "U4", "U6": "U6", "U6atac": "U6", "U5": "U5",
                   "U2": "U2", "U12": "U2", "U1": "U1", "U11": "U1", "U7": "U7"}

path = config.PROCESSED / "variant_features.tsv.gz"
feat = pd.read_csv(path, sep="\t")
genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t")
cat = pd.read_csv(config.PROCESSED / "small_rna_catalogue.tsv", sep="\t").drop_duplicates("gene_name").set_index("gene_name")


def rna_of(name):
    r = cat.loc[name]
    s = fetch_region(r.chrom, int(r.start), int(r.end))
    return s if r.strand == "+" else revcomp(s)


keys = set()
for f in sys.argv[1:]:
    df, errors = load_patients(f, genes, fetch_base)
    keys |= set(df.dropna(subset=["key"]).key) - set(feat.key)
print(f"{len(keys)} patient variants missing from the feature table")

new, skipped_all = [], []
gene_rows = genes.set_index("gene_name")
dep = pd.read_csv(config.RESULTS / "depletion_per_position.tsv", sep="\t")
contacts_cols = [c for c in FEATURE_GROUPS["contacts"] if c in feat]
contact_table = feat.dropna(subset=["family_ref_pos"])[["paralog_family", "family_ref_pos"] + contacts_cols] \
    .drop_duplicates(["paralog_family", "family_ref_pos"])
for g in genes.itertuples():
    gk = [k for k in keys if k.split(":")[0] == g.chrom and g.start - 50 <= int(k.split(":")[1]) <= g.end + 50]
    if not gk:
        continue
    partner = DUPLEX_PARTNER.get(g.gene_name)
    f, skipped = featurize_keys(gk, g, fetch_region(g.chrom, g.start, g.end), rna_of(partner) if partner else None)
    skipped_all += skipped
    if f.empty:
        continue
    grp = g.paralog_group if isinstance(g.paralog_group, str) else None
    f["paralog_family"] = FAMILY_OF_GROUP.get(grp) if grp else None
    if grp:
        m = map_positions(rna_of(g.gene_name), rna_of(FAMILY_REFERENCE[grp]))
        f["family_ref_pos"] = [m.get(int(p)) for p in f.n_pos]
    else:
        f["family_ref_pos"] = np.nan
    d = dep[dep.gene_name == g.gene_name].set_index("n_pos")
    if len(d):
        w = d[["obs", "exp_gene"]].rolling(10, center=True, min_periods=5).sum()
        oe = w.obs / w.exp_gene
        f["rel_oe_w10"] = [oe.get(int(p), np.nan) for p in f.n_pos]
    f["chrom"] = g.chrom
    f["pos"] = f.key.str.split(":").str[1].astype(int)
    f["ref"] = f.key.str.split(":").str[2]
    f["alt"] = f.key.str.split(":").str[3]
    f["hgvs"] = ""
    f["enumerated"] = False
    new.append(f)
    print(f"{g.gene_name:9s} +{len(f)} featurized, {len(skipped)} outside transcript")

if new:
    add = pd.concat(new, ignore_index=True)
    if contacts_cols:
        add = add.merge(contact_table, on=["paralog_family", "family_ref_pos"], how="left")
    if "enumerated" not in feat:
        feat["enumerated"] = True
    feat = pd.concat([feat, add], ignore_index=True)
    feat.to_csv(path, sep="\t", index=False)
print(f"skipped (outside transcript): {len(skipped_all)}; table now {len(feat)} variants")
