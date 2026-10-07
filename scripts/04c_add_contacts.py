"""Add / refresh cryo-EM contact features on an existing feature table (no re-folding).

Structures: data/raw/pdb/*.cif.gz (scripts/00_download_local.py, or the AWS PDB snapshot mirror).
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.align import FAMILY_REFERENCE
from snrna_vep.contacts import family_contact_features
from snrna_vep.coords import revcomp
from snrna_vep.genes import PARALOG_GROUPS
from snrna_vep.model import FEATURE_GROUPS
from snrna_vep.reference import fetch_region

FAMILY_OF_GROUP = {"U4": "U4", "U4atac": "U4", "U6": "U6", "U6atac": "U6", "U5": "U5",
                   "U2": "U2", "U12": "U2", "U1": "U1", "U11": "U1", "U7": "U7"}

cat = pd.read_csv(config.PROCESSED / "small_rna_catalogue.tsv", sep="\t")
members = {g for gs in PARALOG_GROUPS.values() for g in gs}
cat = cat[cat.gene_name.isin(members)].drop_duplicates("gene_name")
rna = {g.gene_name: (lambda s: s if g.strand == "+" else revcomp(s))(fetch_region(g.chrom, int(g.start), int(g.end)))
       for g in cat.itertuples()}

cifs = sorted((config.RAW / "pdb").glob("*.cif*"))
per, agg = family_contact_features(cifs, rna, FAMILY_REFERENCE, FAMILY_OF_GROUP, PARALOG_GROUPS)
per.to_csv(config.PROCESSED / "contacts_per_structure.tsv", sep="\t", index=False)
print(per.groupby(["pdb_id", "state", "family_gene"]).size().rename("nucleotides").to_string())

path = config.PROCESSED / "variant_features.tsv.gz"
feat = pd.read_csv(path, sep="\t")
feat = feat.drop(columns=[c for c in FEATURE_GROUPS["contacts"] if c in feat])
feat = feat.merge(agg, on=["paralog_family", "family_ref_pos"], how="left")
feat.to_csv(path, sep="\t", index=False)
print(f"contact features on {feat.max_protein_res.notna().sum()} / {len(feat)} variants")
