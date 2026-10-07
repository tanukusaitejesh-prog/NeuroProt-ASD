"""Run on a machine with normal internet: downloads spliceosome structures from the PDB.

Supplementary tables (patient variants, SGE scores) must be fetched by hand from the papers listed
in PROPOSAL_snRNA-VEP.md and curated into data/curation/patient_variants.tsv and data/external/sge_rnu4-2.tsv.
"""
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.contacts import SPLICEOSOME_STATES

out = config.RAW / "pdb"
out.mkdir(parents=True, exist_ok=True)
for pdb_id, state in SPLICEOSOME_STATES.items():
    dest = out / f"{pdb_id}.cif.gz"
    if dest.exists():
        continue
    mirrors = [f"https://files.rcsb.org/download/{pdb_id}.cif.gz",
               "https://pdbsnapshots.s3.us-west-2.amazonaws.com/20260101/pub/pdb/data/structures/divided/"
               f"mmCIF/{pdb_id.lower()[1:3]}/{pdb_id.lower()}.cif.gz"]
    for url in mirrors:
        try:
            urllib.request.urlretrieve(url, dest)
            print(f"{pdb_id} ({state}) ok from {url.split('/')[2]}")
            break
        except Exception as e:
            print(f"{pdb_id} ({state}) failed from {url.split('/')[2]}: {e}")
