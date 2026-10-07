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
    url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    try:
        urllib.request.urlretrieve(url, dest)
        print(f"{pdb_id} ({state}) ok")
    except Exception as e:
        print(f"{pdb_id} ({state}) FAILED: {e}  -> check the ID on rcsb.org")
