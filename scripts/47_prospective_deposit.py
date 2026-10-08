"""Build the timestamped prospective prediction set (ANALYSIS_PLAN.md section 4, pre-specified, previously unrun).

Saturation genome editing screens of RNU2-2 and RNU5B-1 are reported as ongoing. Predictions deposited BEFORE
those data exist can be compared against them later; predictions made afterwards cannot. The scientific value of
this file is entirely a function of its timestamp, so it is written once, hashed, and committed unchanged.

Contents: every possible single-nucleotide variant in RNU2-2 and RNU5B-1 with the model's dominant and recessive
scores and within-gene percentiles, plus a declared evaluation protocol fixed in advance so the later comparison
cannot be chosen to flatter the result.

The deposit also records predictions for RNU5A-1 (a candidate gene with VUS only) as a secondary set.
"""
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

TARGETS = {"RNU2-2P": "primary", "RNU5B-1": "primary", "RNU5A-1": "secondary"}
OUT = config.RESULTS / "prospective_deposit_2026-10-09.tsv"

R = pd.read_csv(config.RESULTS / "snrnavep_scores.tsv.gz", sep="\t")
D = R[R.gene_name.isin(TARGETS)].copy()
D["deposit_set"] = D.gene_name.map(TARGETS)
D = D[["gene_name", "deposit_set", "chrom", "pos", "ref", "alt", "key", "hgvs", "vtype",
       "score_dominant", "score_dominant_gene_pct", "score_recessive", "score_recessive_gene_pct",
       "frac_states_protein_contact", "nb_frac_prot", "phylop447", "cadd_phred", "training_label"]]
D = D.sort_values(["gene_name", "pos", "alt"]).reset_index(drop=True)
D.to_csv(OUT, sep="\t", index=False)

h = hashlib.sha256(OUT.read_bytes()).hexdigest()
meta = {
    "deposited_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "file": OUT.name,
    "sha256": h,
    "n_variants": int(len(D)),
    "by_gene": {k: int(v) for k, v in D.gene_name.value_counts().items()},
    "model": "snRNA-VEP v2-nested, label set v5c, frozen in ANALYSIS_PLAN.md",
    "primary_targets": ["RNU2-2P", "RNU5B-1"],
    "secondary_target": ["RNU5A-1"],
    "declared_evaluation_protocol": {
        "primary_metric": "Spearman rho between score_dominant and the negated SGE function score, "
                          "computed over all assayed SNVs with a structure node",
        "comparators": ["cadd_phred", "phylop447", "frac_states_protein_contact"],
        "success_criterion": "score_dominant exceeds BOTH cadd_phred and phylop447, each with a "
                             "position-block bootstrap 95% CI on the paired difference excluding zero",
        "no_refitting": "the model is not retrained, rescored or reselected after the screen is released; "
                        "the scores in this file are final",
        "reporting": "the comparison is reported whatever its outcome, including failure"
    },
    "caveats": [
        "RNU2-2P and RNU5B-1 contributed to model training in the frozen v2-nested system; these are therefore "
        "NOT zero-shot predictions. A zero-shot variant of the comparison is also pre-declared below.",
        "zero_shot_addendum: held-out scores for these genes, produced with the gene excluded from training, "
        "are in results/v2_heldout_predictions_labels_v5c.tsv.gz and should be used for the stricter test.",
        "SGE measures a cell-fitness readout, not the clinical phenotype."
    ],
}
(config.RESULTS / "prospective_deposit_2026-10-09.json").write_text(json.dumps(meta, indent=2))
print(json.dumps(meta, indent=2))
print(f"\nwrote {OUT.name}: {len(D)} variants")
print(D.groupby(["gene_name", "deposit_set"]).size().to_string())
