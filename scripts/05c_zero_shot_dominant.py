"""Zero-shot check on a dominant gene (default RNU4-2): ClinVar P/LP vs gnomAD-observed variants.

No training on the gene. Each feature is scored as-is (oriented so higher = more damaging), and a model
trained on every *other* labelled gene is applied as a transfer test. Depletion is population-derived,
so it is shown but marked circular (controls are gnomAD variants).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import bootstrap_auroc, metrics
from snrna_vep.labels import load_patients, population_controls
from snrna_vep.model import DEFAULT_GROUPS, fit_predict
from snrna_vep.reference import fetch_base

ap = argparse.ArgumentParser()
ap.add_argument("--gene", default="RNU4-2")
ap.add_argument("--patients", default=str(config.CURATION / "patient_variants_clinvar.tsv"))
args = ap.parse_args()

feat = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t")
genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t")
pat, _ = load_patients(args.patients, genes, fetch_base)
gn = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")

pos = set(pat[pat.gene_name == args.gene].key)
ctrl = population_controls(gn[gn.gene_name == args.gene], "dominant", patient_keys=pos)
d = feat[feat.gene_name == args.gene].copy()
d["label"] = np.where(d.key.isin(pos), 1, np.where(d.key.isin(ctrl.key), 0, -1))
d = d[d.label >= 0]
print(f"{args.gene}: {int(d.label.sum())} pathogenic vs {int((d.label == 0).sum())} gnomAD variants")

# transfer model: trained on all other labelled genes (their pathogenic vs their controls)
others = pat[pat.gene_name != args.gene]
tr = feat[feat.gene_name.isin(others.gene_name.unique())].copy()
ctrl_o = pd.concat([population_controls(gn[gn.gene_name == g], "recessive" if m.startswith("AR") else "dominant",
                                        patient_keys=others.key)
                    for g, m in others.groupby("gene_name")["mode"].first().items()])
tr["label"] = np.where(tr.key.isin(others.key), 1, np.where(tr.key.isin(ctrl_o.key), 0, -1))
tr = tr[tr.label >= 0]
d["transfer_model"], cols, _ = fit_predict(tr, d, [g for g in DEFAULT_GROUPS if g != "paralog"], "logistic")

d["abs_ddG_fold"] = d.ddG_fold.abs()
d["ddG_duplex_loss"] = d.ddG_duplex
d["paired_wt"] = 1 - d.p_unpaired_wt
d["depletion_neg_oe (circular)"] = -d.rel_oe_w10
scores = ["phylop_gnomadpos", "max_protein_res", "frac_states_protein_contact", "frac_states_snrna_contact",
          "abs_ddG_fold", "ddG_duplex_loss", "paired_wt", "transfer_model", "depletion_neg_oe (circular)"]
rows = []
for s in scores:
    if s in d and d[s].notna().any():
        m = metrics(d.label, d[s])
        lo, hi = bootstrap_auroc(d.label, d[s], n=500) if m["n_pos"] >= 3 else (np.nan, np.nan)
        rows.append(dict(score=s, **m, auroc_ci=f"[{lo:.2f}, {hi:.2f}]"))
out = pd.DataFrame(rows)
out.to_csv(config.RESULTS / f"zero_shot_{args.gene}_clinvar.tsv", sep="\t", index=False)
print(out.round(3).to_string(index=False))
