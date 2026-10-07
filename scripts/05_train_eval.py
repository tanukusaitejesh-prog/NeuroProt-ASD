"""Benchmark: baselines vs snRNA-VEP, leave-one-gene-out, ablations, time split, SGE zero-shot.

Needs data/curation/patient_variants.tsv (see snrna_vep/labels.py for the schema).
Optional: data/external/sge_rnu4-2.tsv with columns key, sge_score.
Usage: python scripts/05_train_eval.py [--task dominant|recessive] [--cutoff 2025-01]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import bootstrap_auroc, paired_bootstrap_delta, score_table, spearman_vs_sge
from snrna_vep.labels import collapse_patients, load_patients, population_controls
from snrna_vep.model import DEFAULT_GROUPS, fit_predict, leave_one_gene_out, time_split
from snrna_vep.reference import fetch_base

ap = argparse.ArgumentParser()
ap.add_argument("--task", default="dominant", choices=["dominant", "recessive"])
ap.add_argument("--cutoff", default="2025-01")
ap.add_argument("--patients", default=str(config.CURATION / "patient_variants.tsv"))
ap.add_argument("--modes", default=None, help="comma-separated modes to treat as pathogenic")
ap.add_argument("--tag", default="", help="suffix for result files")
args = ap.parse_args()
MODES = args.modes.split(",") if args.modes else {"dominant": ["AD-NDD"], "recessive": ["AR-NDD"]}[args.task]
TAG = f"{args.task}{'_' + args.tag if args.tag else ''}"

if not Path(args.patients).exists():
    raise SystemExit(f"missing {args.patients}: curate it from the paper supplements first "
                     f"(template: data/curation/patient_variants.template.tsv)")

feat = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t")
genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t")
raw, errors = load_patients(args.patients, genes, fetch_base)
for e in errors:
    print("WARN", e)
pat = collapse_patients(raw)
pat = pat[pat["mode"].isin(MODES)]
missing = set(pat.key) - set(feat.key)
if missing:
    print(f"WARN {len(missing)} patient variants are not in the enumerated set (multi-nt indels?): "
          f"{sorted(pat[pat.key.isin(missing)].hgvs_n.fillna(pat.key).replace('', pd.NA).fillna(pat.key))[:10]}")

gn = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")
ctrl = population_controls(gn[gn.gene_name.isin(pat.gene_name.unique())], task=args.task,
                           patient_keys=pat.key)
data = feat[feat.gene_name.isin(pat.gene_name.unique())].copy()
data["label"] = np.where(data.key.isin(pat.key), 1, np.where(data.key.isin(ctrl.key), 0, -1))
data = data[data.label >= 0].merge(pat[["key", "first_published"]], on="key", how="left")
data["first_published"] = pd.to_datetime(data.first_published)
print(data.groupby("gene_name").label.agg(["sum", "count"]).rename(columns={"sum": "pathogenic", "count": "total"}))

# Baselines: oriented so that higher = more damaging
data["abs_ddG_fold"] = data.ddG_fold.abs()
data["ddG_duplex_loss"] = data.ddG_duplex          # positive = duplex destabilised
data["paired_wt"] = 1 - data.p_unpaired_wt
for c, sign in [("rnalm_dll", -1), ("evo2_dll", -1)]:
    if c in data:
        data[c + "_neg"] = sign * data[c]
if "rel_oe_w10" in data:
    data["depletion_neg_oe"] = -data.rel_oe_w10   # population-derived: circular with gnomAD controls
baselines = [c for c in ["phylop", "cadd_phred", "alphagenome_avi", "rnalm_dll_neg", "evo2_dll_neg",
                         "abs_ddG_fold", "ddG_duplex_loss", "paired_wt", "max_protein_res",
                         "frac_states_snrna_contact", "depletion_neg_oe"]
             if c in data and data[c].notna().any()]

out = config.RESULTS
res = []
loso = {}
for kind in ["gbm", "logistic"]:
    p = leave_one_gene_out(data, DEFAULT_GROUPS, kind)
    if len(p):
        loso[kind] = p
        data = data.merge(p[["key", "pred"]].rename(columns={"pred": f"snrnavep_{kind}"}), on="key", how="left")
scores = baselines + [f"snrnavep_{k}" for k in loso]
tab = score_table(data[data[[f"snrnavep_{k}" for k in loso][:1]].notna().any(axis=1)] if loso else data,
                  "label", scores, group_col="gene_name")
tab.to_csv(out / f"benchmark_logo_{TAG}.tsv", sep="\t", index=False)
print("\n== leave-one-gene-out ==")
print(tab[tab.group == "all"].round(3).to_string(index=False))

if "snrnavep_gbm" in data:
    d = data.dropna(subset=["snrnavep_gbm"])
    lo, hi = bootstrap_auroc(d.label, d.snrnavep_gbm)
    print(f"snRNA-VEP (gbm) AUROC 95% CI [{lo:.3f}, {hi:.3f}]")
    for b in baselines:
        delta, l, h, pval = paired_bootstrap_delta(d.label.values, d.snrnavep_gbm.values, d[b].values)
        res.append(dict(vs=b, delta_auroc=delta, ci_lo=l, ci_hi=h, p=pval))
    pd.DataFrame(res).to_csv(out / f"benchmark_paired_{TAG}.tsv", sep="\t", index=False)
    print(pd.DataFrame(res).round(4).to_string(index=False))

# Ablations: drop one feature group at a time
abl = []
for drop in DEFAULT_GROUPS:
    groups = [g for g in DEFAULT_GROUPS if g != drop]
    p = leave_one_gene_out(data, groups, "gbm")
    if len(p):
        m = score_table(p, "label", ["pred"]).query("group == 'all'").iloc[0]
        abl.append(dict(dropped=drop, auroc=m.auroc, auprc=m.auprc))
pd.DataFrame(abl).to_csv(out / f"ablation_{TAG}.tsv", sep="\t", index=False)
print("\n== ablations (LOGO, gbm) ==")
print(pd.DataFrame(abl).round(3).to_string(index=False))

# Time split
if data.first_published.notna().any():
    p, cols = time_split(data, args.cutoff)
    if p.label.sum() > 0:
        t = score_table(p.merge(data[["key"] + baselines], on="key"), "label", ["pred"] + baselines)
        t.to_csv(out / f"benchmark_timesplit_{TAG}.tsv", sep="\t", index=False)
        print(f"\n== time split (train < {args.cutoff}) ==")
        print(t.round(3).to_string(index=False))

# SGE zero-shot on RNU4-2: train on every other gene, never on SGE
sge_path = config.DATA / "external" / "sge_rnu4-2.tsv"
if sge_path.exists():
    sge = pd.read_csv(sge_path, sep="\t")
    allv = feat[feat.gene_name == "RNU4-2"].merge(sge, on="key")
    train = data[data.gene_name != "RNU4-2"]
    pred, _, _ = fit_predict(train, allv.assign(label=0), DEFAULT_GROUPS, "gbm")
    allv["snrnavep_gbm"] = pred
    allv["abs_ddG_fold"] = allv.ddG_fold.abs()
    allv["ddG_duplex_loss"] = allv.ddG_duplex
    s = spearman_vs_sge(allv, ["snrnavep_gbm"] + [b for b in baselines if b in allv])
    s.to_csv(out / "sge_zero_shot_rnu4-2.tsv", sep="\t", index=False)
    print("\n== zero-shot vs RNU4-2 SGE (Spearman with -score) ==")
    print(s.round(3).to_string(index=False))
