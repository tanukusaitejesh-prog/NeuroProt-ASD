"""Main benchmark on curated labels (scripts/05b_curate_supplements.py).

Tasks
  dominant : AD-NDD variants vs controls, held-out gene (RNU4-2, RNU2-2, RNU5B-1).
  anypath  : any pathogenic mode (AD-NDD, AR-NDD, AD-RP) vs controls, held-out gene.
  sge      : model trained without any RNU4-2 label, Spearman with -SGE function score on RNU4-2.
Controls  : RNU4-2 -> UK Biobank / All of Us variants (gnomAD-independent); RP-benign variants;
            other genes -> gnomAD v4.1 PASS AC>=1 (population-derived features are never used).
Output    : results/benchmark_main_*.tsv, printed summary.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import bootstrap_auroc, metrics, paired_bootstrap_delta
from snrna_vep.model import FEATURE_GROUPS, fit_predict

FEATURE_GROUPS["conservation"] = ["phylop447"]          # true track only; no gnomAD-annotated phyloP
GROUPS = ["position", "structure2d", "contacts", "conservation", "paralog"]
SEED = 0

feat = pd.read_csv(config.PROCESSED / "variant_features.tsv.gz", sep="\t", low_memory=False)
pat = pd.read_csv(config.CURATION / "patient_variants_curated.tsv", sep="\t")
ctl = pd.read_csv(config.CURATION / "controls_curated.tsv", sep="\t")
gn = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")
sge = pd.read_csv(config.DATA / "external" / "sge_rnu4-2.tsv", sep="\t")

pat_first = pat.groupby(["key", "mode"]).first_published.min().reset_index()
path_keys = set(pat.key)

# controls
gnc = gn[(gn["filter"] == "PASS") & ~gn.segdup.astype(bool) & (gn.AC >= 1)]
gnc = gnc.assign(key=gnc.chrom + ":" + gnc.pos.astype(str) + ":" + gnc.ref + ":" + gnc.alt)
ctrl_keys = {}
for g in feat.gene_name.unique():
    if g == "RNU4-2":
        ks = set(ctl[(ctl.gene_name == g) & (ctl.control_type == "UKB/AoU")].key)
    else:
        ks = set(gnc[gnc.gene_name == g].key)
    ks |= set(ctl[(ctl.gene_name == g) & (ctl.control_type == "RP benign")].key)
    ctrl_keys[g] = ks - path_keys


def assemble(modes):
    pk = set(pat[pat["mode"].isin(modes)].key)
    genes = sorted(pat[pat["mode"].isin(modes)].gene_name.unique()) if modes != "all" else None
    d = feat.copy()
    d["label"] = np.where(d.key.isin(pk), 1, np.where([k in ctrl_keys.get(g, ()) for g, k in zip(d.gene_name, d.key)], 0, -1))
    d = d[d.label >= 0]
    return d.drop_duplicates("key")


def oriented(d):
    d = d.copy()
    d["abs_ddG_fold"] = d.ddG_fold.abs()
    d["paired_wt"] = 1 - d.p_unpaired_wt
    return d


BASE = ["cadd_phred", "phylop447", "max_protein_res", "frac_states_protein_contact", "frac_states_snrna_contact",
        "abs_ddG_fold", "paired_wt"]


def logo(d, groups, kind, eval_genes, train_genes=None):
    out = []
    for g in eval_genes:
        test = d[d.gene_name == g]
        train = d[(d.gene_name != g) & (d.gene_name.isin(train_genes) if train_genes is not None else True)]
        if test.label.sum() < 3 or test.label.nunique() < 2 or train.label.sum() < 3:
            continue
        p, cols, _ = fit_predict(train, test, groups, kind)
        out.append(test[["gene_name", "key", "label"]].assign(pred=p))
    return pd.concat(out) if out else pd.DataFrame()


def report(name, d, preds):
    rows = []
    d = oriented(d)
    for g in sorted(preds.gene_name.unique()):
        t = d[d.gene_name == g].merge(preds[["key"] + [c for c in preds.columns if c.startswith("model_")]], on="key")
        for s in BASE + [c for c in t.columns if c.startswith("model_")]:
            if s not in t or t[s].notna().sum() < 10:
                continue
            m = metrics(t.label, t[s])
            if m["n_pos"] < 3:
                continue
            lo, hi = bootstrap_auroc(t.label, t[s], n=500) if np.isfinite(m["auroc"]) else (np.nan, np.nan)
            rows.append(dict(task=name, gene=g, score=s, **m, ci_lo=lo, ci_hi=hi))
    r = pd.DataFrame(rows)
    r.to_csv(config.RESULTS / f"benchmark_main_{name}.tsv", sep="\t", index=False)
    piv = r.pivot_table(index="score", columns="gene", values="auroc").round(3)
    npos = r.groupby("gene").n_pos.max()
    print(f"\n=== {name}: AUROC by held-out gene (n pathogenic: {npos.to_dict()}) ===")
    print(piv.to_string())
    return r


def run(name, modes, eval_genes, train_all=True):
    d = assemble(modes)
    preds = None
    for kind in ["logistic", "gbm"]:
        for gname, groups in [("full", GROUPS), ("contacts_only", ["contacts"]), ("no_contacts", [g for g in GROUPS if g != "contacts"])]:
            p = logo(d, groups, kind, eval_genes)
            if p.empty:
                continue
            col = f"model_{kind}_{gname}"
            p = p.rename(columns={"pred": col})
            preds = p if preds is None else preds.merge(p[["key", col]], on="key", how="outer")
            preds["gene_name"] = preds["gene_name"].fillna(preds.key.map(d.set_index("key").gene_name))
            preds["label"] = preds["label"].fillna(preds.key.map(d.set_index("key").label))
    r = report(name, d, preds)
    # paired: best model vs CADD on SNVs CADD scores, pooled over held-out genes
    t = oriented(d).merge(preds.drop(columns=["gene_name", "label"]), on="key")
    snv = t[t.cadd_phred.notna()]
    print(f"paired vs CADD (SNVs only, pooled held-out genes; n {len(snv)}, pathogenic {int(snv.label.sum())}):")
    for m in [c for c in t.columns if c.startswith("model_")] + ["max_protein_res", "phylop447"]:
        if snv[m].notna().sum() > 10 and snv.label.sum() >= 5:
            dlt, lo, hi, p = paired_bootstrap_delta(snv.label.values, snv[m].values, snv.cadd_phred.values, n=1000)
            print(f"  {m:32s} - CADD: dAUROC {dlt:+.3f} [{lo:+.3f}, {hi:+.3f}] p {p:.3f}")
    return d, preds, r


if __name__ == "__main__":
    config.RESULTS.mkdir(exist_ok=True)
    np.random.seed(SEED)
    if "--sge-only" not in sys.argv:
        run("dominant", ["AD-NDD"], ["RNU4-2", "RNU2-2P", "RNU5B-1"])
        d, preds, r = run("anypath", ["AD-NDD", "AR-NDD", "AD-RP"], ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6-9", "RNU6-1", "RNU6-2", "RNU6-8"])

    # SGE zero-shot: train on every labelled gene except RNU4-2, predict all RNU4-2 variants
    print("\n=== RNU4-2 SGE zero-shot (Spearman with -function score; higher = more damaging) ===")
    da = assemble(["AD-NDD", "AR-NDD", "AD-RP"])
    train = da[da.gene_name != "RNU4-2"]
    allv = feat[feat.gene_name == "RNU4-2"]
    allv = oriented(allv if "sge_score" in allv else allv.merge(sge, on="key"))
    allv = allv[allv.sge_score.notna()].copy()
    rows = []
    for kind in ["logistic", "gbm"]:
        p, _, _ = fit_predict(train, allv.assign(label=0), GROUPS, kind)
        allv[f"model_{kind}"] = p
    for s in BASE + ["model_logistic", "model_gbm"]:
        x = allv[[s, "sge_score"]].dropna()
        rho, pv = spearmanr(x[s], -x.sge_score)
        rows.append(dict(score=s, n=len(x), spearman=rho, p=pv))
    sg = pd.DataFrame(rows)
    sg.to_csv(config.RESULTS / "benchmark_main_sge_zeroshot.tsv", sep="\t", index=False)
    print(sg.round(3).to_string(index=False))
