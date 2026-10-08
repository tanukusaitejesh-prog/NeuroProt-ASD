"""Does adding allele-resolved features improve the held-out predictor, or only within-position ranking?

scripts/44 showed the allele features rank alternative alleles correctly WITHIN a position (mean rho +0.346
against -SGE, where CADD reaches -0.077). That is a different question from whether they improve discrimination
of pathogenic from benign variants ACROSS positions, which is what the clinical model does.

This script answers it under the identical held-out-unit protocol as the main model (scripts/11, 16), adding the
allele features to the frozen v2 design:
    frac_states_pair_broken_by_alt, frac_states_wobble_by_alt, frac_states_base_contact, transition
Base-versus-backbone contact is included because only base contacts can be perturbed by a substitution, so it
carries allele-relevant information the position-level contact fraction does not.

STATUS: the allele features themselves were pre-specified (ANALYSIS_PLAN.md section 3). Adding them to the
clinical model and re-evaluating held-out AUROC was NOT pre-specified and is reported as a secondary,
exploratory analysis. The frozen v2 result stands as primary regardless of the outcome here.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
src = open(HERE / "11_v2_model.py").read().split("results, pooled = [], []")[0]
G = {"__file__": str(HERE / "11_v2_model.py"), "__name__": "am"}
exec(compile(src, "11_v2_model", "exec"), G)
config, d, rows_for_training, design, model, R, metrics = (
    G["config"], G["d"], G["rows_for_training"], G["design"], G["model"], G["R"], G["metrics"])
from snrna_vep.evaluate import paired_bootstrap_delta

# ---- attach allele features (SNVs only; indels get the position-level values, see Methods)
A = pd.read_csv(config.PROCESSED / "allele_features.tsv.gz", sep="\t")
A = A.rename(columns={"family_ref_pos": "fpos"})
A["fpos"] = A.fpos.astype(float)
d["alt_rna"] = d.alt.str.replace("T", "U", regex=False)
key = A.set_index(["paralog_family", "fpos", "alt"])
ALLELE = ["frac_states_pair_broken_by_alt", "frac_states_wobble_by_alt", "frac_states_base_contact", "transition"]
for c in ALLELE:
    m = key[c].to_dict()
    d[c] = [m.get((f, p, a), np.nan) for f, p, a in zip(d.paralog_family, d.family_ref_pos.astype(float), d.alt_rna)]
cov = d[ALLELE[0]].notna()
print(f"variants with allele features: {int(cov.sum())} of {len(d)} ({100*cov.mean():.0f}%)")
print(f"  among SNVs: {int((cov & (d.vtype=='snv')).sum())} of {int((d.vtype=='snv').sum())}")
# neutral fill so every variant stays scoreable; within-gene ranks as for the other features
for c in ALLELE:
    d[c] = d[c].fillna(d[c].median())
    d[c + "_r"] = d.groupby("gene_name")[c].rank(pct=True)

BASE = list(R)
AUG = BASE + [c + "_r" for c in ALLELE]
UNITS = ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU6", "RNU4ATAC", "RNU12"]


def evaluate(cols, tag):
    R[:] = cols
    rows, pooled = [], []
    for unit in UNITS:
        tr = rows_for_training(d[d.unit != unit])
        m = model().fit(design(d, tr.i, tr.ctx), tr.y)
        for ctx, pref in [("dominant", "AD-"), ("recessive", "AR-")]:
            t = d[(d.unit == unit) & ((d.y_path == 0) | ((d.y_path == 1) & d["mode"].str.contains(pref)))].copy()
            if (t.y_path == 1).sum() < 3 or (t.y_path == 0).sum() < 3:
                continue
            t["score"] = m.predict_proba(design(d, t.index, [ctx] * len(t)))[:, 1]
            rows.append(dict(model=tag, unit=unit, ctx=ctx, auroc=metrics(t.y_path, t.score)["auroc"],
                             n_pos=int(t.y_path.sum()), n=len(t)))
            t["rk"] = t.score.rank(pct=True)
            pooled.append(t.assign(ctx=ctx))
    R[:] = BASE
    return pd.DataFrame(rows), pd.concat(pooled)


b_rows, b_pool = evaluate(BASE, "v2 (frozen)")
a_rows, a_pool = evaluate(AUG, "v2 + allele features")
RES = pd.concat([b_rows, a_rows])
print("\n=== mean held-out AUROC ===")
print(RES.groupby(["model", "ctx"]).auroc.mean().unstack().round(3).to_string())
print("\n=== per unit ===")
print(RES.pivot_table(index=["ctx", "unit"], columns="model", values="auroc").round(3).to_string())

print("\n=== paired bootstrap, augmented minus frozen (rank-pooled within unit) ===")
for ctx in ["dominant", "recessive"]:
    b = b_pool[b_pool.ctx == ctx].set_index("key"); a = a_pool[a_pool.ctx == ctx].set_index("key")
    common = b.index.intersection(a.index)
    dl, lo, hi, p = paired_bootstrap_delta(b.loc[common].y_path.values, a.loc[common].rk.values,
                                           b.loc[common].rk.values, n=2000)
    print(f"  {ctx:9s} delta {dl:+.3f} [{lo:+.3f}, {hi:+.3f}] p {p:.3f}  (n={len(common)}, pathogenic={int(b.loc[common].y_path.sum())})")

RES.to_csv(config.result("allele_model_heldout.tsv"), sep="\t", index=False)
print("\nNOTE: exploratory. The frozen v2 result remains primary.")
