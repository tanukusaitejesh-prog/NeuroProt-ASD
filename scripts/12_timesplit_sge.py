"""Time split and SGE zero-shot for snRNA-VEP v2 (dominant context).

Time split: train only on what was public before 2024-09 (RNU4-2 dominant ReNU variants, first reported 2024-04,
plus RNU4-2 controls); test on later dominant discoveries in other genes: RNU2-2 (preprint 2024-09), RNU5B-1
(2024-10), RNU6 retinitis pigmentosa genes (2025-01).
SGE zero-shot: v2 (dominant context) trained on every unit except RNU4-2; Spearman with -SGE on RNU4-2.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import bootstrap_auroc, metrics

src = open(Path(__file__).resolve().parent / "11_v2_model.py").read().split("results, pooled = [], []")[0]
G = {"__file__": str(Path(__file__).resolve().parent / "11_v2_model.py"), "__name__": "ts"}
exec(compile(src, "11_v2_model", "exec"), G)
d, rows_for_training, design, model = G["d"], G["rows_for_training"], G["design"], G["model"]

# ---- time split
train = d[(d.gene_name == "RNU4-2") & ((d.y_path == 0) | d["mode"].str.contains("AD-NDD"))].copy()
train["mode"] = np.where(train.y_path == 1, "AD-NDD", "")
tr = rows_for_training(train)
m = model().fit(design(d, tr.i, ["dominant"] * len(tr)), tr.y)
print(f"time split: trained on RNU4-2 only ({int(train.y_path.sum())} ReNU variants, {int((train.y_path == 0).sum())} controls)")
rows = []
for unit in ["RNU2-2P", "RNU5B-1", "RNU6"]:
    t = d[(d.unit == unit) & ((d.y_path == 0) | ((d.y_path == 1) & d["mode"].str.contains("AD-")))].copy()
    t["v2_timesplit"] = m.predict_proba(design(d, t.index, ["dominant"] * len(t)))[:, 1]
    for s in ["v2_timesplit", "frac_states_protein_contact", "phylop447", "cadd_phred"]:
        mm = metrics(t.y_path, t[s])
        lo, hi = bootstrap_auroc(t.y_path, t[s], n=1000) if mm["n_pos"] >= 3 and t[s].notna().sum() > 10 else (np.nan, np.nan)
        rows.append(dict(unit=unit, score=s, auroc=mm["auroc"], ci=f"[{lo:.2f}, {hi:.2f}]", n_pos=mm["n_pos"]))
ts = pd.DataFrame(rows)
ts.to_csv(config.result("v2_timesplit.tsv"), sep="\t", index=False)
print(ts.pivot_table(index="score", columns="unit", values="auroc").round(3).to_string())
print(ts[ts.score == "v2_timesplit"][["unit", "auroc", "ci", "n_pos"]].to_string(index=False))

# ---- SGE zero-shot
test = d[(d.gene_name == "RNU4-2") & d.y_sge.notna()]
out = []
for name, tr_df in [("v2 trained without RNU4-2", d[d.unit != "RNU4-2"]),
                    ("v2 trained without U4 family", d[d.paralog_family != "U4"])]:
    tr = rows_for_training(tr_df)
    mm = model().fit(design(d, tr.i, tr.ctx), tr.y)
    p = mm.predict_proba(design(d, test.index, ["dominant"] * len(test)))[:, 1]
    out.append(dict(score=name, spearman=spearmanr(p, test.y_sge)[0], n=len(test)))
for s in ["nb_frac_prot", "frac_states_protein_contact", "cadd_phred", "phylop447"]:
    x = test[[s, "y_sge"]].dropna()
    out.append(dict(score=s, spearman=spearmanr(x[s], x.y_sge)[0], n=len(x)))
z = pd.DataFrame(out)
z.to_csv(config.result("v2_sge_zeroshot.tsv"), sep="\t", index=False)
print("\nSGE zero-shot (Spearman with -SGE):"); print(z.round(3).to_string(index=False))
