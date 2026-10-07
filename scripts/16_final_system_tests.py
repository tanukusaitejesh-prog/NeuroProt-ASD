"""Paired tests for the final snRNA-VEP system (v2-nested: v2 for dominant, contact+phyloP for recessive, the choices
made by nested selection in scripts/11_v2_model.py) against CADD and phyloP447, per mechanism.
Uses the held-out predictions of scripts/11_v2_model.py; within-unit percentile ranks for pooling.
Usage: SNRNA_LABELS=labels_v4.tsv python scripts/16_final_system_tests.py"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import bootstrap_auroc, metrics, paired_bootstrap_delta

P = pd.read_csv(config.result("v2_heldout_predictions.tsv.gz"), sep="\t", low_memory=False)
N = pd.read_csv(config.result("v2_nested.tsv"), sep="\t")
choice = {(u, c): s for u, c, s in zip(N.unit, N.ctx, N.chosen)}
col = [choice[(u, c)] + "_rk" for u, c in zip(P.unit, P.ctx)]
P["final"] = [P.at[i, c] for i, c in zip(P.index, col)]
rows = []
for ctx in ["dominant", "recessive", "both"]:
    Q = P if ctx == "both" else P[P.ctx == ctx]
    lo, hi = bootstrap_auroc(Q.y_path, Q.final, n=1000)
    rows.append(dict(ctx=ctx, comparison="final (rank-pooled AUROC)", delta=metrics(Q.y_path, Q.final)["auroc"],
                     ci=f"[{lo:.3f}, {hi:.3f}]", p=float("nan"), n=len(Q), n_pos=int(Q.y_path.sum())))
    for b in ["phylop447", "frac_states_protein_contact"]:
        dl, lo, hi, p = paired_bootstrap_delta(Q.y_path.values, Q.final.values, Q[b + "_rk"].values, n=2000)
        rows.append(dict(ctx=ctx, comparison=f"final - {b}", delta=dl, ci=f"[{lo:+.3f}, {hi:+.3f}]", p=p, n=len(Q), n_pos=int(Q.y_path.sum())))
    s = Q[Q.cadd_phred.notna()]
    dl, lo, hi, p = paired_bootstrap_delta(s.y_path.values, s.final.values, s.cadd_phred_rk.values, n=2000)
    rows.append(dict(ctx=ctx, comparison="final - CADD (SNVs)", delta=dl, ci=f"[{lo:+.3f}, {hi:+.3f}]", p=p, n=len(s), n_pos=int(s.y_path.sum())))
out = pd.DataFrame(rows)
out.to_csv(config.result("final_system_tests.tsv"), sep="\t", index=False)
print(out.round(3).to_string(index=False))
