"""Prospective ClinVar test: could snRNA-VEP v2 have predicted variants that became P/LP after the 2025-05-04 snapshot?

Training knowledge (public before 2025-05-04): ClinVar 2025-05-04 P/LP + curated supplement variants first public before
2025-05 (Chen2024, Greene2025 preprint 2024-09, Nava2025 2024-10, Jackson2025 preprint 2024-10, RP2026 2025-01, SGE
ReNU list 2024-04) + the same mode-aware controls (gnomAD v4.1, UKB/AoU).
Test: ClinVar P/LP added after 2025-05-04 (latest release 2026-10-04) that were not in the training knowledge, in a
held-out gene (the model never sees that gene: time split AND gene split), against the gene's controls.
Mode of a new variant from ClinVar's condition text (scripts/05e cv_mode). Units need >= 3 new P/LP and >= 3 controls.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import bootstrap_auroc, metrics, paired_bootstrap_delta

src = open(Path(__file__).resolve().parent / "11_v2_model.py").read().split("results, pooled = [], []")[0]
G = {"__file__": str(Path(__file__).resolve().parent / "11_v2_model.py"), "__name__": "pros"}
exec(compile(src, "11_v2_model", "exec"), G)
d, design, model, R = G["d"], G["design"], G["model"], G["R"]


def load(f):
    c = pd.read_csv(config.PROCESSED / f, sep="\t")
    c["key"] = c.chrom + ":" + c.pos.astype(str) + ":" + c.ref + ":" + c.alt
    return c[~c.revstat.fillna("").str.contains("no_assertion|no_classification")]


PL = ["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"]
old, new = load("clinvar_20250504_core_loci.tsv"), load("clinvar_latest_core_loci.tsv")


def cv_mode(gene, disease):
    t = str(disease).lower()
    return "AD-RP" if "retinitis" in t else "AR-NDD" if ("recessive" in t or gene in ("RNU4ATAC", "RNU12")) else "AD-NDD"


pat = pd.read_csv(config.CURATION / "patient_variants_curated.tsv", sep="\t")
EARLY = ["Chen2024", "Greene2025", "Nava2025", "Jackson2025", "RP2026", "SGE2026"]
pre = pd.concat([pat[pat.source.str.split().str[0].isin(EARLY)][["key", "mode"]],
                 old[old.clnsig.isin(PL)].assign(mode=lambda x: [cv_mode(g, s) for g, s in zip(x.gene_name, x.disease)])[["key", "mode"]]])
pre_modes = pre.groupby("key")["mode"].agg(lambda s: "|".join(sorted(set(s))))
newp = new[new.clnsig.isin(PL) & ~new.key.isin(pre_modes.index)].copy()
newp["mode"] = [cv_mode(g, s) for g, s in zip(newp.gene_name, newp.disease)]
new_modes = newp.groupby("key")["mode"].agg(lambda s: "|".join(sorted(set(s))))

# training labels: pre-2025-05 pathogenic + controls (controls exclude anything later called P/LP)
ctrl = (d.y_path == 0) & ~d.key.isin(new_modes.index)
d["y_pre"] = np.where(d.key.isin(pre_modes.index), 1.0, np.where(ctrl, 0.0, np.nan))
d["mode_pre"] = d.key.map(pre_modes).fillna("")
gm_pre = d[d.y_pre == 1].groupby("gene_name")["mode_pre"].agg(lambda s: set("|".join(s).split("|")))


def rows_pre(df):
    out = []
    for r in df[df.y_pre.notna()].itertuples():
        ms = r.mode_pre.split("|") if r.y_pre == 1 else gm_pre.get(r.gene_name, set())
        ctxs = {("recessive" if m.startswith("AR") else "dominant") for m in ms if m} or {"dominant"}
        out += [(r.Index, c, r.y_pre) for c in ctxs]
    return pd.DataFrame(out, columns=["i", "ctx", "y"])


print(f"training knowledge before 2025-05-04: {int((d.y_pre == 1).sum())} pathogenic; new P/LP after: "
      f"{int(d.key.isin(new_modes.index).sum())} (with a structure node)")
rows, pooled = [], []
for unit in ["RNU4-2", "RNU2-2P", "RNU4ATAC", "RNU12", "RNU6"]:
    tr = rows_pre(d[d.unit != unit])
    m = model().fit(design(d, tr.i, tr.ctx), tr.y)
    for ctx, pref in [("dominant", "AD-"), ("recessive", "AR-")]:
        is_new = d.key.map(new_modes).fillna("").str.contains(pref)
        t = d[(d.unit == unit) & (is_new | ctrl)].copy()
        t["y"] = is_new[t.index].astype(float)
        if t.y.sum() < 3 or (t.y == 0).sum() < 3:
            continue
        t["snrnavep_v2"] = m.predict_proba(design(d, t.index, [ctx] * len(t)))[:, 1]
        t["combo_untrained"] = t[["frac_states_protein_contact_r", "phylop447_r"]].mean(1)
        for s in ["snrnavep_v2", "combo_untrained", "frac_states_protein_contact", "phylop447", "cadd_phred"]:
            mm = metrics(t.y, t[s])
            rows.append(dict(unit=unit, ctx=ctx, score=s, auroc=mm["auroc"], n_new=int(t.y.sum()), n_ctrl=int((t.y == 0).sum())))
            t[s + "_rk"] = t[s].rank(pct=True)
        pooled.append(t.assign(ctx=ctx))
res = pd.DataFrame(rows)
res.to_csv(config.result("prospective_clinvar.tsv"), sep="\t", index=False)
print(res.pivot_table(index="score", columns=["ctx", "unit"], values="auroc").round(3).to_string())
print(res.drop_duplicates(["unit", "ctx"])[["unit", "ctx", "n_new", "n_ctrl"]].to_string(index=False))
P = pd.concat(pooled)
lo, hi = bootstrap_auroc(P.y, P.snrnavep_v2_rk, n=1000)
print(f"\nrank-pooled v2 AUROC {metrics(P.y, P.snrnavep_v2_rk)['auroc']:.3f} [{lo:.3f}, {hi:.3f}]  (new P/LP {int(P.y.sum())})")
for b in ["combo_untrained", "phylop447", "frac_states_protein_contact"]:
    dl, l, h, p = paired_bootstrap_delta(P.y.values, P.snrnavep_v2_rk.values, P[b + "_rk"].values, n=2000)
    print(f"  v2 - {b:28s}: {dl:+.3f} [{l:+.3f}, {h:+.3f}] p {p:.3f}")
s = P[P.cadd_phred.notna()]
dl, l, h, p = paired_bootstrap_delta(s.y.values, s.snrnavep_v2_rk.values, s.cadd_phred_rk.values, n=2000)
print(f"  v2 - CADD (SNVs n {len(s)}, new {int(s.y.sum())}): {dl:+.3f} [{l:+.3f}, {h:+.3f}] p {p:.3f}")
