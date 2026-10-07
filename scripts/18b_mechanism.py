"""Mechanistic analysis: which splicing states and which named proteins distinguish pathogenic snRNA variants
from controls, separately for dominant and recessive mechanisms.

(a) State: for each gene unit and mechanism, AUROC of 'protein contact at this position in state s' among variants
    whose position is resolved in state s.
(b) Protein: Mantel-Haenszel odds ratio (strata = gene units) of contacting protein P in any state, pathogenic vs
    control; Benjamini-Hochberg FDR. Variant level (several variants share a position) -> also position level.
Usage: SNRNA_LABELS=labels_v5.tsv python scripts/18b_mechanism.py
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.evaluate import metrics

src = open(Path(__file__).resolve().parent / "11_v2_model.py").read().split("results, pooled = [], []")[0]
G = {"__file__": str(Path(__file__).resolve().parent / "11_v2_model.py"), "__name__": "mech"}
exec(compile(src, "11_v2_model", "exec"), G)
d = G["d"]

SYMBOL = [(r"SPLICING FACTOR 8|HPRP8", "PRPF8"), (r"SM D1|SM D2|SM D3|RIBONUCLEOPROTEIN [EFG]\b|PROTEINS B AND B'", "Sm ring"),
          (r"PRP31", "PRPF31"), (r"PRP3\b", "PRPF3"), (r"NHP2-LIKE", "SNU13"), (r"RNA-BINDING PROTEIN 42", "RBM42"),
          (r"CELL DIVISION CYCLE 5", "CDC5L"), (r"BUD31", "BUD31"), (r"CWC15", "CWC15"), (r"PROCESSING FACTOR 6\b", "PRPF6"),
          (r"SNW DOMAIN|^SKIP$", "SNW1"), (r"CROOKED NECK", "CRNKL1"), (r"SYF2", "SYF2"), (r"PROCESSING FACTOR 17", "CDC40"),
          (r"RBM22", "RBM22"), (r"40 KDA", "SNRNP40"), (r"3B SUBUNIT 1\b", "SF3B1"), (r"3B SUBUNIT 2\b", "SF3B2"),
          (r"3A SUBUNIT", "SF3A"), (r"PRKR-INTERACTING", "PRKRIP1"), (r"PLEIOTROPIC", "PLRG1"),
          (r"TRI-SNRNP-ASSOCIATED PROTEIN 1", "SART1"), (r"U11/U12", "U11/U12 proteins"), (r"REPETITIVE MATRIX PROTEIN 2", "SRRM2"),
          (r"U1 SMALL NUCLEAR RIBONUCLEOPROTEIN|^U1 ", "U1 snRNP proteins"), (r"27 KDA", "SNRNP27"), (r"116 KDA", "EFTUD2"),
          (r"DDX23", "DDX23"), (r"PPIL2", "PPIL2"), (r"200 KDA", "SNRNP200"), (r"SODIUM CHANNEL MODIFIER", "SCNM1"),
          (r"PROGRAMMED CELL DEATH PROTEIN 7", "PDCD7"), (r"PHD FINGER-LIKE", "PHF5A"), (r"PRP4 HOMOLOG", "PRPF4B"),
          (r"SM-LIKE PROTEIN LSM", "LSm ring"), (r"U2 SMALL NUCLEAR RIBONUCLEOPROTEIN A'", "SNRPA1"),
          (r"B''", "SNRPB2"), (r"CENTROSOMAL AT-AC", "CENATAC"), (r"THIOREDOXIN-LIKE PROTEIN 4A", "TXNL4A"),
          (r"HIV TAT-SPECIFIC", "HTATSF1"), (r"RNA-BINDING PROTEIN 48", "RBM48")]


def symbol(name):
    for pat, s in SYMBOL:
        if re.search(pat, name):
            return s
    return name.title()


n = pd.read_csv(config.PROCESSED / "graph_nodes_named.tsv.gz", sep="\t")
n["node"] = n.family + ":" + n.ref_pos.astype(int).astype(str)
n["prot_list"] = n.partners.fillna("").str.split(";").apply(lambda xs: sorted({symbol(x) for x in xs if x}))
prot_any = n.explode("prot_list").dropna(subset=["prot_list"]).groupby("node").prot_list.agg(set)
state_contact = n.assign(c=(n.n_protein_atoms > 0).astype(float)).pivot_table(index="node", columns="state", values="c", aggfunc="max")


def ctx_rows(unit, ctx):
    pref = "AD-" if ctx == "dominant" else "AR-"
    return d[(d.unit == unit) & ((d.y_path == 0) | ((d.y_path == 1) & d["mode"].str.contains(pref)))]


SETS = [(u, c) for u, c in [("RNU4-2", "dominant"), ("RNU2-2P", "dominant"), ("RNU5B-1", "dominant"), ("RNU6", "dominant"),
                            ("RNU4-2", "recessive"), ("RNU2-2P", "recessive"), ("RNU4ATAC", "recessive"), ("RNU12", "recessive")]]
# (a) states
rows = []
for u, c in SETS:
    t = ctx_rows(u, c)
    if (t.y_path == 1).sum() < 3 or (t.y_path == 0).sum() < 3:
        continue
    for s in state_contact.columns:
        x = t.node_name.map(state_contact[s])
        ok = x.notna()
        if ok.sum() < 10 or t.y_path[ok].sum() < 3 or (t.y_path[ok] == 0).sum() < 3:
            continue
        rows.append(dict(unit=u, ctx=c, state=s, auroc=metrics(t.y_path[ok], x[ok])["auroc"], n_pos=int(t.y_path[ok].sum()),
                         frac_path_contact=x[ok][t.y_path[ok] == 1].mean(), frac_ctrl_contact=x[ok][t.y_path[ok] == 0].mean()))
S = pd.DataFrame(rows)
S.to_csv(config.result("mechanism_states.tsv"), sep="\t", index=False)
print("=== (a) AUROC of protein contact in each splicing state ===")
print(S.pivot_table(index="state", columns=["ctx", "unit"], values="auroc").round(2).to_string())

# (b) proteins, MH across units, variant and position level
allprot = sorted({p for s in prot_any for p in s})
out = []
for ctx in ["dominant", "recessive"]:
    for level in ["variant", "position"]:
        for p in allprot:
            tables, npc = [], 0
            for u, c in SETS:
                if c != ctx:
                    continue
                t = ctx_rows(u, c)
                if level == "position":
                    t = t.groupby("node_name").y_path.max().rename("y_path").reset_index()
                hit = t.node_name.map(lambda nd: p in prot_any.get(nd, set()))
                a, b = int(((t.y_path == 1) & hit).sum()), int(((t.y_path == 1) & ~hit).sum())
                cc, dd = int(((t.y_path == 0) & hit).sum()), int(((t.y_path == 0) & ~hit).sum())
                if a + b >= 3 and cc + dd >= 3:
                    tables.append(np.array([[a, b], [cc, dd]], dtype=float) + 0.5)
                    npc += a
            if len(tables) == 0 or npc < 3:
                continue
            st = StratifiedTable(tables)
            out.append(dict(ctx=ctx, level=level, protein=p, OR_MH=st.oddsratio_pooled, p=st.test_null_odds().pvalue,
                            n_path_contact=npc, n_units=len(tables)))
P = pd.DataFrame(out)
P["fdr"] = P.groupby(["ctx", "level"]).p.transform(lambda x: multipletests(x, method="fdr_bh")[1])
P.to_csv(config.result("mechanism_proteins.tsv"), sep="\t", index=False)
print("\n=== (b) proteins whose contact is enriched (MH OR, strata = genes) ===")
for (ctx, level), g in P.groupby(["ctx", "level"]):
    print(f"\n[{ctx}, {level} level]")
    print(g.sort_values("p").head(12)[["protein", "OR_MH", "p", "fdr", "n_path_contact", "n_units"]].round(4).to_string(index=False))
