"""Day-one confound test for the tissue-specificity question.

Hypothesis: tissues that fail first when the spliceosome is perturbed (brain, retina) are enriched for transcripts
that depend on WEAK 5' splice sites, so a global perturbation hits them first. For the retina this is the
demonstrated mechanism (Nat Commun 2024, PRPF8/Brr2: weak/suboptimal 5'SS selection fails in retina-specific
transcripts). For the brain it has never been tested.

THE TEST THAT COULD KILL IT: brain-expressed genes are longer, have more introns, and are more highly expressed
than average. Any raw "brain has weaker 5' splice sites" signal may be those covariates. This script asks whether
the difference survives matching on intron count, gene length (log) and expression level.

Inputs (public): GENCODE v44 annotation, GRCh38 primary assembly, GTEx v8 median TPM by tissue.
5'SS score: position weight matrix over the 9-mer (3 exonic + 6 intronic) built from ALL annotated internal donor
sites, scored as a log-odds against position-specific background. Higher = stronger / closer to consensus.

Outputs: results/splice_strength_by_tissue.tsv, results/splice_strength_matched.tsv
"""
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

TISS = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("tiss")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


# ---------- 1. donor sites from GENCODE
print("parsing GENCODE ...", flush=True)
exons = defaultdict(list)          # transcript -> [(chrom, start, end, strand, gene, gene_name)]
with gzip.open(TISS / "gencode.gtf.gz", "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "exon":
            continue
        a = f[8]
        if 'transcript_type "protein_coding"' not in a:
            continue
        tid = a.split('transcript_id "')[1].split('"')[0]
        gid = a.split('gene_id "')[1].split('"')[0].split(".")[0]
        gname = a.split('gene_name "')[1].split('"')[0]
        exons[tid].append((f[0], int(f[3]), int(f[4]), f[6], gid, gname))
print(f"  protein-coding transcripts: {len(exons)}", flush=True)

# donor site = last base of an exon that is followed by an intron (internal/5' exons)
donors = []                        # (chrom, pos0, strand, gene, gene_name, transcript)
tx_meta = {}
for tid, ex in exons.items():
    ex.sort(key=lambda e: e[1])
    chrom, strand, gid, gname = ex[0][0], ex[0][3], ex[0][4], ex[0][5]
    span = ex[-1][2] - ex[0][1]
    tx_meta[tid] = dict(gene=gid, gene_name=gname, n_exons=len(ex), tx_len=span)
    if len(ex) < 2:
        continue
    if strand == "+":
        for c, s, e, st, g, gn in ex[:-1]:
            donors.append((chrom, e, "+", gid, gn, tid))
    else:
        for c, s, e, st, g, gn in ex[1:]:
            donors.append((chrom, s, "-", gid, gn, tid))
D = pd.DataFrame(donors, columns=["chrom", "pos", "strand", "gene", "gene_name", "tx"]).drop_duplicates(
    ["chrom", "pos", "strand"]).reset_index(drop=True)
print(f"  unique donor sites: {len(D)}", flush=True)

# ---------- 2. pull the 9-mers by streaming the genome once
need = defaultdict(list)
for i, (c, p, s) in enumerate(zip(D.chrom, D.pos, D.strand)):
    need[c].append((i, p, s))
seqs = [None] * len(D)
print("streaming genome ...", flush=True)
cur, buf = None, []


def flush(chrom, seq):
    if chrom not in need:
        return
    for i, p, s in need[chrom]:
        if s == "+":                       # exon ...| intron -> 3 exonic + 6 intronic around p (1-based last exon base)
            w = seq[p - 3:p + 6]
        else:                              # minus strand: p is the exon's leftmost base
            w = rc(seq[p - 7:p + 2])
        if len(w) == 9 and "N" not in w.upper():
            seqs[i] = w.upper()


with gzip.open(TISS / "genome.fa.gz", "rt") as fh:
    for line in fh:
        if line[0] == ">":
            if cur is not None:
                flush(cur, "".join(buf))
            cur, buf = line[1:].split()[0], []
            print(f"  {cur}", end=" ", flush=True)
        else:
            buf.append(line.strip())
    if cur is not None:
        flush(cur, "".join(buf))
print()
D["site"] = seqs
D = D[D.site.notna()].reset_index(drop=True)
print(f"  donor 9-mers extracted: {len(D)}")
print("  consensus check (positions 4-5 should be GT):", pd.Series([s[3:5] for s in D.site]).value_counts().head(3).to_dict())

# ---------- 3. PWM log-odds score
M = np.zeros((9, 4))
IDX = {c: i for i, c in enumerate("ACGT")}
for s in D.site:
    for j, ch in enumerate(s):
        M[j, IDX[ch]] += 1
M = (M + 1) / (M + 1).sum(1, keepdims=True)
bg = np.full(4, 0.25)
LO = np.log2(M / bg)
D["ss5"] = [sum(LO[j, IDX[ch]] for j, ch in enumerate(s)) for s in D.site]
print(f"  5'SS score: mean {D.ss5.mean():.2f}, sd {D.ss5.std():.2f}")

# ---------- 4. per-gene summaries
TX = pd.DataFrame(tx_meta).T.reset_index(names="tx")
g = D.groupby("gene")
G = pd.DataFrame({"n_donors": g.size(), "ss5_mean": g.ss5.mean(), "ss5_min": g.ss5.min(),
                  "ss5_q10": g.ss5.quantile(0.1), "frac_weak": g.ss5.apply(lambda v: (v < D.ss5.quantile(0.1)).mean()),
                  "gene_name": g.gene_name.first()}).reset_index()
TX["n_exons"] = TX.n_exons.astype(int)
TX["tx_len"] = TX.tx_len.astype(int)
tx_by_gene = TX.groupby("gene").agg(n_exons=("n_exons", "max"), tx_len=("tx_len", "max")).reset_index()
G = G.merge(tx_by_gene, on="gene", how="left")

# ---------- 5. GTEx tissue specificity
print("reading GTEx ...", flush=True)
X = pd.read_csv(TISS / "gtex_med.gct.gz", sep="\t", skiprows=2)
X["gene"] = X.Name.str.split(".").str[0]
tis = [c for c in X.columns if c not in ("Name", "Description", "gene")]
brain_cols = [c for c in tis if c.startswith("Brain")]
other_cols = [c for c in tis if not c.startswith("Brain")]
X = X.set_index("gene")
bm, om = X[brain_cols].max(1), X[other_cols].max(1)
E = pd.DataFrame({"brain_max": bm, "other_max": om, "expr": X[tis].max(1)}).reset_index()
E["brain_spec"] = np.log2((E.brain_max + 1) / (E.other_max + 1))
G = G.merge(E, on="gene", how="inner")
G = G[(G.expr > 1) & (G.n_donors >= 2)].copy()
for c in ["n_exons", "tx_len", "expr", "ss5_mean", "ss5_min", "frac_weak"]:
    G[c] = pd.to_numeric(G[c], errors="coerce")
G = G.dropna(subset=["n_exons", "tx_len", "expr"])
print(f"  genes with donors + GTEx: {len(G)}")

G["class"] = np.where(G.brain_spec >= 1, "brain-enriched", np.where(G.brain_spec <= -1, "other-enriched", "ubiquitous"))
print("\n=== RAW: 5'SS strength by tissue class ===")
print(G.groupby("class")[["ss5_mean", "ss5_min", "frac_weak", "n_exons", "tx_len", "expr"]].median().round(3).to_string())
print(G["class"].value_counts().to_string())

# ---------- 6. the confound test: match on intron count, length, expression
from scipy.stats import mannwhitneyu
B = G[G["class"] == "brain-enriched"].copy()
O = G[G["class"] == "other-enriched"].copy()
print(f"\n=== MATCHED: each brain-enriched gene paired to an other-enriched gene with similar n_exons, log length, log expression ===")
O = O.reset_index(drop=True)
key = np.column_stack([np.log2(O.n_exons.clip(1)), np.log2(O.tx_len.clip(1)), np.log2(O.expr.clip(0.1))])
pairs = []
used = set()
for r in B.itertuples():
    q = np.array([np.log2(max(r.n_exons, 1)), np.log2(max(r.tx_len, 1)), np.log2(max(r.expr, 0.1))])
    d = np.abs(key - q).sum(1)
    d[list(used)] = np.inf
    j = int(np.argmin(d))
    if np.isfinite(d[j]) and d[j] < 1.5:
        used.add(j)
        pairs.append((r.ss5_mean, O.ss5_mean.iloc[j], r.ss5_min, O.ss5_min.iloc[j], r.frac_weak, O.frac_weak.iloc[j]))
P = pd.DataFrame(pairs, columns=["b_mean", "o_mean", "b_min", "o_min", "b_weak", "o_weak"])
print(f"  matched pairs: {len(P)}")
rows = []
for lab, bcol, ocol in [("ss5_mean", "b_mean", "o_mean"), ("ss5_min", "b_min", "o_min"), ("frac_weak", "b_weak", "o_weak")]:
    diff = P[bcol] - P[ocol]
    from scipy.stats import wilcoxon
    try:
        st, p = wilcoxon(P[bcol], P[ocol])
    except Exception:
        st, p = np.nan, np.nan
    rows.append(dict(metric=lab, brain_median=P[bcol].median(), other_median=P[ocol].median(),
                     median_diff=diff.median(), wilcoxon_p=p))
R = pd.DataFrame(rows)
print(R.round(4).to_string(index=False))
verdict = "SURVIVES matching" if (R.loc[R.metric == "ss5_mean", "wilcoxon_p"].iloc[0] < 0.05 and
                                  R.loc[R.metric == "ss5_mean", "median_diff"].iloc[0] < 0) else "DOES NOT survive matching"
print(f"\nDAY-ONE VERDICT: brain-enriched genes have weaker 5' splice sites -> {verdict}")
G.to_csv(config.RESULTS / "splice_strength_by_tissue.tsv", sep="\t", index=False)
R.to_csv(config.RESULTS / "splice_strength_matched.tsv", sep="\t", index=False)
D[["chrom", "pos", "strand", "gene", "gene_name", "site", "ss5"]].to_csv(
    config.PROCESSED / "donor_sites_ss5.tsv.gz", sep="\t", index=False)
