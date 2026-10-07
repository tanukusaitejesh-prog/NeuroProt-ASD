"""Baseline / external scores: CADD, phyloP, AlphaGenome Atlas, Evo 2, RNA language models.

Each loader returns a DataFrame keyed on the variant `key` (chrom:pos:ref:alt, left-normalized).
Remote files are read by range requests, so nothing large is downloaded. These hosts are blocked in
the cloud session; run locally.
"""
import numpy as np
import pandas as pd

CADD_SNV = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
CADD_INDEL = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz"
PHYLOP_BW = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP447way/hg38.phyloP447way.bw"


def cadd_scores(variants, url=CADD_SNV):
    """PHRED for every SNV key. Indels need CADD's web scorer (not precomputed for unseen indels)."""
    import pysam
    tb = pysam.TabixFile(url)
    want = set(variants.key)
    out = []
    for (chrom, lo, hi) in variants.groupby("chrom").pos.agg(["min", "max"]).itertuples():
        for rec in tb.fetch(chrom.replace("chr", ""), lo - 1, hi):
            c, p, r, a, raw, phred = rec.split("\t")[:6]
            key = f"chr{c}:{p}:{r}:{a}"
            if key in want:
                out.append((key, float(phred)))
    return pd.DataFrame(out, columns=["key", "cadd_phred"])


def phylop_scores(variants, url=PHYLOP_BW):
    """Per-base phyloP at the variant position (for indels, the anchor+1 base)."""
    import pyBigWig
    bw = pyBigWig.open(url)
    vals = []
    for v in variants.itertuples():
        p = v.pos if len(v.ref) == len(v.alt) else v.pos + 1
        x = bw.values(v.chrom, p - 1, p)[0]
        vals.append(np.nan if x is None else x)
    return pd.DataFrame({"key": variants.key.values, "phylop": vals})


def alphagenome_atlas(path):
    """Load a locally downloaded AlphaGenome Atlas slice (columns: chrom,pos,ref,alt,<score...>).

    The Atlas file layout is not known from the cloud session; adjust the column names here once.
    """
    df = pd.read_csv(path, sep=None, engine="python")
    df["key"] = df.chrom.astype(str).str.replace("^(?!chr)", "chr", regex=True) + ":" + \
        df.pos.astype(str) + ":" + df.ref + ":" + df.alt
    score_col = [c for c in df.columns if c.lower() in ("avi", "avi_score", "alphagenome_avi", "score")][0]
    return df[["key", score_col]].rename(columns={score_col: "alphagenome_avi"})


def rna_lm_delta_ll(variants, rna_by_gene, model_name="multimolecule/rnafm"):
    """Masked-LM log-likelihood ratio (mutant vs WT) from an RNA language model (needs HF access)."""
    import torch
    from multimolecule import RnaTokenizer, RnaFmForMaskedLM
    from .variants import apply_variant
    tok = RnaTokenizer.from_pretrained(model_name)
    lm = RnaFmForMaskedLM.from_pretrained(model_name).eval()

    def pll(seq):
        ids = tok(seq.replace("T", "U"), return_tensors="pt")["input_ids"]
        with torch.no_grad():
            logp = torch.log_softmax(lm(input_ids=ids).logits, -1)
        return logp[0, torch.arange(1, ids.shape[1] - 1), ids[0, 1:-1]].sum().item()

    cache, rows = {}, []
    for v in variants.itertuples():
        wt = rna_by_gene[v.gene_name]
        if wt not in cache:
            cache[wt] = pll(wt)
        rows.append((v.key, pll(apply_variant(wt, v.hgvs)) - cache[wt]))
    return pd.DataFrame(rows, columns=["key", "rnalm_dll"])


def evo2_delta_ll(variants, fetch_region, model_name="evo2_7b", context=4096):
    """Evo 2 log-likelihood ratio on genomic context around each variant (needs a GPU)."""
    from evo2 import Evo2
    model = Evo2(model_name)
    rows = []
    for v in variants.itertuples():
        lo = v.pos - context // 2
        wt = fetch_region(v.chrom, lo, lo + context - 1)
        i = v.pos - lo
        mut = wt[:i] + v.alt + wt[i + len(v.ref):]
        s = model.score_sequences([wt, mut])
        rows.append((v.key, float(s[1] - s[0])))
    return pd.DataFrame(rows, columns=["key", "evo2_dll"])
