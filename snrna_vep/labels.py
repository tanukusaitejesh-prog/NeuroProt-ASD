"""Labels: curated patient variants, population (benign-proxy) variants, and SGE scores.

Patient table schema (data/curation/patient_variants.tsv, one row per variant per paper):
    gene_name        e.g. RNU4-2 (use GENCODE names: RNU2-2P for RNU2-2)
    hgvs_n           transcript-level HGVS on the GENCODE gene, e.g. n.64_65insT
    mode             AD-NDD | AR-NDD | AD-RP | AR-other
    inheritance      de_novo | inherited | biallelic_hom | biallelic_comphet | unknown
    n_probands       probands carrying it in that paper
    first_published  YYYY-MM of the first paper/preprint reporting it (drives the time split)
    source           DOI or PMID
    cohort           e.g. GEL-100kGP, DDD, SPARK (used to drop the same patient across papers)
    notes            free text
"""
import pandas as pd

from .coords import hgvs_n_to_vcf

PATIENT_COLUMNS = ["gene_name", "hgvs_n", "mode", "inheritance", "n_probands",
                   "first_published", "source", "cohort", "notes"]


def load_patients(path, genes, fetch):
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    missing = set(PATIENT_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"patient table missing columns: {sorted(missing)}")
    g = genes.set_index("gene_name")
    keys, errors = [], []
    for r in df.itertuples():
        try:
            gi = g.loc[r.gene_name]
            c, p, ref, alt = hgvs_n_to_vcf(r.hgvs_n.strip(), gi.chrom, int(gi.start), int(gi.end), gi.strand, fetch)
            keys.append(f"{c}:{p}:{ref}:{alt}")
        except Exception as e:  # keep going; report every bad row at once
            keys.append(None)
            errors.append(f"{r.gene_name} {r.hgvs_n}: {e}")
    df["key"] = keys
    df["n_probands"] = pd.to_numeric(df.n_probands, errors="coerce").fillna(1).astype(int)
    df["first_published"] = pd.to_datetime(df.first_published, format="%Y-%m", errors="coerce")
    return df, errors


def collapse_patients(df):
    """One row per unique variant x mode, keeping the earliest publication date."""
    agg = (df.dropna(subset=["key"])
             .groupby(["gene_name", "key", "mode"], as_index=False)
             .agg(hgvs_n=("hgvs_n", "first"), n_probands=("n_probands", "sum"),
                  first_published=("first_published", "min"), n_papers=("source", "nunique")))
    return agg


def population_controls(gnomad, task="dominant", min_ac=1, min_hom=1, patient_keys=()):
    """Benign proxies from gnomAD v4.1 genomes (PASS, not segdup).

    dominant task: any observed allele (severe dominant NDD is ~absent from gnomAD adults).
    recessive task: homozygotes observed (a carrier allele alone says nothing about AR disease).
    """
    x = gnomad[(gnomad["filter"] == "PASS") & ~gnomad.segdup.astype(bool)]
    x = x[(x.n_pos >= 1)]
    if task == "dominant":
        x = x[x.AC >= min_ac]
    elif task == "recessive":
        x = x[x.nhomalt >= min_hom]
    else:
        raise ValueError(task)
    x = x.assign(key=x.chrom + ":" + x.pos.astype(str) + ":" + x.ref + ":" + x.alt)
    return x[~x.key.isin(set(patient_keys))][["gene_name", "key", "AC", "AF", "nhomalt"]]


def load_sge(path, gene, fetch, score_col="score", hgvs_col=None):
    """SGE table with either an HGVS n. column or chrom/pos/ref/alt columns."""
    df = pd.read_csv(path, sep=None, engine="python")
    if hgvs_col:
        df["key"] = [":".join(map(str, hgvs_n_to_vcf(h, gene.chrom, gene.start, gene.end, gene.strand, fetch)))
                     for h in df[hgvs_col]]
    else:
        df["key"] = df.chrom.astype(str) + ":" + df.pos.astype(str) + ":" + df.ref + ":" + df.alt
    return df[["key", score_col]].rename(columns={score_col: "sge_score"})
