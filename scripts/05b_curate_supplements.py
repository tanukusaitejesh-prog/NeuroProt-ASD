"""Curate patient variants, independent controls and SGE scores from published supplements.

Supplements are downloaded from the article pages (Nature / Nature Genetics ESM) into SUPP (see
download notes in RESULTS_TONIGHT.md). Every row keeps its source table so each label is auditable.

Rules (fixed before looking at any model output):
  RNU4-2 SGE (De Jonghe et al., Nature 2026; preprint 2025-04): category ReNU_syndrome -> AD-NDD;
      UKBB/AllofUs -> population control (independent of gnomAD); function_score_mean -> SGE table.
  Chen et al. 2024 (Nature; preprint 2024-04) ST2: all individuals -> AD-NDD (de novo ReNU).
  Nava et al. 2025 (Nat Genet; preprint 2024-10) ST1 (RNU4-2) and ST6 (RNU5A-1/RNU5B-1): P/LP -> AD-NDD;
      ST2 'Patients LP/P' -> AD-NDD.
  Greene et al. 2025 (Nat Genet; preprint 2024-09) Feuil1: RNU2-2 n.4G>A, n.35A>G -> AD-NDD.
  Jackson et al. 2025 (Nat Genet 2025-05; date of publication used, conservative) ST5 (RNU2-2), ST7 (RNU5B-1):
      de novo -> AD-NDD.
  Systematic RNU2-2 (Nat Genet 2026; preprint 2025-09) ST8: P/LP; dominant if de novo / PS2, recessive if PM3 /
      segregation; ST6 inheritance column breaks ties.
  Biallelic RNU2-2 (Nat Genet 2026; preprint 2025-08) ST2: 'Dominant' -> AD-NDD; 'Tier 1' / 'Stronger evidence'
      -> AR-NDD; 'Weaker evidence' excluded.
  Biallelic RNU4-2 (Nat Genet 2026; preprint 2025-08) ST7: P/LP -> AR-NDD; ST2 (UK Biobank biallelic carriers)
      -> recessive-task controls.
  Retinitis pigmentosa (Nat Genet 2026-01; preprint 2025-01) S1, S3: P/LP -> AD-RP; S3 Benign -> benign control.
Variant coordinates are GRCh38, left-normalized to the pipeline's key (chrN:pos:ref:alt).
"""
import re
import sys
import warnings
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.coords import left_normalize
from snrna_vep.reference import fetch_base

warnings.filterwarnings("ignore")
SUPP = Path(sys.argv[1]) if len(sys.argv) > 1 else config.RAW / "supplements"
COORD = re.compile(r"(?:chr)?(\d{1,2}|X)\s*[:\-_ ]\s*(?:g\.)?(\d{6,9})\s*[:\-_ ]?\s*([ACGT]+)\s*(?:>|-|:|/|_)\s*([ACGT]+)", re.I)
G_INS = re.compile(r"(?:chr)?(\d{1,2})\s*:\s*g\.(\d+)_(\d+)ins([ACGT]+)", re.I)

genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t")
GI = [(r.chrom.replace("chr", ""), r.start - 5, r.end + 5, r.gene_name) for r in genes.itertuples()]


def gene_at(c, p):
    for gc, s, e, n in GI:
        if gc == c and s <= p <= e:
            return n


def keys_in(text):
    """All GRCh38 variants in a string, as (gene, normalized key)."""
    out = []
    for m in COORD.finditer(str(text)):
        c, p, ref, alt = m.group(1), int(m.group(2)), m.group(3).upper(), m.group(4).upper()
        out.append((c, p, ref, alt))
    for m in G_INS.finditer(str(text)):   # chr15:g.65304718_65304719insA
        c, a, ins = m.group(1), int(m.group(2)), m.group(4).upper()
        base = fetch_base("chr" + c, a)
        out.append((c, a, base, base + ins))
    res = []
    for c, p, ref, alt in out:
        g = gene_at(c, p)
        if not g:
            continue
        chrom = "chr" + c
        try:
            if fetch_base(chrom, p) != ref[0]:
                continue                      # reference mismatch: not GRCh38 or mis-parsed
            p2, r2, a2 = left_normalize(chrom, p, ref, alt, fetch_base) if len(ref) != len(alt) else (p, ref, alt)
        except Exception:
            continue
        res.append((g, f"{chrom}:{p2}:{r2}:{a2}"))
    return res


def xl(name, sheet, header):
    f = next(SUPP.glob(name))
    return pd.read_excel(f, sheet_name=sheet, header=header, dtype=str).fillna("")


P, C = [], []   # patients, controls


def add(rows, mode, inheritance, date, source, notes=""):
    for g, k in rows:
        P.append(dict(gene_name=g, key=k, hgvs_n="", mode=mode, inheritance=inheritance, n_probands=1,
                      first_published=date, source=source, cohort=source.split(" ")[0], notes=notes))


# --- RNU4-2 SGE table
sge = xl("41586_2026_10334_MOESM2*", "ST1_20251212_RNU4_2_merged_gf", 1)
sge_rows = []
for r in sge.itertuples():
    ks = keys_in(r.ID)
    if not ks:
        continue
    g, k = ks[0]
    if r.category == "ReNU_syndrome":
        add([(g, k)], "AD-NDD", "de_novo", "2024-04", "SGE2026 ST1", "ReNU category")
    elif r.category == "UKBB/AllofUs":
        C.append(dict(gene_name=g, key=k, control_type="UKB/AoU", task="dominant",
                      AC=pd.to_numeric(r.AoU_AC, errors="coerce"), hom=pd.to_numeric(r.UKBiobank_hom, errors="coerce"),
                      source="SGE2026 ST1"))
    s = pd.to_numeric(r.function_score_mean, errors="coerce")
    if pd.notna(s):
        sge_rows.append(dict(key=k, sge_score=s, sge_type=r.Type))

# --- Chen 2024 ST2
chen = xl("41586_2024_7773_MOESM3*", "Supplementary Table 2", 2)
for v in chen["Variant (GRCh38)"]:
    add(keys_in(v), "AD-NDD", "de_novo", "2024-04", "Chen2024 ST2")

# --- Nava 2025
for sh, lastcol in [("Supp. Table 1", -1), ("Supp. Table 6", -1)]:
    d = xl("41588_2025_2184_MOESM4*", sh, 2)
    for _, r in d.iterrows():
        cls = str(r.iloc[lastcol]).strip()
        if cls in ("P", "LP"):
            add(keys_in(" ".join(r.astype(str))), "AD-NDD", "de_novo", "2024-10", f"Nava2025 {sh}", cls)
d = xl("41588_2025_2184_MOESM4*", "Supp. Table 2", 2)
for _, r in d[d.group1 == "Patients LP/P"].iterrows():
    add(keys_in(r["Variant g."]), "AD-NDD", "de_novo", "2024-10", "Nava2025 Supp. Table 2", "Patients LP/P")

# --- Greene 2025 RNU2-2 (dominant)
d = xl("41588_2025_2159_MOESM4*", "Feuil1", None)
add(list(set(keys_in(" ".join(d.iloc[1].astype(str))))), "AD-NDD", "de_novo", "2024-09", "Greene2025 Feuil1")

# --- Jackson 2025 R-loop
for sh in ["Supplementary Table 5", "Supplementary Table 7"]:
    d = xl("41588_2025_2209_MOESM3*", sh, 2)
    inh = [c for c in d.columns if "nherit" in c][0]
    for _, r in d.iterrows():
        if "novo" in str(r[inh]).lower():
            add(keys_in(" ".join(r.astype(str))), "AD-NDD", "de_novo", "2025-05", f"Jackson2025 {sh}")

# --- Systematic RNU2-2 2026: ST8 classifications, ST6 inheritance
st6 = xl("41588_2026_2547_MOESM4*", "Supp. Table 6", 2)
mode6 = {}
for _, r in st6.iterrows():
    txt = " ".join(r.astype(str))
    m = "AD-NDD" if "Dominant" in txt else ("AR-NDD" if "Recessive" in txt else None)
    for g, k in keys_in(txt):
        if m:
            mode6.setdefault(k, set()).add(m)
st8 = xl("41588_2026_2547_MOESM4*", "Supp. Table 8", 2)
for _, r in st8.iterrows():
    if r["Classification"].strip() not in ("P", "LP"):
        continue
    crit = r["ACMG criteria"]
    for g, k in keys_in(r["variant g."]):
        if r["de novo (our study)"] == "dn" or "PS2" in crit:
            m = "AD-NDD"
        elif "PM3" in crit or r["segregation (our study)"]:
            m = "AR-NDD"
        elif len(mode6.get(k, ())) == 1:
            m = next(iter(mode6[k]))
        else:
            continue
        add([(g, k)], m, "de_novo" if m == "AD-NDD" else "biallelic", "2025-09", "Systematic2026 ST8", r["Classification"])

# --- Biallelic RNU2-2 2026 ST2
d = xl("41588_2026_2539_MOESM3*", "ST2", 0)
for _, r in d.iterrows():
    grp = r["Fig. 5 group"]
    m = {"Dominant": "AD-NDD", "Tier 1": "AR-NDD", "Stronger evidence": "AR-NDD"}.get(grp)
    if m:
        add(keys_in(r["Variant (GRCh38)"]), m, "de_novo" if m == "AD-NDD" else "biallelic", "2025-08",
            "BiallelicRNU2-2 2026 ST2", grp)

# --- Biallelic RNU4-2 2026
d = xl("41588_2026_2554_MOESM3*", "Sup. Table 7", 2)
for _, r in d.iterrows():
    if r["Classification"].strip().lower() in ("pathogenic", "likely pathogenic"):
        add(keys_in(r["variant (GRCh38)"]), "AR-NDD", "biallelic", "2025-08", "BiallelicRNU4-2 2026 ST7",
            r["Classification"])
d = xl("41588_2026_2554_MOESM3*", "Sup. Table 2", 2)
for _, r in d.iterrows():
    for g, k in keys_in(r["variant (GRCh38)"]):
        C.append(dict(gene_name=g, key=k, control_type="UKB biallelic", task="recessive", AC=None, hom=None,
                      source="BiallelicRNU4-2 2026 ST2"))

# --- Retinitis pigmentosa
for sh in ["S1 Recurrent variants", "S3 Variants identified Sanger"]:
    d = xl("41588_2025_2451_MOESM3*", sh, 2)
    for _, r in d.iterrows():
        cls = str(r.get("ACMG classification", "")).strip().lower()
        try:
            v = f"{int(float(r['Chr (hg38)']))}:{int(float(r['Position (hg38)']))} {r['Reference (hg38)']}>{r['Alternative (hg38)']}"
        except (ValueError, TypeError):
            continue
        if cls in ("pathogenic", "likely pathogenic"):
            add(keys_in(v), "AD-RP", "dominant", "2025-01", f"RP2026 {sh}", cls)
        elif cls == "benign":
            for g, k in keys_in(v):
                C.append(dict(gene_name=g, key=k, control_type="RP benign", task="dominant", AC=None, hom=None,
                              source=f"RP2026 {sh}"))

pat = pd.DataFrame(P)
ctl = pd.DataFrame(C).drop_duplicates(["key", "task"])
pat = pat.drop_duplicates(["key", "mode", "source"])
config.CURATION.mkdir(parents=True, exist_ok=True)
pat.to_csv(config.CURATION / "patient_variants_curated.tsv", sep="\t", index=False)
ctl.to_csv(config.CURATION / "controls_curated.tsv", sep="\t", index=False)
(config.DATA / "external").mkdir(exist_ok=True)
pd.DataFrame(sge_rows).drop_duplicates("key").to_csv(config.DATA / "external" / "sge_rnu4-2.tsv", sep="\t", index=False)

u = pat.drop_duplicates(["gene_name", "key", "mode"])
print("patient rows", len(pat), "| unique variant x mode", len(u))
print(u.groupby(["gene_name", "mode"]).size().unstack(fill_value=0).to_string())
print("\ncontrols", ctl.groupby(["gene_name", "control_type"]).size().to_string())
print("\nSGE scores", len(sge_rows))
both = u.groupby("key")["mode"].nunique()
print("variants with >1 mode:", int((both > 1).sum()))
