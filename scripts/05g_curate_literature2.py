"""Curate 2026 minor-spliceosome papers (RNU4ATAC, RNU6ATAC) from GRCh38 genomic coordinates given in the papers.

Sources (supplements/tables downloaded beforehand into <dir>; Elsevier mmc files, Europe PMC full text):
  GIM2026 RNU4ATAC cohort  (Matalon et al., Genet Med 2026, doi 10.1016/j.gim.2026.102633), Supplemental Table 2
      'All variants': column 'Genomic coordinates (NC_000002.12) [GRCh38]'; Final ACMG classification P/LP kept, VUS dropped.
  AJHG2026 RNU6ATAC/RNU4ATAC (doi 10.1016/j.ajhg.2026.02.017), Tables 1-2: hg38 'Genomic co-ordinate' rows of
      affected individuals with bi-allelic 'disease-causing variants'.
  HGGA2026 RNU6ATAC (Bi-allelic RNU6ATAC variants ..., HGG Adv 2026, PMC13049632), Table 1: hg38 '9-pos-ref-alt'.
All three report bi-allelic (recessive) disease; mode 'AR-MS' (recessive, multisystem / NDD).
  iSci2026 compiled (Whole-genome discovery of pathogenic snRNA variants..., iScience 2026, doi
      10.1016/j.isci.2026.116814) Table S1: variants reported in patients by 65 publications (gene, n. variant, AD/AR,
      source). No per-variant classification, so it is kept as a separate tier ('compiled'); variants explicitly
      classified VUS by GIM2026 are excluded; n_sources = number of distinct first-author publications.
Every allele is checked against the GRCh38 reference base. Output: data/curation/literature2_curated.tsv
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.reference import fetch_base, fetch_region

root = Path(sys.argv[1])
rows = []


def add(chrom, pos, ref, alt, gene, source, date, cls="", zyg=""):
    if fetch_base(chrom, int(pos)) != ref[0]:
        print("  reference mismatch, skipped:", chrom, pos, ref, alt, source)
        return
    rows.append(dict(gene_name=gene, key=f"{chrom}:{int(pos)}:{ref}:{alt}", mode="AR-MS", inheritance="biallelic",
                     first_published=date, source=source, notes=f"{cls} {zyg}".strip()))


def g_del(chrom, a, b):
    """HGVS g.a_b del -> left-anchored VCF (ref = base before + deleted bases)."""
    seq = fetch_region(chrom, int(a) - 1, int(b))
    return int(a) - 1, seq, seq[0]


# ---- GIM2026 Supplemental Table 2
d = pd.read_excel(root / "els/S1098360026009512_mmc3.xlsx", sheet_name="All variants", header=1)
gcol = [c for c in d.columns if "Genomic coordinates" in c][0]
ccol = [c for c in d.columns if "Final ACMG" in c][0]
for g, cls, z in zip(d[gcol].astype(str), d[ccol].astype(str), d["Zygosity"].astype(str)):
    if not re.search("pathogenic", cls, re.I):
        continue
    m = re.match(r"g\.(\d+)([ACGT])>([ACGT])", g.strip())
    if m:
        add("chr2", m.group(1), m.group(2), m.group(3), "RNU4ATAC", "GIM2026 RNU4ATAC ST2", "2026-06", cls, z)
        continue
    m = re.match(r"g\.(\d+)_(\d+)del", g.strip())
    if m:
        p, r, a = g_del("chr2", m.group(1), m.group(2))
        add("chr2", p, r, a, "RNU4ATAC", "GIM2026 RNU4ATAC ST2", "2026-06", cls, z)
    else:
        print("  unparsed GIM2026:", g, cls)

# ---- AJHG2026 Tables 1-2 (main-text PDF, layout text): every hg38 SNV coordinate in the two tables
t = (root / "els/b.txt").read_text()
lines = t.splitlines()
s = next(i for i, l in enumerate(lines) if "Table 1." in l and "RNU6ATAC" in l)
e = next(i for i, l in enumerate(lines) if "Furthermore, all individuals" in l)
for pos, ref, alt in set(re.findall(r"(134164\d{3}|1215\d{5})\s?([ACGT])>([ACGT])", "\n".join(lines[s:e]))):
    chrom, gene = ("chr9", "RNU6ATAC") if pos.startswith("134164") else ("chr2", "RNU4ATAC")
    add(chrom, pos, ref, alt, gene, "AJHG2026 RNU6ATAC/RNU4ATAC Tables 1-2", "2026-03")

# ---- HGGA2026 Table 1 ('9-134164537-G-A')
x = (root / "PMC13049632.xml").read_text(encoding="utf8", errors="ignore")
for pos, ref, alt in set(re.findall(r"\b9-(1341645\d\d)-([ACGT])-([ACGT])\b", x)):
    add("chr9", pos, ref, alt, "RNU6ATAC", "HGGA2026 RNU6ATAC Table 1", "2026-03")

# ---- iSci2026 Table S1 (compiled literature)
from snrna_vep.literature import to_key
gim_vus = set()
for g, cls in zip(d[gcol].astype(str), d[ccol].astype(str)):
    m = re.match(r"g\.(\d+)([ACGT])>([ACGT])", g.strip())
    if m and "VUS" in cls:
        gim_vus.add(f"chr2:{m.group(1)}:{m.group(2)}:{m.group(3)}")
c = pd.read_excel(root / "els/S2589004226021929_mmc2.xlsx", header=2).dropna(subset=["Gene", "Variant"])
c["first_author"] = c.Source.astype(str).str.split(",").str[0].str.strip()
gene_map = {"RNU2-2": "RNU2-2P"}
for (gene, var, moi), grp in c.groupby(["Gene", "Variant", "Mode of inheritance"]):
    gene = gene_map.get(gene, gene)
    m = re.match(r"n\.(\d+(?:_\d+)?)([ACGTU]>[ACGTU]|del[ACGTU]*|ins[ACGTU]+|dup[ACGTU]*)$", str(var).strip())
    key = to_key(gene, m.group(1), m.group(2)) if m else None
    if key is None:
        print("  iSci2026 unconverted:", gene, var)
        continue
    if key in gim_vus:
        continue
    rows.append(dict(gene_name=gene, key=key, mode="AD-NDD" if moi == "AD" else "AR-NDD",
                     inheritance="de_novo" if moi == "AD" else "biallelic", first_published="2026-08",
                     source="iSci2026 compiled", notes=f"n_sources={grp.first_author.nunique()}"))

out = pd.DataFrame(rows)
out.to_csv(config.CURATION / "literature2_curated.tsv", sep="\t", index=False)
u = out.drop_duplicates(["key"])
print(f"{len(out)} alleles, {len(u)} distinct variants")
print(u.groupby(["gene_name"]).size().to_string())
print(out.drop_duplicates(["key", "source"]).groupby(["source"]).size().to_string())
