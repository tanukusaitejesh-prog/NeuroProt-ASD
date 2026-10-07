"""Labels v2: curated supplements + latest ClinVar, with mode-aware controls.

Pathogenic (label 1):
  - curated supplement variants (scripts/05b_curate_supplements.py), any mode;
  - ClinVar 2026-10-04 P/LP (any review status >= 1 star) in spliceosomal snRNA genes. Mode from gene and the
    ClinVar condition text: 'retinitis' -> AD-RP; 'recessive' or genes with only recessive disease (RNU4ATAC, RNU12)
    -> AR; otherwise AD-NDD.
Controls (label 0), mode-aware:
  - RNU4-2: UK Biobank / All of Us variants (SGE supplement) + ClinVar B/LB.
  - genes with dominant disease only (RNU5B-1, RNU5A-1, RNU6 RP genes): gnomAD v4.1 PASS AC >= 1 + B/LB.
  - genes with recessive disease (RNU2-2, RNU4ATAC, RNU12): gnomAD homozygotes (nhomalt >= 1) + ClinVar B/LB +
    UKB biallelic carriers; heterozygous population carriage is not evidence of benignity for recessive disease.
  - RP-benign variants (RP supplement) for the U4/U6 genes.
A variant labelled both ways is dropped from controls. Output: data/curation/labels_v2.tsv
--v3: additionally the 2026 minor-spliceosome papers (scripts/05g_curate_literature2.py: RNU4ATAC cohort, RNU6ATAC);
RNU6ATAC (recessive disease) gets recessive-mode controls. Output: data/curation/labels_v3.tsv
--v4: v3 + the compiled literature tier (iScience 2026 Table S1).            Output: labels_v4.tsv
--v4s: as v4, but compiled variants need >= 2 independent publications.        Output: labels_v4s.tsv
(v3 labels are unchanged by the compiled tier: it is excluded unless --v4/--v4s is given.)
--v5: v4 + UK Biobank WGS homozygote controls for recessive genes (scripts/05h). Output: labels_v5.tsv
--v5c: as v5, but variants both labelled pathogenic and homozygous in UKB are removed entirely. Output: labels_v5c.tsv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config

SPLICEOSOMAL = {"RNU4-2", "RNU4-1", "RNU2-2P", "RNU5B-1", "RNU5A-1", "RNU6-1", "RNU6-2", "RNU6-7", "RNU6-8", "RNU6-9",
                "RNU4ATAC", "RNU6ATAC", "RNU12", "RNU11", "RNU5D-1", "RNU5E-1", "RNU5F-1"}
V5C = "--v5c" in sys.argv
V5 = "--v5" in sys.argv or V5C
V4S = "--v4s" in sys.argv
V4 = "--v4" in sys.argv or V4S or V5
V3 = "--v3" in sys.argv or V4
RECESSIVE_GENES = {"RNU2-2P", "RNU4ATAC", "RNU12"} | ({"RNU6ATAC"} if V3 else set())
DOMINANT_ONLY = {"RNU5B-1", "RNU5A-1", "RNU6-1", "RNU6-2", "RNU6-8", "RNU6-9"}

pat = pd.read_csv(config.CURATION / "patient_variants_curated.tsv", sep="\t")
ctl = pd.read_csv(config.CURATION / "controls_curated.tsv", sep="\t")
cv = pd.read_csv(config.PROCESSED / "clinvar_latest_core_loci.tsv", sep="\t")
gn = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")

cv = cv[cv.gene_name.isin(SPLICEOSOMAL) & ~cv.revstat.fillna("").str.contains("no_assertion|no_classification")]
cv["key"] = cv.chrom + ":" + cv.pos.astype(str) + ":" + cv.ref + ":" + cv.alt
plp = cv[cv.clnsig.isin(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"])].copy()
blb = cv[cv.clnsig.isin(["Benign", "Likely_benign", "Benign/Likely_benign"])].copy()


def cv_mode(r):
    t = str(r.disease).lower()
    if "retinitis" in t:
        return "AD-RP"
    if "recessive" in t or r.gene_name in ("RNU4ATAC", "RNU12"):
        return "AR-NDD"
    return "AD-NDD"


plp["mode"] = plp.apply(cv_mode, axis=1)
P = pd.concat([pat[["gene_name", "key", "mode", "source"]],
               plp.assign(source="ClinVar 2026-10")[["gene_name", "key", "mode", "source"]]])
if V3:
    lit = pd.read_csv(config.CURATION / "literature2_curated.tsv", sep="\t")
    comp = lit.source == "iSci2026 compiled"
    nsrc = lit.notes.astype(str).str.extract(r"n_sources=(\d+)")[0].astype(float)
    keep_comp = (nsrc >= 2) if V4S else pd.Series(V4, index=lit.index)
    lit = lit[~comp | keep_comp]
    P = pd.concat([P, lit[["gene_name", "key", "mode", "source"]]])
P = P.drop_duplicates(["key", "mode"])
pkeys = set(P.key)

LEN = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t").set_index("gene_name").length
gp = gn[(gn["filter"] == "PASS") & ~gn.segdup.astype(bool)].copy()
gp = gp[(gp.n_pos >= 0) & (gp.n_pos <= gp.gene_name.map(LEN) + 1)]          # transcribed region only (indel anchors)
gp["key"] = gp.chrom + ":" + gp.pos.astype(str) + ":" + gp.ref + ":" + gp.alt
C = []
for g in sorted(set(gp.gene_name) | set(ctl.gene_name)):
    if g == "RNU4-2":
        ks = set(ctl[(ctl.gene_name == g) & (ctl.control_type == "UKB/AoU")].key)
        ks |= set(ctl[(ctl.gene_name == g) & (ctl.control_type == "UKB biallelic")].key)
    elif g in RECESSIVE_GENES:
        ks = set(gp[(gp.gene_name == g) & (gp.nhomalt >= 1)].key)
    else:
        ks = set(gp[(gp.gene_name == g) & (gp.AC >= 1)].key)
    if V5 and g in RECESSIVE_GENES and (config.CURATION / "ukb_afb_homozygotes.tsv").exists():
        ukb = pd.read_csv(config.CURATION / "ukb_afb_homozygotes.tsv", sep="\t")
        ks |= set(ukb[ukb.gene_name == g].key)
    ks |= set(blb[blb.gene_name == g].key)
    ks |= set(ctl[(ctl.gene_name == g) & (ctl.control_type == "RP benign")].key)
    for k in ks - pkeys:
        C.append(dict(gene_name=g, key=k, mode="", source="control"))

out = pd.concat([P.assign(label=1), pd.DataFrame(C).assign(label=0)], ignore_index=True)
if V5C:
    ukb = pd.read_csv(config.CURATION / "ukb_afb_homozygotes.tsv", sep="\t")
    out = out[~out.key.isin(set(ukb[ukb.conflict].key))]
out.to_csv(config.CURATION / ("labels_v5c.tsv" if V5C else "labels_v5.tsv" if V5 else "labels_v4s.tsv" if V4S else "labels_v4.tsv" if V4 else "labels_v3.tsv" if V3 else "labels_v2.tsv"), sep="\t", index=False)
t = out.groupby(["gene_name", "label"]).size().unstack(fill_value=0)
t.columns = ["controls", "pathogenic"]
print(t[t.pathogenic > 0].to_string())
print(out[out.label == 1].groupby(["gene_name", "mode"]).size().unstack(fill_value=0).to_string())
