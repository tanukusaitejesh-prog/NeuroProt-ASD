"""Small-RNA gene catalogue (GENCODE v27, GRCh38) and disease-gene annotations."""
import re

import pandas as pd

# Disease associations established in the literature. Each row cites the paper it comes from.
# "modes" uses: AD-NDD dominant neurodevelopmental, AR-NDD recessive neurodevelopmental,
# AD-RP dominant retinitis pigmentosa, AR-other recessive non-NDD-predominant disorder.
DISEASE_GENES = [
    ("RNU4-2", "AD-NDD;AR-NDD;AD-RP", 2024,
     "Chen 2024 Nature s41586-024-07773-7; SGE Nature 2026 s41586-026-10334-9; "
     "biallelic Nat Genet 2026 s41588-026-02554-6; RP Nat Genet 2025 s41588-025-02451-4"),
    ("RNU2-2P", "AD-NDD;AR-NDD", 2025,
     "Nat Genet 2025 s41588-025-02159-5; Nat Genet 2025 s41588-025-02209-y; "
     "Nat Genet 2026 s41588-026-02539-5, s41588-026-02551-9, s41588-026-02547-5"),
    ("RNU5B-1", "AD-NDD", 2025, "Nat Genet 2025 s41588-025-02209-y; Nat Genet 2025 s41588-025-02184-4"),
    ("RNU5A-1", "AD-NDD", 2025, "Nat Genet 2025 s41588-025-02184-4"),
    ("RNU4ATAC", "AR-NDD", 2011, "MOPD1/Roifman/Lowry-Wood syndromes"),
    ("RNU12", "AR-other", 2015, "early-onset cerebellar ataxia; CDAGS"),
    ("RNU7-1", "AR-other", 2020, "Aicardi-Goutieres syndrome"),
    ("RMRP", "AR-other", 2001, "cartilage-hair hypoplasia"),
]

# Functional (expressed) paralog groups used for paralog-aware transfer.
PARALOG_GROUPS = {
    "U4": ["RNU4-1", "RNU4-2"],
    "U4atac": ["RNU4ATAC"],
    "U2": ["RNU2-1", "RNU2-2P"],
    "U5": ["RNU5A-1", "RNU5B-1", "RNU5D-1", "RNU5E-1", "RNU5F-1"],
    "U6": ["RNU6-1", "RNU6-2", "RNU6-7", "RNU6-8", "RNU6-9"],
    "U6atac": ["RNU6ATAC"],
    "U1": ["RNU1-1", "RNU1-2", "RNU1-3", "RNU1-4"],
    "U11": ["RNU11"],
    "U12": ["RNU12"],
    "U7": ["RNU7-1"],
}

_ATTR = re.compile(r'(\S+) "([^"]*)"')


def parse_gtf(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(_ATTR.findall(f[8]))
            rows.append({
                "gene_id": attrs.get("gene_id"),
                "gene_name": attrs.get("gene_name"),
                "gene_type": attrs.get("gene_type"),
                "chrom": f[0],
                "start": int(f[3]),
                "end": int(f[4]),
                "strand": f[6],
            })
    df = pd.DataFrame(rows)
    df["length"] = df["end"] - df["start"] + 1
    return df


def build_catalogue(gtf_path):
    df = parse_gtf(gtf_path)
    # GENCODE carries duplicate RMRP records; keep the longest per name+chrom.
    df = (df.sort_values("length", ascending=False)
            .drop_duplicates(["gene_name", "chrom"]).reset_index(drop=True))
    dis = pd.DataFrame(DISEASE_GENES, columns=["gene_name", "disease_modes", "first_reported", "evidence"])
    df = df.merge(dis, on="gene_name", how="left")
    df["is_disease_gene"] = df["disease_modes"].notna()
    group_of = {g: grp for grp, genes in PARALOG_GROUPS.items() for g in genes}
    df["paralog_group"] = df["gene_name"].map(group_of)
    df["likely_pseudogene"] = df["gene_name"].str.contains(r"P$|^RNU6-\d{2,}|^RNVU", regex=True) & df["paralog_group"].isna()
    return df.sort_values(["chrom", "start"]).reset_index(drop=True)
