"""Context-aware observed/expected SNV depletion along each core snRNA gene (gnomAD v4.1 genomes).

Sanity check: without any patient labels, the RNU4-2 ReNU critical region (T-loop / stem III,
around n.62-79) should come out as the most depleted window.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from snrna_vep import config
from snrna_vep.coords import g_to_n
from snrna_vep.mutation_model import fit, observed_expected, site_table, windows_depleted

WINDOW = 10
EXCLUDE_LOW_QUALITY = {"RNU1-3", "RNU1-4", "RNU2-1"}   # <10 PASS calls in gene: multicopy arrays
BUFFER = 100                                         # keep promoter/3' box out of the neutral flank

genes = pd.read_csv(config.PROCESSED / "core_genes.tsv", sep="\t")
gn = pd.read_csv(config.PROCESSED / "gnomad_v4.1_core_loci.tsv.gz", sep="\t")
refs = pd.read_csv(config.PROCESSED / "reference_core_loci.tsv.gz", sep="\t").set_index("gene_name")

snv = gn[(gn["filter"] == "PASS") & (gn.ref.str.len() == 1) & (gn.alt.str.len() == 1) & ~gn.segdup.astype(bool)]

flank_sites, gene_sites = [], {}
for g in genes.itertuples():
    if g.gene_name in EXCLUDE_LOW_QUALITY:
        continue
    r = refs.loc[g.gene_name]
    obs = set(zip(snv[snv.gene_name == g.gene_name].pos, snv[snv.gene_name == g.gene_name].alt))
    s = site_table(r.chrom, r.start, r.end, r.seq, obs)
    in_gene = (s.pos >= g.start) & (s.pos <= g.end)
    near = (s.pos >= g.start - BUFFER) & (s.pos <= g.end + BUFFER)
    flank_sites.append(s[~near])
    gene_sites[g.gene_name] = (g, s[in_gene].copy())

model = fit(pd.concat(flank_sites))
model.to_csv(config.PROCESSED / "context_model_flanks.tsv", sep="\t")
print(f"context model: {len(model)} classes, mean P(observed)={model.mean():.3f}")

summary, per_rows = [], []
for name, (g, s) in gene_sites.items():
    per = observed_expected(s, model)
    per["n_pos"] = [g_to_n(p, g.start, g.end, g.strand) for p in per.index]
    per = per.sort_values("n_pos").set_index("n_pos")
    tot_o, tot_e = per.obs.sum(), per.exp.sum()
    if tot_o < 20:   # not callable (segdup / multicopy); nothing to rescale against
        print(f"skip {name}: {int(tot_o)} usable PASS SNVs")
        continue
    # snRNA genes are hypermutable relative to their flanks, so rescale expectation to the
    # gene's own average rate; depletion is then *within-gene* relative constraint.
    per["exp_gene"] = per.exp * tot_o / tot_e
    w = windows_depleted(per.rename(columns={"exp": "exp_flank", "exp_gene": "exp"}), WINDOW)
    per_rows.append(per.assign(gene_name=name).reset_index())
    best = w.oe.idxmin() if len(w) else None
    summary.append(dict(gene_name=name, length=g.length, obs=int(tot_o), exp=round(tot_e, 1),
                        gene_oe=round(tot_o / tot_e, 3),
                        min_window_center=best, min_window_rel_oe=round(w.oe.min(), 3) if len(w) else None,
                        min_window_p=w.p.min() if len(w) else None,
                        n_sig_windows=int(w.significant.sum()),
                        disease_modes=g.disease_modes))

summ = pd.DataFrame(summary).sort_values("min_window_p")
config.RESULTS.mkdir(exist_ok=True)
summ.to_csv(config.RESULTS / "depletion_summary.tsv", sep="\t", index=False)
pd.concat(per_rows).to_csv(config.RESULTS / "depletion_per_position.tsv", sep="\t", index=False)
print(summ.to_string(index=False))

# Figure: per-position O/E (10-nt window) for disease genes
show = [n for n in ["RNU4-2", "RNU2-2P", "RNU5B-1", "RNU5A-1", "RNU4ATAC", "RNU12", "RNU4-1", "RNU6-1"] if n in gene_sites]
fig, axes = plt.subplots(len(show), 1, figsize=(9, 1.6 * len(show)), sharex=False)
for ax, name in zip(axes, show):
    per = pd.concat(per_rows).query("gene_name == @name").set_index("n_pos")
    w = per[["obs", "exp_gene"]].rolling(WINDOW, center=True, min_periods=WINDOW).sum()
    ax.plot(w.index, w.obs / w.exp_gene, lw=1.5)
    ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.set_ylabel("rel. O/E", fontsize=8)
    ax.set_ylim(0, 1.6)
    ax.set_title(name, fontsize=9, loc="left")
    if name == "RNU4-2":
        ax.axvspan(62, 79, color="tab:red", alpha=0.15, lw=0)
axes[-1].set_xlabel("transcript position (n.)")
fig.tight_layout()
config.FIGURES.mkdir(exist_ok=True)
fig.savefig(config.FIGURES / "depletion_profiles.png", dpi=150)
print("figure saved")
