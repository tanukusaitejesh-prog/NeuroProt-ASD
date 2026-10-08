# Results (draft, 8 October 2026)

Every number traces to a committed script; the script number is in brackets. Primary label set is v5c
(v5 with the 16 UK Biobank homozygote conflicts removed); v5 is retained as a sensitivity analysis. Features,
labels, model and the two gate criteria were frozen in a timestamped plan before any of the evaluations below
were run (ANALYSIS_PLAN.md, commit 350165e).

---

## 1. A mechanism-aware benchmark of snRNA variants

We assembled 266 pathogenic variants across 10 spliceosomal snRNA genes from 14 sources, with mechanism-aware
controls: for recessive genes only gnomAD homozygotes and ClinVar B/LB count as controls, because heterozygous
carriage carries no information about recessive benignity; RNU4-2 uses UK Biobank and All of Us variants, which
are independent of gnomAD; dominant-only genes use gnomAD PASS alleles. Homozygous variants from ~490,000 UK
Biobank genomes [05h] added 71 recessive controls in RNU2-2, 16 in RNU12 and 14 in RNU4ATAC. Sixteen variants
annotated as recessive-pathogenic were carried homozygously by at least one healthy UK Biobank adult; these are
retained in v5 and removed in v5c, and we report both.

## 2. Multi-state structure outperforms any single splicing state

Contact profiles were recomputed from scratch for each structure subset using the main-pipeline definitions, so
that configurations differ only in which structures contribute [20, 22]. Combining all 11 cryo-EM structures gave
a mean held-out AUROC of 0.815 for dominant variants, against 0.631 for the best single splicing state chosen by
nested leave-one-unit-out over the training genes only (position-block ΔAUROC +0.184, 95% CI 0.126–0.292;
better in 3 of 4 genes). No single state suffices because the genes need different ones: RNU4-2 requires a
pre-activation structure, since every Bact/C*/P structure lacks U4 entirely, whereas RNU2-2 and RNU5B-1 are best
resolved by activated and catalytic states. Only the multi-state model exceeds 0.73 on all four dominant genes.

Against a degree-preserving null — per-position contact profiles permuted among resolved positions within each
structure and family, and the RNA–RNA contact graph rewired to preserve every node's degree, with the full model
refitted on each of 200 replicates — the observed value of 0.815 exceeded every replicate (null mean 0.557, 95th
percentile 0.689; empirical p = 0.005) [22]. The gate's three pre-specified criteria were all met.

Leave-one-state-out showed no single state is indispensable (0.796–0.831), and graph smoothing contributed
+0.018. Removing the minor-spliceosome genes from training left the dominant result unchanged (0.810 → 0.813),
excluding a family-projection artefact [24].

## 3. Held-out performance against established predictors

On genes held out entirely from training, and with the scorer for each mechanism chosen by an inner
leave-one-unit-out loop over training genes only [11, 16]:

| held-out gene | n pathogenic | snRNA-VEP | contacts | phyloP447 | CADD |
|---|---|---|---|---|---|
| RNU2-2 (dominant) | 24 | **0.807** [0.681, 0.910] | 0.640 | 0.674 | 0.462 |
| RNU4-2 (dominant) | 37 | **0.785** [0.700, 0.860] | 0.704 | 0.502 | 0.584 |
| RNU5B-1 (dominant) | 7 | **0.736** [0.580, 0.950] | 0.477 | 0.705 | 0.920 |
| RNU6 (dominant, RP) | 4 | **0.912** [0.840, 0.970] | 0.930 | 0.491 | — |

Rank-pooled across units, the final system exceeded CADD by +0.205 AUROC (95% CI 0.106–0.299, p < 0.001) and
phyloP447 by +0.198 (0.113–0.282, p < 0.001) for dominant variants; it also exceeded the single best contact
feature by +0.106 (0.054–0.159, p < 0.001), showing the trained combination adds to its strongest input [16].
In a logistic model with a gene random intercept, the odds ratio per standard deviation of score was 3.50
(2.74–4.48), against 1.32 (0.98–1.77) for CADD [24].

For recessive variants the system did not beat conservation (ΔAUROC vs phyloP +0.010, p = 0.74; vs CADD +0.023,
p = 0.51), and no method exceeded 0.69 once healthy homozygous carriers were used as controls [16]. We report
this as a limit of position-based prediction rather than a property of any one method.

## 4. Agreement with saturation genome editing

With RNU4-2 and the entire U4 family excluded from training, scores were compared with saturation-genome-editing
function scores for 485 RNU4-2 variants [24, 27]. Spearman ρ with −SGE was 0.405 (position-block 95% CI
0.255–0.535) for the model and 0.437 (0.293–0.561) for the graph-smoothed contact feature alone, against 0.072
for CADD and −0.004 for phyloP447. Within the dominant critical region the model reached 0.575 (0.310–0.732).

## 5. Population constraint and molecular function diverge

A pre-specified gate failed. gnomAD-derived population constraint matched the model on clinical labels (mean
held-out AUROC 0.766 vs 0.810; position-block CI of the difference −0.123 to +0.112), so on curated labels the
two are not separable [23].

Against measured function they separate decisively [27]. On the same 485 SGE variants, constraint reached
ρ 0.206 (0.035–0.344) and per-position variant density 0.258, against 0.405 for the model; the paired difference
was +0.231 (0.093–0.363, p = 0.002), and every paired comparison excluded zero. Partial correlations are
asymmetric: structure retains ρ 0.394 controlling for constraint, whereas constraint retains 0.025 controlling
for structure. Among the 121 most functionally damaging variants by assay, the structural score places 59% in its
own top quartile, against 45% for density, 40% for constraint, 21% for CADD and 11% for phyloP.

Curated pathogenicity labels are themselves derived in part from population data, which explains why a depletion
score matches a mechanistic model on labels while carrying almost no independent information about function.
For non-coding variants the two are therefore not interchangeable benchmarks.

## 6. Allele-specific structural scoring

The contact model scores positions, so all three substitutions at a nucleotide receive one score. We added
allele-resolved features from the same structures [44]: Watson–Crick pairing geometry per state with partner
identity, so each alternative allele can be classed as preserving a canonical pair, forming a wobble, or breaking
the pair; and protein contacts split into base contacts (allele-sensitive) versus backbone contacts (not).
Across structures, 54% of modelled nucleotides are paired in at least one state, 31% have a base–protein contact
and 22% have backbone contacts only.

These were tested where the answer is measurable: within positions of RNU4-2 at which SGE assayed all three
alternative alleles (137 positions with structural features). A position-level score has, by construction, no
within-position signal, so any correlation is attributable to allele identity. The allele score reached a mean
within-position Spearman of +0.346 with −SGE (Wilcoxon p = 0.0010), and pair disruption alone +0.374; CADD, which
is allele-resolved, reached −0.077. The most damaging allele was ranked first 44–46% of the time against 33%
expected by chance, and 32% for CADD. Both position-level controls were informative at zero positions, as
expected.

The allele score varies only where the nucleotide is base-paired, which is 45 of 137 tested positions. The method
is therefore variant-level at paired positions and position-level elsewhere.

## 7. Prospective and temporal validation

Trained only on evidence public before the ClinVar release of 4 May 2025 and only on other genes, the system
ranked the 41 variants that subsequently became P/LP (pooled ΔAUROC vs CADD +0.173, 0.039–0.304, p = 0.009;
vs phyloP +0.111, p = 0.11) [17]. Trained on RNU4-2 alone, it transferred above chance to every later-discovered
gene [12].

## 8. Mechanism

Mantel–Haenszel odds ratios stratified by gene, at position level with BH-FDR [18a, 18b], placed dominant
pathogenic variants at contacts with the activation and catalytic machinery — PRPF8 OR 6.4 (FDR 5×10⁻⁴),
RBM42 12.0 (0.003), SART1 5.8 (0.003), SNRNP200/Brr2 4.4 (0.006), SNRNP27 4.9 (0.009) — with the Sm ring not
enriched (0.8). The RBM42 signal independently recovers the region reported for severe ReNU. Recessive variants
were enriched at snRNP-assembly contacts (SNU13 k-turn OR 4.4, FDR 0.036).

A pre-specified interaction test of this contrast did not meet its criterion: the joint Wald test over the three
feature × mechanism terms gave p = 0.022, but the contact term ran opposite to the predicted direction and
RNU2-2 alone was flat (p = 0.895) [24]. We therefore report a gradient rather than a categorical separation.

## 9. What the evidence base cannot see

Population data cannot assess part of this gene family at all. RNU1-3 and RNU1-4 have zero callable alleles
across 76,215 gnomAD genomes, while their near-identical siblings RNU1-1 and RNU1-2 are at full cohort depth
[e372a5f]. In the HPRC Release 2 pangenome, built from long-read haplotype assemblies, the same two genes carry
22 and 15 variants across 146 and 62 of 464 haplotypes — RNU1-3 matching RNU1-2 exactly at 134 variants per kb
[41]. The absence is ascertainment, not invariance, and neither constraint nor de novo enrichment can be
computed there.

The same method distinguishes real absence from invisibility: of five genes returning zero variants of every
class in a published systematic screen, three (RNVU1-1 and two U7 genes) are genuinely near-invariant in long-read
assemblies, and only RNU1-3 and RNU1-4 are not. RNU4-2 itself shows 6 variants across 6 haplotypes, consistent
with genuine constraint at a confirmed disease gene.

Across coding sequence a comparable artefact exists but is small and already recognised: 89 of 34,556 MANE genes
(0.26%) are missing at least half their expected synonymous variation, and only 2 carry any ClinVar P/LP variant
[42]. Coding sequence has an internal control — synonymous variation — that non-coding sequence lacks.

## 10. Resource

Scores for 13,280 variants across 17 spliceosomal snRNA genes are provided, with within-gene percentiles,
features, and an annotation of which loci population sequencing can and cannot assess [19, 41].

---

## Still to add
- AlphaGenome as a baseline (requires an API key; the UCSC tracks are licence-restricted)
- Novelty table against the SGE, RNU2-2 and RNU4ATAC papers
- Figure plan (lead with §4, the SGE validation)
- Methods: how indels are scored
