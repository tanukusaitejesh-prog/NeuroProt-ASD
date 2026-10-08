# snRNA-VEP manuscript (assembled draft, 9 October 2026)

Assembled from the component drafts; each remains the editable source. Every number traces to a
committed script. Primary label set v5c. Pre-registration: ANALYSIS_PLAN.md.

---

# Abstract v2 — reframed for AJHG (8 October 2026)

Supersedes the abstract in MANUSCRIPT_DRAFT.md (v0.1). Primary label set is now v5c. Every number below is traceable
to a committed script; sources in brackets. Changes from v0.1 are listed at the end.

---

**Title**

Multi-state spliceosome structure predicts the measured functional impact of snRNA variants that conservation,
CADD and population constraint do not

*(alternatives: "A structure-based method for interpreting variants in spliceosomal snRNA genes, calibrated against
saturation genome editing"; "Population constraint and molecular function diverge in snRNA variant interpretation")*

---

**Background.** Variants in spliceosomal small nuclear RNA (snRNA) genes have emerged as a frequent cause of
neurodevelopmental disorders — ReNU syndrome (*RNU4-2*), *RNU2-2*, *RNU5B-1* — and of recessive multisystem disease
(*RNU4ATAC*, *RNU12*, *RNU6ATAC*). Current guidance for classifying these variants states that existing in silico
tools were not designed or calibrated for non-coding RNA genes, and recommends against using them for PP3/BP4 until
they are; it further notes that snRNA structural relationships change through the splicing cycle in ways that "may not
be captured by the analysis of a single structure" (D'Souza et al., 2026).

**Methods.** We built snRNA-VEP, which represents each snRNA nucleotide by its protein and RNA contacts across 11
cryo-EM structures spanning the splicing cycle, projects paralogs and minor-spliceosome analogues onto shared
coordinates, and propagates contact information over a multi-state 3D contact graph. Evidence is combined per disease
mechanism, with the scorer chosen by nested leave-one-gene-out cross-validation. Features, labels, model and gate
criteria were frozen in a timestamped analysis plan before the evaluations below were run [ANALYSIS_PLAN.md].

**Results.** Combining all 11 structures outperformed any single splicing state selected without access to the test
gene (mean held-out AUROC 0.815 vs 0.631; position-block ΔAUROC +0.184, 95% CI 0.126–0.292) and exceeded a
degree-preserving shuffled-contact null in 200/200 replicates (p = 0.005) [22]. On genes held out entirely from
training, snRNA-VEP exceeded CADD by +0.205 AUROC (95% CI 0.106–0.299) and phyloP447 by +0.198 for dominant variants
[16]; in a mixed model with a gene random effect the odds ratio per standard deviation was 3.50 (2.74–4.48) against
1.32 (0.98–1.77) for CADD [24].

Without any RNU4-2 or U4-family labels in training, the model recovered saturation-genome-editing damage in *RNU4-2*
(ρ = 0.405, 95% CI 0.255–0.535; graph-smoothed contacts 0.437, 0.293–0.561), where CADD reached 0.072 and phyloP447
−0.004 [24].

A pre-specified gate failed: gnomAD population constraint matched the model on clinical labels (0.766 vs 0.810; CI of
the difference −0.123 to 0.112) [23]. Against measured function, however, the two separated decisively — constraint
reached ρ 0.206 against the model's 0.405, the paired difference was +0.231 (0.093–0.363, p = 0.002), and partial
correlations were asymmetric: structure retained ρ 0.394 controlling for constraint, while constraint retained 0.025
controlling for structure [27]. Curated pathogenicity labels and molecular function are therefore not interchangeable
benchmarks for non-coding variants, and depletion-based scores track the former without carrying independent
information about the latter.

Trained only on evidence public before May 2025, the model ranked variants that subsequently became ClinVar P/LP
(pooled ΔAUROC vs CADD +0.173, 0.039–0.304, p = 0.009) [17]. Dominant pathogenic variants were enriched at contacts
with the activation and catalytic machinery (PRPF8 OR 6.4, FDR 5×10⁻⁴; RBM42 12.0; SART1 5.8; SNRNP200/Brr2 4.4),
recessive variants at snRNP-assembly contacts (SNU13 k-turn OR 4.4, FDR 0.036) [18a, 18b]; a formal interaction test
of this contrast did not meet its pre-specified criterion and the distinction is reported as a gradient. Recessive
prediction remained limited for every method (AUROC 0.63–0.69) once healthy homozygous carriers from 490,000 UK
Biobank genomes were used as controls.

Finally, population data cannot assess part of this gene family at all. *RNU1-3* and *RNU1-4* have zero callable
alleles across 76,215 gnomAD genomes, yet carry 22 and 15 variants across 146 and 62 of 464 long-read haplotypes in
the HPRC pangenome — as variable as their callable sibling *RNU1-2* [41]. Constraint and de novo enrichment, the two
lines of evidence used to nominate snRNA disease genes, are uncomputable there.

**Conclusions.** Multi-state structural context identifies dominant snRNA variants that conservation and genome-wide
predictors miss, and is the only evidence type tested here that tracks measured molecular function. We provide
precomputed scores for 13,280 variants across 17 spliceosomal snRNA genes, annotated with which loci population
sequencing can and cannot assess.

---

---

---

# Introduction

Variants in spliceosomal small nuclear RNA genes have moved, in under three years, from a curiosity to one of
the more frequent identified causes of neurodevelopmental disorder. RNU4-2 accounts for roughly 0.4–0.5% of
undiagnosed developmental disorder in large cohorts; RNU2-2 has been described as the most prevalent known
recessive neurodevelopmental disorder; RNU5B-1, RNU4ATAC, RNU12 and RNU6ATAC follow. Together these genes now
explain more than one percent of cases that genome sequencing had left unsolved.

Interpreting new variants in them is unresolved. The 2026 guidance for classifying snRNA variants states that
existing in silico tools were not designed or calibrated for non-coding RNA genes and should not be used for
PP3/BP4 until they are. It also notes that the structural relationships within and between snRNAs change through
the splicing cycle, imposing constraints that "may not be captured by the analysis of a single structure".

That second observation is the starting point here. An snRNA is not a static molecule with a fixed set of
contacts; U4 is base-paired to U6 in the tri-snRNP, unwound during activation, and absent from the catalytic
spliceosome entirely. Which nucleotides matter depends on when you look. We therefore represent each nucleotide
by its contacts across eleven cryo-EM structures spanning the cycle, and ask whether that representation
identifies pathogenic variants, predicts measured molecular function, and distinguishes between the alternative
alleles at a position.

A second question runs alongside. Evidence for non-coding variant pathogenicity rests heavily on one source:
population sequencing. Curated labels are assembled partly from depletion; the comparator scores are themselves
built on conservation and allele frequency. Whether those two things — selection against a variant, and what
the variant does to the molecule — are interchangeable has not been tested where both can be measured. snRNA
genes are one of the few places where they can.

---

---

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

---

# Discussion

## What the structural model does and does not do

Across held-out genes, multi-state structural context identifies dominant snRNA variants that conservation and
genome-wide predictors miss, and does so with a margin that survives gene-clustered uncertainty. The multi-state
representation is not decorative: a single splicing state chosen without access to the test gene performs far
worse, and the full model exceeds every one of 200 degree-preserving shuffled-contact nulls. This supports the
structural claim the classification guidance makes on biological grounds.

For recessive variants it does not work. No method tested exceeded 0.69 once healthy homozygous carriers were
used as controls, and the model is no better than conservation. We state this as a limit of position-based
prediction rather than a property of one method: recessive snRNA variants are hypomorphic, cluster in regions
that tolerate heterozygous variation, and are measurably milder in the one functional assay available. The
pre-specified test of a categorical dominant/recessive mechanism split did not meet its criterion, and we
report a gradient instead of a separation.

## Labels and function are not the same benchmark

The most consequential result was a failed gate. Population constraint matched the structural model on clinical
labels, which by the pre-specified rule was a failure. Against measured function the two separate decisively,
and the partial correlations are asymmetric: structure retains almost all of its signal controlling for
constraint, constraint retains almost none controlling for structure.

The explanation is not subtle. Curated pathogenicity labels are assembled partly from population depletion, so a
depletion score predicts them partly for free. A benchmark built on such labels will rank a population-derived
score alongside a mechanistic one while the two carry different information about the molecule. This has been
shown for missense variants, where gene identity dominates ClinVar benchmarks and function-based evaluation
reorders predictors; we arrive second to that framing and say so. What we add is the non-coding case, where the
effect is larger and where, unusually, both label types exist for the same variants.

The practical consequence is narrow and worth stating precisely: for non-coding variants, computational evidence
drawn from conservation or depletion is not independent of population-frequency evidence. Where a classification
framework treats them as separate lines, it is at risk of counting one observation twice.

## Allele resolution is a separate problem from classification

Base-pair disruption predicts which of the three substitutions at a position is most damaging — the first
demonstration of allele-level structural prediction for an snRNA, and something CADD, despite being
allele-resolved, does not achieve on the same data. But adding those features to the clinical model made it
worse. Ranking alleles within a position and discriminating pathogenic from benign across positions are
different tasks, and conflating them costs accuracy. We therefore report the allele score as a complementary
output for the specific clinical question of which substitution at a known position is worst, not as a component
of the classifier.

This result rests on a single dataset. RNU4-2 saturation genome editing is the only snRNA experiment that
measured more than one allele per position; the RNU4ATAC assay has none. Until a second such dataset exists, the
allele result should be treated as promising and unreplicated.

## Part of this gene family cannot be assessed at all

RNU1-3 and RNU1-4 have no callable alleles across 76,215 gnomAD genomes yet carry abundant variation in
long-read assemblies. The cause is not poor coverage but non-identifiability: the four RNU1 copies are identical
across all 164 nucleotides, so no short read can ever be assigned among them. The same holds for U6. Neither
population depletion nor de novo enrichment — the two lines of evidence used to nominate snRNA disease genes —
can be computed at those loci, and a published systematic screen reports zeros there that read as negative
results.

U1 is the only major spliceosomal snRNA family with no established disease gene. We do not claim this is why.
We claim only that the question has not been asked of those loci, cannot be asked with short reads, and can be
asked with long reads.

## Limitations

Two of four dominant genes contribute 7 and 4 pathogenic variants, so the dominant result rests substantially on
RNU4-2 and RNU2-2. The function validation is one assay in one gene. The allele result is unreplicated. The
pangenome comparison uses 232 samples and graph-based genotypes that can err in repetitive sequence. GENCODE v27
models RNU4-2 as 141 nt against the current 145, leaving four 3′ nucleotides unscored. And the model scores
positions for indels, which matters because the single most common ReNU allele is an insertion.

## Conclusion

Structural context across the splicing cycle is, among the evidence types tested here, the only one that tracks
what variants do to the molecule. Population-derived scores track what selection has done to the position, and
on curated labels the two are hard to tell apart. For a gene class where clinical interpretation is explicitly
uncalibrated, that distinction decides what computational evidence is worth.

---

# Methods (draft)

## Pre-registration
Features, label set, model specification and the two gate criteria were fixed in a timestamped analysis plan
(`ANALYSIS_PLAN.md`) committed before any of the evaluations reported here were run. Analyses added afterwards
are labelled exploratory in the text and in the commit history. No hyperparameter, feature or label was changed
after that commit; results that failed their criteria are reported as failures.

## Genes, coordinates and paralog projection
Seventeen spliceosomal snRNA genes were taken from GENCODE v27 (`01_build_catalogue.py`). Paralogs and
minor-spliceosome analogues were projected onto a shared family coordinate by global pairwise alignment to a
family reference (U4 → RNU4-2, U6 → RNU6-1, U5 → RNU5A-1, U2 → RNU2-2P, U1 → RNU1-1), so evidence can transfer
between copies (`snrna_vep/align.py`). Reference sequences are stored with 2,000 bp flanks and are
reverse-complemented for minus-strand genes before use.

## Structural features
Eleven cryo-EM structures spanning the splicing cycle were used (tri-snRNP 3JCR, 6QW6; pre-B 6QX9; B 5O9Z;
Bact 6FF7; Bact-mature 5Z56; C* 5XJC; P 6QDV; U2-snRNP 6Y5Q; minor-Bact 7DVQ; minor-preB 8Y6O). snRNA chains are
identified by sequence alignment against the family references rather than by hand-entered chain IDs, so partially
modelled chains still match (≥85% identity over modelled residues). For each nucleotide in each structure we
record contacting protein residues, minimum protein distance, snRNA–snRNA contacts and contacts with non-snRNA
RNA, using a 4.5 Å heavy-atom cutoff (`snrna_vep/contacts.py`).

Aggregating over states gives: fraction of resolved states with a protein contact, maximum protein residues,
fraction of states with an snRNA contact, and a graph-smoothed contact fraction obtained by one round of
message passing over the multi-state RNA–RNA contact graph (row-normalised adjacency including self-loops and
backbone neighbours).

## Allele-specific features
Watson–Crick pairing is detected geometrically per structure: a purine N1 within 3.5 Å of a pyrimidine N3 is
called a pair and the partner base recorded, so each alternative allele can be classed as preserving a canonical
pair (A-U, G-C), forming a wobble (G-U) or breaking the pair. Protein contacts are split into **base** contacts
(ring and exocyclic atoms, allele-sensitive) and **backbone** contacts (phosphate and ribose, allele-insensitive),
since only base contacts can be perturbed by a substitution (`snrna_vep/allele_features.py`).

**Indels.** Insertions and deletions are scored at the affected node — the first deleted base, or the anchor+1
base for an insertion — and inherit the position-level features of that node. Allele-specific pairing features
are not defined for indels and are not applied to them; indels therefore receive position-level scores only.
This is stated because 70–77% of ReNU cases are a single-base insertion, so the limitation is material.

## Labels
Pathogenic variants were curated from nine supplements plus ClinVar (≥1 star, 2026-10-04), with every row
traceable to its source. Controls are mechanism-aware: for recessive genes only gnomAD homozygotes
(nhomalt ≥ 1) and ClinVar B/LB qualify, because heterozygous carriage is uninformative about recessive
benignity; RNU4-2 uses UK Biobank and All of Us variants, which are independent of gnomAD; dominant-only genes
use gnomAD PASS alleles restricted to the transcribed region. Homozygous variants from ~490,000 UK Biobank
genomes were added as recessive controls (`05h`). The primary label set is **v5c** — v5 with the 16 variants
annotated as recessive-pathogenic yet carried homozygously by a healthy adult removed; v5 is retained as a
sensitivity analysis.

## Model and evaluation
Features enter as unsupervised within-gene percentile ranks computed over all possible variants of each gene, so
no labels are used in feature construction. The model is L2 logistic regression (C = 0.3, class-weighted) on
[features, features × recessive-context, recessive-context]. Evaluation holds out an entire gene (RNU6-1/2/8/9
treated as one unit); a unit and context is evaluated only with ≥3 pathogenic and ≥3 controls. For each held-out
unit the scorer is chosen from four candidates by an inner leave-one-unit-out loop over the **training** units
only, so selection never sees the test gene.

## Uncertainty
Position-block bootstrap resamples snRNA positions within units, so all variants at a position move together;
the hierarchical bootstrap resamples units and then positions. Per-gene AUROCs carry position-block CIs. A
logistic GLMM with a gene random intercept is reported as a secondary analysis. Paired comparisons use
rank-pooling within unit. **P-values are not used in the abstract.**

## Null model for the multi-state gate
Per-structure, per-family contact profiles are permuted among resolved positions, preserving each structure's
degree sequence; the RNA–RNA contact graph is rewired within each structure so every node keeps its degree. The
full model is refitted on each of 200 replicates.

## External data
Saturation genome editing scores for RNU4-2 from De Jonghe et al. (2026); the RNU4ATAC cellular assay and
published RNAstructure scores from Benoit-Pilven et al. (2020), whose numbering (NR_023343.1) is 1 nt upstream
of the GENCODE model and is shifted accordingly, validated on reference bases at 234/241 SNVs. HPRC Release 2
(v2.1) Minigraph-Cactus pangenome, queried by HTTP range requests against the tabix index
(`snrna_vep/remote_tabix.py`).

## Code and data
All analyses are scripted and committed; each result in the text names its script. Precomputed scores for
13,280 variants in 17 genes are provided, annotated with which loci population sequencing can assess.

---

# Figure plan

The strongest result is the agreement with measured function, so it leads. Benchmark tables move to supplement;
a reviewer who wants AUROCs will find them, but they are not the argument.

---

## Figure 1 — The multi-state contact model
**a.** Schematic: U4/U6 duplex through the splicing cycle, with the 11 cryo-EM structures placed on the cycle
(tri-snRNP → pre-B → B → Bact → C* → P, plus minor-spliceosome states). Shows directly why no single structure
covers every gene: Bact/C*/P lack U4 entirely.
**b.** Contact profile of RNU4-2 across states — per-nucleotide fraction of states with a protein contact, with
the 18-nt critical region marked. Pathogenic variants overlaid.
**c.** Gate G1: mean held-out AUROC for ALL 11 structures, each single state, nested single state, and the
degree-preserving shuffled null (violin of 200 replicates with the observed value marked). *[22]*

## Figure 2 — Agreement with saturation genome editing *(the lead result)*
**a.** Scatter: model score (U4 family withheld from training) vs SGE function score, 485 variants, with ρ and
position-block CI. *[24, 27]*
**b.** The same for CADD and phyloP447 on identical variants — visually flat.
**c.** Bar: fraction of the 121 most functionally damaging variants placed in each score's own top quartile —
structure 59%, density 45%, constraint 40%, CADD 21%, phyloP 11%. *[27]*

## Figure 3 — Population constraint and molecular function diverge
**a.** Paired: constraint vs the model on *clinical labels* (overlapping, gate G2 failure) and on *measured
function* (separated). Makes the dissociation visible in one panel.
**b.** Partial correlations: structure given constraint +0.394; constraint given structure +0.025. *[27]*
**c.** Interpretation panel: curated labels derive partly from population data, so depletion scores them for
free. Schematic, not data.

## Figure 4 — Allele-specific resolution
**a.** Base-pairing geometry across states for one worked position, showing canonical / wobble / broken
outcomes for each alternative allele.
**b.** Within-position ranking against −SGE, 137 positions with all three alleles measured: allele score
+0.346, pair-disruption +0.374, CADD −0.077, position-level controls informative at zero positions. *[44]*
**c.** Coverage honesty: the allele score varies at 45 of 137 positions; elsewhere the method is position-level.

## Figure 5 — Mechanism
**a.** Mantel–Haenszel odds ratios by contacted protein, dominant vs recessive, BH-FDR. PRPF8 6.4, RBM42 12.0,
SART1 5.8, SNRNP200 4.4 for dominant; SNU13 4.4 for recessive. *[18a, 18b]*
**b.** The same variants mapped onto the tri-snRNP structure, coloured by mechanism.
**c.** Stated plainly in the legend: the pre-specified interaction test did not meet its criterion, so this is a
gradient, not a categorical separation.

## Figure 6 — The limits of the evidence base
**a.** gnomAD callable allele number per snRNA gene; RNU1-3 and RNU1-4 at zero, RNU2-1 at 15%. *[e372a5f]*
**b.** The same loci in the HPRC pangenome: 22 and 15 variants across 146 and 62 of 464 haplotypes. *[41]*
**c.** Pairwise paralog identity within families with the fraction of 101-nt reads that can be assigned — U4/U5
100%, U2 92%, U6 12%, **U1 0%**. The point: U1 copies are sequence-identical, so this is non-identifiability,
not low coverage. *[43]*

---

## Supplementary figures
- S1 held-out AUROC per gene with position-block CIs, all scorers *[16, 24]*
- S2 extended baselines: ViennaRNA, trained secondary structure, Rfam elements, population density,
  paralog density, distance to critical region *[23]*
- S3 prospective ClinVar and the RNU4-2-only time split *[17, 12]*
- S4 leave-one-state-out and the graph-smoothing ablation *[22]*
- S5 RNU4ATAC external comparison, including that conservation beats structure there *[25]*
- S6 recessive performance across all methods (the negative) *[16]*
- S7 allele features added to the clinical model — the −0.031 negative *[46]*
- S8 GLMM with gene random intercept *[24]*
- S9 minor-spliceosome exclusion *[24]*
- S10 coding-genome control: 89 of 34,556 MANE genes with synonymous deficit, 2 clinically relevant *[42]*

## Tables
- T1 label set composition by gene, mode and source
- T2 pre-specified gates: criterion, result, pass/fail
- T3 the resource: 13,280 variants, 17 genes, with per-locus assessability flags

---

# What is new here, against the directly competing literature

Every claim in the right-hand column is traceable to a committed script. Where a competing paper already
establishes something, that is stated plainly rather than minimised.

## 1. Against the papers that defined these disease genes

| paper | what it established | what it did not do | this work |
|---|---|---|---|
| Chen et al., *Nature* 2024 (ReNU / RNU4-2) | RNU4-2 as a frequent NDD gene; an 18-nt critical region defined by UK Biobank depletion; RNU4-2 more highly expressed in brain than RNU4-1 | no variant-effect model; screened 28 snRNA genes with de novo enrichment but could not assess loci with no callable alleles | a structure-based score for all variants in 17 genes; quantifies which loci their screen could not see (RNU1-3, RNU1-4: zero callable alleles in gnomAD but 22 and 15 variants in long-read assemblies) |
| Nava et al., *Nat Genet* 2025 (U4/U5 NDD) | RNU5B-1 as an NDD gene; variants cluster in stem III and T-loop/quasi-pseudoknot; severity differs by domain; 5′SS usage altered | domains assigned by inspection of a single structure; no predictive model; the SVM is a methylation episignature classifier, not a variant predictor | contacts computed across 11 cryo-EM states; multi-state beats any single state chosen without the test gene (0.815 vs 0.631, p = 0.005 against a degree-preserving null) |
| De Jonghe et al., *Nature* 2026 (RNU4-2 SGE) | saturation functional scores for 501 variants; function tracks severity; a distinct recessive disorder in the 5′ stem-loop | function measured, not predicted; no model that generalises to other genes | predicts those scores **zero-shot**, with RNU4-2 and the whole U4 family withheld (ρ 0.405; CADD 0.072, phyloP −0.004) |
| Greene et al. / Nat Genet 2026 (RNU2-2) | RNU2-2 as a frequent recessive NDD gene; systematic screen of 200 snRNA genes | excluded 43 genes as unresolvable by short reads and called for long-read analysis; that analysis was not done | provides it: which snRNA loci are resolvable *in principle*, derived from sequence identity and confirmed against the HPRC pangenome |
| Benoit-Pilven et al., *PLoS One* 2020 (RNU4ATAC) | the only previously published structural score for an snRNA gene (RNAstructure bifold on the U4atac:U6atac duplex) | poor agreement with its own cellular assay (ρ 0.046) | direct comparison on the same 23 variants; reported honestly, including that **conservation beats structure for this recessive gene** (CADD ρ 0.713) |
| D'Souza et al., medRxiv 2026 (classification guidance) | states that no tool was "designed or calibrated specifically for variants in non-coding RNA genes", and that snRNA structure changes through the cycle in ways that "may not be captured by the analysis of a single structure" | issues the requirement; supplies no tool | the multi-state model is a direct response, with the gate G1 result as evidence for the guidance's own structural claim |

## 2. Against the general variant-interpretation literature

| prior work | relation |
|---|---|
| "Gene identity, not variant effect, dominates ClinVar benchmarks" (bioRxiv, Aug 2026) | **precedent, acknowledged.** Shows for *missense* variants that ClinVar benchmarks reflect gene identity and that function-based evaluation reorders predictors. This work is the non-coding analogue, and arrives second to the framing. What is added is the mechanism — population constraint matches a mechanistic model on labels while contributing ρ 0.025 beyond it on measured function — and the demonstration in a class where labels and function coexist. |
| MitoTIP (2017) | precedent for paralog/position transfer in a non-coding RNA; implemented here as a baseline (paralog density), which the model exceeds by +0.195 AUROC |
| 5ULTRA (*AJHG* 2026) | closest structural analogue: an underserved non-coding variant class, a principled feature model, validated against measured effects. Establishes the venue and format rather than competing on content |

## 3. Claims that are new here, with their limits

| claim | evidence | limit |
|---|---|---|
| Multi-state structural context beats any single splicing state | 0.815 vs 0.631; exceeds 200/200 shuffled-contact nulls, p = 0.005 | hierarchical CI crosses zero with only 4 dominant genes; ties the *oracle* single state |
| Structure predicts measured function where conservation, CADD and population constraint do not | ρ 0.405 vs 0.072 / −0.004 / 0.206; partials +0.394 vs +0.025 | one gene, one assay |
| Within-position allele ranking from base-pair disruption | mean within-position ρ +0.346, p = 0.001; CADD −0.077 | works only at base-paired positions (45 of 137); **no second dataset exists** — RNU4ATAC has 0 of 23 positions with multiple alleles |
| Parts of the gene family are formally unmeasurable | RNU1-1/2/3/4 are 100% identical across 164 nt; gnomAD AN = 0; 22 and 15 variants in long-read assemblies | 232 pangenome samples; graphs can err in repeats |
| Allele features do not improve the clinical model | −0.031 AUROC dominant, p < 0.001 | reported as a negative; separates allele ranking from case/control discrimination |

## 4. Pre-specified analyses that failed, and are reported as failures

- **Gate G2.** Population constraint matched the model on clinical labels (0.766 vs 0.810, CI spanning zero). Reported, with the function analysis as its resolution.
- **Dominant/recessive interaction.** Joint Wald p = 0.022 but with the contact term in the wrong direction and RNU2-2 flat (p = 0.895). Criterion not met; reported as a gradient.
- **Recessive prediction.** 0.591 rank-pooled, no better than CADD or phyloP. Stated as a limit, not a result.
- **BP4_supporting.** Did not survive the hierarchical bootstrap; dropped.
