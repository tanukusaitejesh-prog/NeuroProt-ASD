# MuST-VEP: multi-state spliceosome structure for variant effect prediction in spliceosomal snRNA genes

*Draft 1, 9 October 2026. Target: Briefings in Bioinformatics, Problem Solving Protocol.*
*Every number traces to a committed script; the script number is given in brackets. Primary label set v5c.*
*Pre-registration: `ANALYSIS_PLAN.md`, commit 350165e, timestamped before any evaluation reported here.*

---

## Abstract

Variants in spliceosomal small nuclear RNA (snRNA) genes have become a frequent identified cause of
neurodevelopmental disorder, yet the 2026 classification guidance states that existing in silico tools were not
designed or calibrated for non-coding RNA genes and should not be used for PP3/BP4 until they are. The same
guidance notes that snRNA structural relationships change through the splicing cycle in ways that may not be
captured by analysing a single structure. We address that observation directly. MuST-VEP represents each snRNA
nucleotide by its protein and RNA contacts across eleven cryo-EM structures spanning the splicing cycle, projects
paralogues and minor-spliceosome analogues onto shared family coordinates, and propagates contact information
over a multi-state three-dimensional contact graph; evidence is combined per disease mechanism with the scorer
selected by nested leave-one-gene-out cross-validation. Against thirteen baselines under an identical held-out-gene
protocol, MuST-VEP reached a mean AUROC of 0.810 for dominant variants, exceeding CADD (0.655), phyloP447 (0.593)
and the RNA foundation model RNA-FM (0.460, below chance). Combining all eleven structures outperformed any single
splicing state chosen without access to the test gene (0.815 versus 0.631) and exceeded all 200 degree-preserving
shuffled-contact nulls. A pre-specified gate failed: gnomAD population constraint was not separable from the model
on curated clinical labels (0.766 versus 0.810, CI of the difference −0.123 to +0.112). Against measured function,
however, the two separated decisively. On 485 RNU4-2 saturation-genome-editing variants scored zero-shot with the
entire U4 family withheld, MuST-VEP reached Spearman 0.405 against 0.206 for constraint, 0.072 for CADD and −0.004
for phyloP447; partial correlations were asymmetric, structure retaining 0.394 controlling for constraint while
constraint retained 0.025 controlling for structure. Curated pathogenicity labels and molecular function are
therefore not interchangeable benchmarks for non-coding variants. Two obvious extensions — conformational-trajectory
features and protein-partner identity — were specified with kill criteria, tested, and rejected; both are reported.
We provide scores for 13,280 variants across 17 genes, annotated with which loci population sequencing can and
cannot assess.

**Keywords:** spliceosome, small nuclear RNA, variant effect prediction, cryo-EM, neurodevelopmental disorder,
saturation genome editing, ACMG calibration

---

## Introduction

Variants in spliceosomal small nuclear RNA genes have moved, in under three years, from a curiosity to one of the
more frequent identified causes of neurodevelopmental disorder. *RNU4-2* accounts for roughly 0.4–0.5% of
undiagnosed developmental disorder in large cohorts; *RNU2-2* has been described as the most prevalent known
recessive neurodevelopmental disorder; *RNU5B-1*, *RNU4ATAC*, *RNU12* and *RNU6ATAC* follow. Together these genes
now explain more than one percent of cases that genome sequencing had left unsolved.

Interpreting new variants in them is unresolved. The 2026 guidance for classifying snRNA variants states that
existing in silico tools were not designed or calibrated for non-coding RNA genes and should not be used for
PP3/BP4 until they are. It also notes that the structural relationships within and between snRNAs change through
the splicing cycle, imposing constraints that "may not be captured by the analysis of a single structure".

That second observation is the starting point here. An snRNA is not a static molecule with a fixed set of contacts.
U4 is base-paired to U6 in the tri-snRNP, unwound during activation, and absent from the catalytic spliceosome
entirely; U5 loop I holds the exons through both transesterifications; U2 recognises the branch point early. Which
nucleotides matter depends on when you look. The spliceosome is, moreover, the one large ribonucleoprotein machine
whose entire catalytic cycle has been resolved at near-atomic resolution, so this is a question that can actually
be asked. We therefore represent each nucleotide by its contacts across eleven cryo-EM structures spanning the
cycle, and ask whether that representation identifies pathogenic variants, predicts measured molecular function,
and distinguishes between the alternative alleles at a position.

A second question runs alongside, and it turned out to matter more. Evidence for non-coding variant pathogenicity
rests heavily on one source: population sequencing. Curated labels are assembled partly from depletion; the
comparator scores are themselves built on conservation and allele frequency. Whether those two things — selection
against a variant, and what the variant does to the molecule — are interchangeable has not been tested where both
can be measured. snRNA genes are one of the few places where they can, because *RNU4-2* has a saturation genome
editing screen that assays function without consulting allele frequency at all.

Three design commitments distinguish this work from the comparable methods literature. First, features, labels,
model specification and two gate criteria were frozen in a timestamped analysis plan committed before any
evaluation was run; results that failed their criteria are reported as failures, including the one that failed.
Second, evaluation holds out an entire gene — strictly harder than the cluster splits usual in this literature,
because no variant, position or paralogue of the test gene is seen in training. Third, every claim about structure
is tested against a degree-preserving null in which contact profiles are permuted but every node keeps its degree,
so a result cannot be produced by the amount of contact alone.

---

**Figures.** Figure 1 — the MuST-VEP workflow, the state × family coverage matrix, and gate G1 against the
degree-preserving null (`figures/fig1_workflow.png`, `scripts/54`). Figure 2 — the dissociation between clinical
labels and measured function (`figures/fig2_dissociation.png`, `scripts/55`). Supplementary figures: held-out
AUROC per gene, the *RNU4-2* contact profile, and the per-protein mechanism analysis.

---

## Materials and methods

### Pre-registration

Features, label set, model specification and the two gate criteria were fixed in a timestamped analysis plan
(`ANALYSIS_PLAN.md`) committed before any of the evaluations reported here were run. Analyses added afterwards are
labelled exploratory in the text and in the commit history. No hyperparameter, feature or label was changed after
that commit.

### Genes, coordinates and paralogue projection

Seventeen spliceosomal snRNA genes were taken from GENCODE v27 [01]. Paralogues and minor-spliceosome analogues
were projected onto a shared family coordinate by global pairwise alignment to a family reference (U4 → *RNU4-2*,
U6 → *RNU6-1*, U5 → *RNU5A-1*, U2 → *RNU2-2P*, U1 → *RNU1-1*), so that evidence can transfer between copies
(`snrna_vep/align.py`). Reference sequences are stored with 2,000 bp flanks and are reverse-complemented for
minus-strand genes before use.

### Structural features

Eleven cryo-EM structures spanning the splicing cycle were used: tri-snRNP 3JCR and 6QW6; pre-B 6QX9; B 5O9Z;
Bact 6FF7; Bact-mature 5Z56; C* 5XJC; P 6QDV; U2-snRNP 6Y5Q; minor-Bact 7DVQ; minor-preB 8Y6O. snRNA chains are
identified by sequence alignment against the family references rather than by hand-entered chain identifiers, so
partially modelled chains still match (≥85% identity over modelled residues). For each nucleotide in each structure
we record contacting protein residues, minimum protein distance, snRNA–snRNA contacts and contacts with non-snRNA
RNA, using a 4.5 Å heavy-atom cutoff (`snrna_vep/contacts.py`).

Aggregating over states gives four features: the fraction of resolved states with a protein contact, the maximum
number of protein residues contacted, the fraction of states with an snRNA contact, and a graph-smoothed contact
fraction obtained by one round of message passing over the multi-state RNA–RNA contact graph (row-normalised
adjacency including self-loops and backbone neighbours).

### Allele-specific features

Watson–Crick pairing is detected geometrically per structure: a purine N1 within 3.5 Å of a pyrimidine N3 is called
a pair and the partner base recorded, so each alternative allele can be classed as preserving a canonical pair
(A–U, G–C), forming a wobble (G–U) or breaking the pair. Protein contacts are split into **base** contacts (ring
and exocyclic atoms, allele-sensitive) and **backbone** contacts (phosphate and ribose, allele-insensitive), since
only base contacts can be perturbed by a substitution (`snrna_vep/allele_features.py`).

Insertions and deletions are scored at the affected node — the first deleted base, or the anchor+1 base for an
insertion — and inherit the position-level features of that node. Allele-specific pairing features are not defined
for indels and are not applied to them. This is stated because 70–77% of ReNU cases are a single-base insertion, so
the limitation is material.

### Labels

Pathogenic variants were curated from nine supplements plus ClinVar (≥1 star, 2026-10-04): 266 variants across 10
genes from 14 sources, every row traceable to its source. Controls are mechanism-aware. For recessive genes only
gnomAD homozygotes (nhomalt ≥ 1) and ClinVar B/LB qualify, because heterozygous carriage is uninformative about
recessive benignity; *RNU4-2* uses UK Biobank and All of Us variants, which are independent of gnomAD; dominant-only
genes use gnomAD PASS alleles restricted to the transcribed region. Homozygous variants from approximately 490,000
UK Biobank genomes were added as recessive controls [05h]. The primary label set is **v5c** — v5 with the 16
variants annotated as recessive-pathogenic yet carried homozygously by a healthy adult removed; v5 is retained as a
sensitivity analysis.

### Model and evaluation

Features enter as unsupervised within-gene percentile ranks computed over all possible variants of each gene, so no
labels are used in feature construction. The model is L2 logistic regression (C = 0.3, class-weighted) on
[features, features × recessive-context, recessive-context]. Evaluation holds out an entire gene (*RNU6-1/2/8/9*
treated as one unit); a unit and context is evaluated only with ≥3 pathogenic and ≥3 controls. For each held-out
unit the scorer is chosen from four candidates by an inner leave-one-unit-out loop over the **training** units only,
so selection never sees the test gene.

### Baselines

Thirteen competitors were scored under the identical protocol [23, 53]: gnomAD constraint (1 − observed/expected in
10-nt windows), gnomAD allele density at the position and within ±2 nt, a trained secondary-structure model (same
design and regularisation, contact features replaced by ViennaRNA features), ViennaRNA fold ΔΔG, pairing probability
and duplex ΔΔG, Rfam element-type enrichment, paralogue pathogenic density, distance to the nearest known pathogenic
position, CADD (PHRED), phyloP447, and RNA-FM scored zero-shot under the standard language-model convention
(damage = −[log P(alt) − log P(ref)]).

### Uncertainty

The position-block bootstrap resamples snRNA positions within units, so all variants at a position move together;
the hierarchical bootstrap resamples units and then positions. Per-gene AUROCs carry position-block CIs. A logistic
GLMM with a gene random intercept is reported as a secondary analysis. Paired comparisons use rank-pooling within
unit.

### Null model

Per-structure, per-family contact profiles are permuted among resolved positions, preserving each structure's degree
sequence; the RNA–RNA contact graph is rewired within each structure so every node keeps its degree. The full model
is refitted on each of 200 replicates.

### External data

Saturation genome editing scores for *RNU4-2* from De Jonghe et al. (2026); the *RNU4ATAC* cellular assay and
published RNAstructure scores from Benoit-Pilven et al. (2020), whose numbering (NR_023343.1) is 1 nt upstream of
the GENCODE model and is shifted accordingly, validated on reference bases at 234/241 SNVs. HPRC Release 2 (v2.1)
Minigraph-Cactus pangenome, queried by HTTP range requests against the tabix index (`snrna_vep/remote_tabix.py`).

---

## Results and discussion

### Multi-state structure outperforms any single splicing state

Contact profiles were recomputed from scratch for each structure subset using the main-pipeline definitions, so
that configurations differ only in which structures contribute [20, 22]. Combining all 11 cryo-EM structures gave a
mean held-out AUROC of 0.815 for dominant variants, against 0.631 for the best single splicing state chosen by
nested leave-one-unit-out over the training genes only (position-block ΔAUROC +0.184, 95% CI 0.126–0.292; better in
3 of 4 genes).

No single state suffices because the genes need different ones (Figure 1b). *RNU4-2* requires a pre-activation
structure, since every Bact/C*/P structure lacks U4 entirely, whereas *RNU2-2* and *RNU5B-1* are best resolved by activated and
catalytic states. Only the multi-state model exceeds 0.73 on all four dominant genes.

Against the degree-preserving null the observed 0.815 exceeded every one of 200 replicates (Figure 1c) (null mean 0.557, 95th
percentile 0.689; empirical p = 0.005) [22]. All three pre-specified criteria for this gate were met. Removing the
minor-spliceosome genes from training left the dominant result unchanged (0.810 → 0.813), excluding a
family-projection artefact [24].

### Held-out performance against thirteen baselines

**Table 1. MuST-VEP against 13 baselines under the pre-registered held-out-gene protocol.**

| Baseline | Dominant AUROC | Dominant ΔAUROC [95% CI] | Recessive AUROC | Recessive ΔAUROC [95% CI] |
|---|---|---|---|---|
| **MuST-VEP (this work)** | **0.810** | – | **0.662** | – |
| gnomAD constraint (1 − o/e, 10-nt window) | 0.766 | +0.044 [−0.123, +0.112] | 0.658 | +0.004 [−0.112, +0.111] |
| gnomAD allele density (±2 nt) | 0.686 | +0.124 [−0.137, +0.182] | 0.597 | +0.065 [−0.038, +0.172] |
| Trained secondary-structure model | 0.665 | +0.145 [+0.068, +0.213]* | 0.621 | +0.041 [−0.026, +0.116] |
| CADD (PHRED) | 0.655 | +0.121 [+0.021, +0.241]* | 0.625 | +0.037 [−0.083, +0.162] |
| Paralogue pathogenic density | 0.615 | +0.195 [+0.136, +0.310]* | 0.743 | −0.081 [−0.183, +0.013] |
| gnomAD allele density (position) | 0.609 | +0.201 [−0.064, +0.259] | 0.588 | +0.073 [−0.021, +0.169] |
| phyloP (447-way) | 0.593 | +0.217 [+0.079, +0.282]* | 0.639 | +0.023 [−0.035, +0.081] |
| Distance to nearest known pathogenic position | 0.554 | +0.256 [+0.185, +0.438]* | 0.669 | −0.007 [−0.109, +0.086] |
| Rfam element-type enrichment | 0.518 | +0.292 [+0.216, +0.363]* | 0.618 | +0.044 [−0.057, +0.146] |
| RNA-FM (foundation model, zero-shot) | 0.460 | +0.350 [+0.143, +0.408]* | 0.623 | +0.007 [−0.085, +0.095] |
| ViennaRNA pairing probability | 0.434 | +0.376 [+0.299, +0.534]* | 0.515 | +0.147 [+0.017, +0.284]* |
| ViennaRNA fold ΔΔG | 0.260 | +0.549 [+0.353, +0.646]* | 0.519 | +0.143 [+0.009, +0.277]* |
| ViennaRNA duplex ΔΔG | 0.251 | +0.598 [+0.194, +0.672]* | 0.636 | +0.056 [−0.077, +0.191] |

*Mean AUROC across held-out units. ΔAUROC is MuST-VEP minus the baseline, averaged over units, with a paired
position-block bootstrap 95% CI (2,000 resamples). \* marks intervals excluding zero. Hierarchical
(unit-then-position) intervals are wider throughout, reflecting that only four units contribute to each context.*

Two results in this table deserve comment. First, **RNA-FM performs below chance on dominant variants** (0.460 at
99.9% coverage). An RNA foundation model trained on sequence cannot rank pathogenic snRNA variants, which is
consistent with the guidance's position that general non-coding tools are not calibrated for this gene class, and
is worth recording because such models are increasingly offered as general-purpose non-coding predictors.

Second, **gnomAD constraint is not separable from MuST-VEP on clinical labels** (0.766 versus 0.810, CI −0.123 to
+0.112). This was a pre-specified gate and it failed. We return to it below, because the explanation turned out to
be the most consequential result in the study.

In a logistic model with a gene random intercept the odds ratio per standard deviation of score was 3.50
(2.74–4.48), against 1.32 (0.98–1.77) for CADD [24]. For recessive variants MuST-VEP did not beat conservation
(ΔAUROC versus phyloP +0.010, p = 0.74; versus CADD +0.023, p = 0.51), and no method exceeded 0.69 once healthy
homozygous carriers were used as controls. We report this as a limit of position-based prediction rather than a
property of any one method: recessive snRNA variants are hypomorphic, cluster in regions that tolerate heterozygous
variation, and are measurably milder in the one functional assay available.

### Ablation

**Table 2. Ablation of every component and of each conformational state.**

| Ablation | Dominant | Recessive | ΔDominant |
|---|---|---|---|
| **None (full MuST-VEP, 11 structures)** | **0.806** | **0.636** | – |
| Remove the contact-graph smoothing term | 0.788 | 0.635 | −0.018 |
| Single structure, chosen without the test gene | 0.597 | 0.654 | −0.209 |
| Single structure, chosen WITH the test gene (oracle) | 0.795 | 0.675 | −0.011 |
| Only 5Z56 (Bact-mature) | 0.795 | 0.644 | −0.011 |
| Only 6FF7 (Bact) | 0.792 | 0.634 | −0.015 |
| Only 5XJC (C*) | 0.752 | 0.608 | −0.054 |
| Only 7DVQ (minor-Bact) | 0.747 | 0.675 | −0.060 |
| Only 6QDV (P) | 0.687 | 0.599 | −0.120 |
| Only 6QW6 (tri-snRNP) | 0.655 | 0.637 | −0.152 |
| Only 6QX9 (pre-B) | 0.608 | 0.604 | −0.198 |
| Only 6Y5Q (U2-snRNP) | 0.587 | 0.625 | −0.219 |
| Only 8Y6O (minor-preB) | 0.582 | 0.654 | −0.224 |
| Only 3JCR (tri-snRNP) | 0.581 | 0.629 | −0.226 |
| Only 5O9Z (B) | 0.551 | 0.607 | −0.255 |
| ADD conformational-trajectory features [51] | n.s. | n.s. | +0.030, fails its degree-matched null (p = 1.00) |
| ADD protein-partner identity [52] | n.s. | n.s. | −0.026, does not transfer across genes |

The gap between the nested single structure (0.597) and the oracle single structure (0.795) is the whole argument
for the multi-state representation: a single state can be nearly as good as the ensemble, but only if you already
know which gene you are testing, which in deployment you do not.

**Two extensions were specified with kill criteria, tested and rejected.** Both are reported because a negative
result about an obvious extension is evidence about the problem, not an omission.

The first was a conformational-hinge hypothesis: that pathogenic variants sit at nucleotides whose protein partners
turn over between consecutive states, rather than at the most buried ones. On the saturation screen, trajectory
rewiring reached Spearman 0.150, but its partial correlation given burial was only +0.092 — below the +0.10 line
declared in advance — while partner count given rewiring retained +0.299. Decisively, under a degree-matched
permutation in which every position keeps its exact per-state partner count and only partner identity is shuffled,
the null produced a *stronger* correlation than the real data (+0.281 ± 0.031 against +0.150, p = 1.00). The
apparent +0.030 held-out gain resolved into a single gene with seven positives moving from 0.526 to 0.960 while two
other genes got worse [51].

The second was partner identity, motivated by large per-protein Mantel–Haenszel odds ratios (PPIL2 31.8, SF3B1 22.2,
PRPF8 7.7). Those are stratified by gene only, and pathogenic variants cluster *within* a gene, so any protein whose
footprint covers the cluster earns a large ratio for free. Requiring the information to transfer to a held-out gene
removes the confound, and it does not: held-out AUROC fell from 0.574 to 0.548, and identity alone reached 0.437,
below chance. Only nine partners appear at labelled positions in three or more of the five labelled genes, and eight
of those nine are the Sm ring, which binds every snRNA by construction [52]. We note that embedding partner proteins
with a protein language model is the load-bearing component of comparable geometric-deep-learning methods in
adjacent problems; it does not port here, and the reason is identifiable.

### Four evaluation regimes

**Table 3. Four evaluation regimes of increasing difficulty.**

| Regime | What | n | Positives | Performance |
|---|---|---|---|---|
| Held-out gene, dominant | Every variant of one gene withheld entirely | 1,211 | 72 | AUROC 0.766 [0.713, 0.811] |
| Held-out gene, recessive | Same, recessive positives | 566 | 177 | AUROC 0.591 [0.540, 0.639] |
| Held-out gene, pooled | Both mechanisms, rank-pooled | 1,777 | 249 | AUROC 0.636 [0.599, 0.672] |
| Zero-shot functional transfer | *RNU4-2* saturation genome editing; no functional measurement in training | 485 | continuous | Spearman 0.405 [0.255, 0.535] |
| Prospective, timestamped | Deposited 2026-10-09 with SHA-256, evaluation protocol fixed in advance | 2,777 | pending | to be evaluated on future reports |

Held-out-gene is strictly harder than a cluster split: no variant, position or paralogue of the test gene is seen in
training. In a temporal evaluation, trained only on evidence public before the ClinVar release of 4 May 2025 and
only on other genes, the system ranked the 41 variants that subsequently became P/LP above CADD (pooled ΔAUROC
+0.173, 0.039–0.304, p = 0.009; versus phyloP +0.111, p = 0.11) [17].

### Population constraint and molecular function diverge

This is the study's principal finding, it began as a failed gate, and it is summarised in Figure 2.

**Table 4. Zero-shot Spearman correlation with measured function (*RNU4-2* saturation genome editing, 485 variants).**

| Score | n | ρ | 95% CI |
|---|---|---|---|
| Graph-smoothed contacts | 485 | **0.437** | [0.293, 0.561] |
| MuST-VEP | 485 | **0.405** | [0.255, 0.535] |
| Protein-contact fraction | 485 | 0.347 | [0.198, 0.482] |
| gnomAD allele density (±2 nt) | 485 | 0.258 | [0.084, 0.410] |
| gnomAD constraint (1 − o/e) | 485 | 0.206 | [0.035, 0.344] |
| gnomAD allele density (position) | 485 | 0.149 | [−0.006, 0.287] |
| CADD (PHRED) | 411 | 0.072 | [−0.070, 0.208] |
| phyloP (447-way) | 485 | −0.004 | [−0.129, 0.128] |

Paired differences, all excluding zero: structure − constraint +0.231 (0.093–0.363, p = 0.002); structure − CADD
+0.334 (0.144–0.525); structure − phyloP +0.440 (0.271–0.604). Partial correlations are asymmetric in both
directions tested (Figure 2c): **graph-smoothed contacts retain ρ 0.394 controlling for constraint, whereas
constraint retains 0.025 controlling for contacts**; for the combined MuST-VEP score the pair is 0.358 against
0.038. Among the
121 most functionally damaging variants by assay, the structural score places 59% in its own top quartile, against
45% for density, 40% for constraint, 21% for CADD and 11% for phyloP. Within the dominant critical region
MuST-VEP reaches 0.575 (0.310–0.732).

The explanation is not subtle. Curated pathogenicity labels are assembled partly from population depletion, so a
depletion score predicts them partly for free. A benchmark built on such labels will rank a population-derived score
alongside a mechanistic one while the two carry different information about the molecule. This has been shown for
missense variants, where gene identity dominates ClinVar benchmarks and function-based evaluation reorders
predictors; we arrive second to that framing and say so. What we add is the non-coding case, where the effect is
larger, and where — unusually — both label types exist for the same variants.

The practical consequence is narrow and worth stating precisely: for non-coding variants, computational evidence
drawn from conservation or depletion is not independent of population-frequency evidence. Where a classification
framework treats them as separate lines, it risks counting one observation twice.

### Calibration against measured function

The 2026 guidance states that in silico tools should inform PP3/BP4 only once calibrated for snRNAs. Using the
*RNU4-2* saturation screen as functional truth (damaging = lowest quartile of function score; 121 damaging, 364
tolerated across 137 positions) and zero-shot scores with the U4 family withheld, MuST-VEP reaches a positive
likelihood ratio of 12.03 at its 90th percentile, with a position-block 95% CI lower bound of 4.87 — meeting the
ClinGen Bayesian threshold for PP3_moderate on the conservative bound — at 33% sensitivity. The graph-smoothed
contact feature alone reaches PP3_supporting at its 80th percentile (LR+ 7.41, lower bound 4.11, sensitivity 57%).
CADD and phyloP447 reach no evidence band at any threshold tested (LR+ 0.50–1.04).

This is calibration against a cell-fitness readout, not against clinical outcome, and should not be substituted for
clinical calibration. It is reported because it is the only functional truth set available for any snRNA gene, and
because the alternative — calibrating on curated labels — is precisely the circularity documented above.

### Allele-resolved scoring is a separate problem

The contact model scores positions, so all three substitutions at a nucleotide receive one score. We added
allele-resolved features from the same structures [44] and tested them where the answer is measurable: within the
137 *RNU4-2* positions at which the screen assayed all three alternative alleles. A position-level score has, by
construction, no within-position signal, so any correlation is attributable to allele identity. The allele score
reached a mean within-position Spearman of +0.346 with −SGE (Wilcoxon p = 0.0010), and pair disruption alone +0.374;
CADD, which is allele-resolved, reached −0.077. The most damaging allele was ranked first 44–46% of the time against
33% expected by chance, and 32% for CADD.

But adding those features to the clinical classifier made it worse: mean held-out AUROC for dominant variants fell
from 0.810 to 0.785 (−0.026), driven by *RNU2-2* (0.807 → 0.716), with no consistent gain for recessive variants
(+0.009) [46]. Ranking alleles within a position and discriminating pathogenic from benign across positions are different
tasks, and conflating them costs accuracy. We therefore report the allele score as a complementary output for the
specific clinical question of which substitution at a known position is worst, not as a component of the classifier.
The allele score varies only where the nucleotide is base-paired, which is 45 of 137 tested positions. This result
rests on a single dataset and should be treated as promising and unreplicated.

### Mechanism

Mantel–Haenszel odds ratios stratified by gene, at position level with BH-FDR [18a, 18b], placed dominant pathogenic
variants at contacts with the activation and catalytic machinery — PRPF8 OR 6.4 (FDR 5×10⁻⁴), RBM42 12.0 (0.003),
SART1 5.8 (0.003), SNRNP200/Brr2 4.4 (0.006), SNRNP27 4.9 (0.009) — with the Sm ring not enriched (0.8). Recessive
variants were enriched at snRNP-assembly contacts (SNU13 k-turn OR 4.4, FDR 0.036).

A pre-specified interaction test of this contrast did not meet its criterion: the joint Wald test over the three
feature × mechanism terms gave p = 0.022, but the contact term ran opposite to the predicted direction and *RNU2-2*
alone was flat (p = 0.895) [24]. We therefore report a gradient rather than a categorical separation. As shown in
Table 2, these per-protein associations do not survive the requirement to transfer across genes, so they describe
where pathogenic variants sit rather than providing independent predictive information.

### Part of this gene family cannot be assessed at all

*RNU1-3* and *RNU1-4* have zero callable alleles across 76,215 gnomAD genomes, while their near-identical siblings
*RNU1-1* and *RNU1-2* are at full cohort depth. In the HPRC Release 2 pangenome, built from long-read haplotype
assemblies, the same two genes carry 22 and 15 variants across 146 and 62 of 464 haplotypes — *RNU1-3* matching
*RNU1-2* exactly at 134 variants per kb [41]. The cause is not poor coverage but non-identifiability: the four RNU1
copies are identical across all 164 nucleotides, so no short read can ever be assigned among them. The same holds
for U6.

The same method distinguishes real absence from invisibility: of five genes returning zero variants of every class
in a published systematic screen, three (*RNVU1-1* and two U7 genes) are genuinely near-invariant in long-read
assemblies, and only *RNU1-3* and *RNU1-4* are not. *RNU4-2* itself shows 6 variants across 6 haplotypes, consistent
with genuine constraint at a confirmed disease gene. Across coding sequence a comparable artefact exists but is
small and already recognised: 89 of 34,556 MANE genes (0.26%) are missing at least half their expected synonymous
variation, and only 2 carry any ClinVar P/LP variant [42]. Coding sequence has an internal control — synonymous
variation — that non-coding sequence lacks.

U1 is the only major spliceosomal snRNA family with no established disease gene. We do not claim this is why. We
claim only that the question has not been asked of those loci, cannot be asked with short reads, and can be asked
with long reads.

### Limitations

Two of four dominant genes contribute 7 and 4 pathogenic variants, so the dominant result rests substantially on
*RNU4-2* and *RNU2-2*; the hierarchical intervals in Table 1 are correspondingly wide and several cross zero. The
function validation is one assay in one gene, because *RNU4-2* is the only snRNA with a saturation screen. The
allele result is unreplicated. The pangenome comparison uses 232 samples and graph-based genotypes that can err in
repetitive sequence. GENCODE v27 models *RNU4-2* as 141 nt against the current 145, leaving four 3′ nucleotides
unscored. The model scores positions for indels, which matters because the single most common ReNU allele is an
insertion. Recessive prediction does not work, for any method tested.

---

## Conclusion

Structural context across the splicing cycle is, among the evidence types tested here, the only one that tracks what
variants do to the molecule. Population-derived scores track what selection has done to the position, and on curated
labels the two are hard to tell apart — which is why a pre-specified gate comparing them failed, and why that
failure, once resolved against an ascertainment-free readout, became the study's main result. For a gene class where
clinical interpretation is explicitly uncalibrated, that distinction decides what computational evidence is worth.

MuST-VEP provides scores for 13,280 variants across 17 spliceosomal snRNA genes, with within-gene percentiles,
features, and an annotation of which loci population sequencing can and cannot assess. A set of 2,777 predictions
was deposited with a cryptographic hash and a fixed evaluation protocol before any of the variants concerned were
reported, and will be evaluated as new cases accumulate.

---

## Figures

**Figure 1.** The MuST-VEP workflow and the case for a multi-state representation. (a) Pipeline: eleven cryo-EM
states and 17 GENCODE snRNA genes feed contact extraction and family projection; multi-state aggregation with
one-hop graph smoothing produces the features; a mechanism-aware L2 logistic model with nested leave-one-gene-out
scorer selection produces position scores, allele scores where the nucleotide is paired, and a ClinGen Bayesian
calibration. (b) Nucleotides resolved per snRNA family in each state. U4 is absent from every structure after
activation, which is why no single state can serve *RNU4-2* and *RNU2-2* simultaneously. (c) Gate G1: mean
held-out AUROC for dominant variants, each single state, the nested single-state selection, and the
degree-preserving shuffled-contact null (200 replicates). The full ensemble reaches 0.815 and exceeds every
replicate.

**Figure 2.** Population constraint and molecular function diverge. (a) The same comparison under two readouts,
as paired differences rather than on a shared axis, because AUROC and Spearman are not commensurable: on
clinical labels MuST-VEP − constraint is +0.044 [−0.123, +0.112] and covers zero; on measured function it is
+0.199 [+0.057, +0.341] and does not. (b) Zero-shot Spearman correlation with saturation-editing damage for
every score, with position-block 95% intervals; structural features in blue, population-derived in orange.
(c) Partial correlations, showing the information is not shared: graph-smoothed contacts retain 0.394
controlling for constraint, while constraint retains 0.025 controlling for contacts. (d) The underlying scatter
for MuST-VEP and CADD on the same 485 variants, with the U4 family withheld from training.

---

## Key points

- Spliceosomal snRNA variants are a frequent cause of neurodevelopmental disorder, but current guidance states that
  existing in silico tools are not calibrated for non-coding RNA genes and notes that a single structure cannot
  capture how snRNA contacts change through the splicing cycle.
- MuST-VEP represents each nucleotide by its contacts across eleven cryo-EM structures spanning the cycle. The
  multi-state representation beats any single state chosen without access to the test gene (0.815 versus 0.631) and
  exceeds all 200 degree-preserving shuffled-contact nulls.
- Against thirteen baselines under a held-out-gene protocol, MuST-VEP reaches 0.810 for dominant variants; the RNA
  foundation model RNA-FM reaches 0.460, below chance.
- Population constraint matches the model on curated clinical labels but not on measured function (partial
  correlations 0.025 versus 0.394 on saturation genome editing). Curated labels and molecular function are not
  interchangeable benchmarks for non-coding variants, because the labels are partly derived from population data.
- Two obvious extensions — conformational-trajectory features and protein-partner identity — were specified with
  kill criteria in advance, tested, and rejected; both are reported rather than omitted.

---

## Author actions outstanding

1. Deposit the prospective prediction set (`results/prospective_deposit_2026-10-09.tsv`, SHA-256 `b5b2249a…`) on
   Zenodo and cite the DOI.
2. Generate the figures listed in `FIGURE_PLAN.md`, in particular **Figure 1**, the workflow schematic — both
   comparable papers in this venue lead with one and it does substantial work.
3. Send outreach email #1 in `OUTREACH_DRAFTS.md` to the saturation-editing group; time-sensitive.
4. Optional: AlphaGenome as a fourteenth baseline, which needs a free DeepMind API key.
