# Introduction and Discussion (draft)

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
