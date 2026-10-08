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

Calibrated against measured function, the score reaches a positive likelihood ratio of 12.0 at its 90th
percentile (position-block 95% CI lower bound 4.87), meeting the ClinGen Bayesian threshold for PP3_moderate on
the conservative bound, at 33% sensitivity; CADD and phyloP447 reach no evidence band at any threshold (LR+
0.50-1.04) [sge_calibration]. This is a functional rather than clinical calibration and is labelled as such.

**Conclusions.** Multi-state structural context identifies dominant snRNA variants that conservation and genome-wide
predictors miss, and is the only evidence type tested here that tracks measured molecular function. We provide
precomputed scores for 13,280 variants across 17 spliceosomal snRNA genes, annotated with which loci population
sequencing can and cannot assess.

---

## Changes from v0.1 and why

- **Title no longer claims to "separate dominant from recessive."** The pre-specified interaction test (stacked
  logistic model, position-clustered SEs, joint Wald) gave p = 0.022 but with the contact×recessive term in the wrong
  direction, and RNU2-2 alone was flat (p = 0.895). Criterion not met [24]. Reported as a gradient.
- **The guidance claim is corrected.** v0.1 said guidance states "no in silico tool captures features relevant to
  snRNA variant deleteriousness." The source actually says tools were not "designed or calibrated specifically for
  variants in non-coding RNA genes." Verified against the full text.
- **The constraint/function dissociation is new and now carries the paper.** It did not exist in v0.1; the G2 gate
  failure was unexplained there.
- **The callability result is new** [e372a5f, 09fe17c].
- **BP4_supporting is dropped** — it did not survive the hierarchical bootstrap.
- **Label set is v5c throughout** (v5 with the 16 UKB-homozygote conflicts removed); v5 retained as sensitivity.
- **RNA-FM negative result moved to the main text body**, not the abstract — it is a secondary point.

## Still to write

- Methods: allele-specific scoring and how indels are scored (pre-specified, not yet run)
- Novelty table vs the SGE, RNU2-2 and RNU4ATAC papers
- Figure plan
- Data/code availability; Zenodo deposit of the timestamped prediction set
