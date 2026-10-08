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
