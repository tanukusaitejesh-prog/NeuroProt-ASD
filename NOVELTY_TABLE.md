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
