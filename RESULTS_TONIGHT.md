# snRNA-VEP: results (night of 7 Oct 2026; not committed until 8 Oct)

## Headline so far
1. **Cryo-EM contact profiles predict wet-lab variant damage in RNU4-2 without any RNU4-2 training data;
   CADD and conservation do not.** Spearman with −SGE function score (501 variants, De Jonghe et al. Nature 2026):
   fraction of splicing states with protein contact **0.35**, max protein contacts **0.32**, snRNA–snRNA contacts 0.23,
   CADD 0.08 (n.s.), phyloP447 −0.01 (n.s.), ViennaRNA folding ≈0.
2. **Held-out-gene pathogenicity (any mode), SNVs that CADD scores (907 variants, 63 pathogenic):**
   max protein contacts − CADD = **+0.18 AUROC [0.08, 0.28], p<0.001**; logistic model (contacts only) +0.16
   [0.07, 0.26]; logistic full model +0.13 [0.01, 0.25], p=0.02.
3. **RNU4-2 zero-shot (ClinVar P/LP vs gnomAD):** protein contacts 0.73 vs CADD 0.63 vs phyloP447 0.55.
4. **Honest caveat:** the current trained models (logistic / gradient boosting on all features) do not beat the single
   best contact feature. Few labels per gene; extra features add noise. Next: the multi-state spliceosome network
   trained on dense SGE labels plus patient labels.

## Labels (scripts/05b_curate_supplements.py; 9 papers, every row traceable)
| Gene | Dominant NDD | Recessive NDD | RP |
|---|---|---|---|
| RNU4-2 | 20 | 9 | 2 |
| RNU2-2 | 12 | 38 | – |
| RNU5B-1 | 5 | – | – |
| RNU6-1/2/8/9 | – | – | 4 |
Controls: 353 UK Biobank/All of Us variants (RNU4-2, gnomAD-independent), 14 UKB biallelic carriers,
50 RP-benign variants, gnomAD v4.1 PASS variants for other genes. 524 RNU4-2 SGE scores.

## Files
- results/benchmark_main_dominant.tsv, benchmark_main_anypath.tsv, benchmark_main_sge_zeroshot.tsv
- results/zero_shot_RNU4-2_clinvar.tsv, depletion_*.tsv, label_free_*.tsv

## Model development log (all on identical held-out test sets; nothing tuned on test genes)
1. **Multi-state graph network v0 (per-state node features, 2-layer message passing, multi-task):** inverted on held-out
   genes (pooled AUROC 0.37; SGE ρ −0.34). Diagnosis: (a) population variants in non-disease paralog copies
   (RNU4-1, RNU5D/E/F-1, …) sit at the same structural positions as patient variants in the expressed copy, so
   labelling them benign teaches "structurally critical → benign" → training restricted to disease genes;
   (b) per-state and partner-identity node features are family-specific (U2 nucleotides appear in other complexes than
   U4), so rules learned on U2/U5/U6 do not transfer to U4 → replaced by family-agnostic node features
   (state-averaged contact summaries + within-family percentile ranks). After both fixes, v0 is above chance on every
   held-out gene (0.58–0.78) but not better than single contact features.
2. **snRNA-VEP v1 = Simplified Graph Convolution (k=2) + L2 logistic regression** (low capacity for ~90 positives):

| held-out unit | RNU2-2 | RNU4-2 | RNU5B-1 | RNU6 (RP) | mean |
|---|---|---|---|---|---|
| snRNA-VEP v1 | 0.674 | 0.628 | 0.960 | 0.953 | **0.804** |
| same, no graph propagation | 0.671 | 0.676 | 0.844 | 0.949 | 0.785 |
| frac. states with protein contact | 0.633 | 0.698 | 0.481 | 0.933 | 0.686 |
| CADD | 0.525 | 0.601 | 0.915 | – | 0.680 |
| phyloP447 | 0.609 | 0.467 | 0.739 | 0.491 | 0.576 |

   Rank-pooled paired tests: v1 − CADD (SNVs) **+0.10 [0.00, 0.20], p=0.047**; v1 − phyloP **+0.12 [0.03, 0.21],
   p=0.003**; v1 − contacts alone +0.03 (n.s.); graph propagation vs none −0.01 (n.s.).
   SGE zero-shot (RNU4-2): v1 ρ 0.16, raw contact fraction **0.35**, CADD 0.07, phyloP 0.00.
   Note: gene-out and family-out coincide here because each family has one disease gene (U6 genes form one unit).

## Labels v2 (scripts/05e_labels_v2.py) — supersedes the table above for everything below
Curated supplements + ClinVar 2026-10-04 P/LP (≥1 star). **Mode-aware controls**: for recessive genes (RNU2-2,
RNU4ATAC, RNU12) only gnomAD homozygotes (nhomalt ≥ 1) + ClinVar B/LB count as controls, because heterozygous carriage
says nothing about recessive benignity; RNU4-2 uses UKB/AoU (gnomAD-independent); dominant-only genes use gnomAD AC ≥ 1;
gnomAD restricted to the transcribed region.

| gene | pathogenic | controls |
|---|---|---|
| RNU2-2 | 64 (dominant 26, recessive 38) | 24 |
| RNU4-2 | 39 (dominant 30, recessive 9) | 345 |
| RNU4ATAC | 27 (recessive) | 13 |
| RNU12 | 7 (recessive) | 8 |
| RNU5B-1 | 5 (dominant) | 223 |
| RNU6-1/2/8/9 | 4 (dominant RP) | gnomAD |

With these stricter labels the v1 SGC dropped to mean 0.667 (below raw contacts 0.739, CADD 0.726): v1's lead was partly
helped by easy heterozygous-carrier controls in recessive genes. Reported, not hidden.

## Mechanism finding (drives the v2 design)
Dominant (gain/altered-function) snRNA variants are **seen by contact profiles and missed by conservation**; recessive
(loss-of-function) variants are **seen by conservation**. Untrained mean of the two ranks: dominant 0.741, recessive 0.772.
This mirrors ACMG logic (evidence depends on the disease mechanism) and is, to our knowledge, not described for snRNAs.

## snRNA-VEP v2 (scripts/11_v2_model.py): mechanism-aware evidence model
Unsupervised within-gene ranks of 5 features (protein-contact fraction across 11 cryo-EM states, max protein contacts,
snRNA–snRNA contacts, **graph-smoothed contact fraction over 3D neighbours in the multi-state contact graph**, phyloP447),
interacted with the mechanism context (dominant / recessive); L2 logistic regression (11 coefficients).
Held-out unit (the whole gene is absent from training), AUROC:

| dominant context | RNU2-2 | RNU4-2 | RNU5B-1 | RNU6 (RP) | **mean** |
|---|---|---|---|---|---|
| **snRNA-VEP v2** | 0.802 | **0.819** | 0.750 | 0.926 | **0.824** |
| contacts + phyloP (untrained) | 0.935 | 0.642 | 0.613 | 0.802 | 0.748 |
| protein-contact fraction | 0.786 | 0.716 | 0.486 | 0.930 | 0.730 |
| CADD | 0.674 | 0.586 | 0.915 | – | 0.725 |
| phyloP447 | 0.860 | 0.477 | 0.741 | 0.491 | 0.642 |

| recessive context | RNU2-2 | RNU4-2 | RNU4ATAC | RNU12 | mean |
|---|---|---|---|---|---|
| snRNA-VEP v2 | 0.736 | 0.533 | 0.642 | 0.878 | 0.697 |
| contacts + phyloP (untrained) | 0.888 | 0.508 | 0.844 | 0.837 | **0.769** |
| phyloP447 | 0.853 | 0.361 | 0.858 | 0.704 | 0.694 |
| CADD | 0.692 | 0.555 | 0.865 | 0.571 | 0.671 |

Rank-pooled paired bootstrap over all held-out units and contexts: v2 − phyloP **+0.075, p=0.025**; v2 − CADD (SNVs)
+0.066, p=0.096. Final-model coefficients: graph-smoothed contact fraction **+0.85 (largest)**, snRNA contacts +0.56,
recessive context +1.34; phyloP is used mainly in the recessive context. The graph contributes the single strongest feature.

**Nested selection (v2-nested):** for each held-out unit, the scorer (v2 / untrained combo / contacts / phyloP) is chosen
by an inner leave-one-unit-out loop over the training units only. Dominant: v2 chosen in every unit, mean **0.824**.
Recessive: mean 0.701 (RNU4-2 0.508, RNU2-2 0.752, RNU4ATAC 0.708, RNU12 0.837).
**Honest status:** v2 is the best dominant-mechanism predictor; for recessive variants it does not beat simple
conservation-based scores yet (RNU4-2 recessive, 9 variants, is hard for every method).

## Time split (scripts/12_timesplit_sge.py): could v2 have predicted the later discoveries?
Trained only on RNU4-2 (ReNU, first reported 2024-04; 28 dominant variants + 322 controls), dominant context, then
applied to genes discovered later:

| later-discovered gene | v2 (RNU4-2 only) | contacts | phyloP447 | CADD |
|---|---|---|---|---|
| RNU2-2 (2024-09) | 0.763 [0.61, 0.89] | 0.786 | **0.860** | 0.674 |
| RNU5B-1 (2024-10) | 0.645 [0.51, 0.78] | 0.486 | 0.741 | **0.915** |
| RNU6 RP genes (2025-01) | **0.915** [0.86, 0.96] | 0.930 | 0.491 | – |

Honest reading: training on one gene transfers (all above chance), but no single method wins every later gene with
one-gene training; the leave-one-unit-out v2 (trained on all other genes) is better and is the main result.

## SGE zero-shot, v2 (no RNU4-2 labels in training)
Spearman with −SGE function score (485 variants with a structure node): v2 trained without RNU4-2 **0.412**; without the
whole U4 family 0.413; graph-smoothed contact fraction alone **0.437**; raw contact fraction 0.347; CADD 0.072;
phyloP447 −0.004. The multi-state graph smoothing improves wet-lab correlation over raw contacts (0.44 vs 0.35).

## Clinical calibration (scripts/14_calibration.py; held-out predictions only, dominant context)
ACMG/ClinGen Bayesian thresholds (OddsPath 2.08 / 4.33 / 18.7; benign 0.48 / 0.23), requiring the 95% bootstrap bound to
meet the threshold:
- **BP4_supporting reached** by v2 < 0.238 (LR− 0.30, bound 0.47; 41% of controls, 12% of pathogenic variants).
  CADD reaches BP4_supporting only for 5% of controls; phyloP and raw contacts reach nothing.
- **PP3 not reached when units are pooled**: held-out models from different folds put scores on different scales.
  Within units the point LR+ is high (RNU2-2 9.5 at 50% sensitivity; RNU6 10.3; RNU4-2 21 at very low sensitivity), but
  with 4–30 positives per gene the bounds are too wide. **Honest status: v2 is a prioritisation score; PP3-grade
  calibration needs more labelled variants (e.g. the ongoing RNU2-2 / RNU5B-1 SGE screens) and is future work.**

## More data (labels v3 / v4; scripts/05f, 05g, 05e --v3/--v4/--v4s; snrna_vep/literature.py)
New sources (all GRCh38 reference-checked):
- Genet Med 2026 RNU4ATAC cohort (Matalon et al.) Supp. Table 2: 105 patient alleles, ACMG P/LP kept (VUS dropped).
- AJHG 2026 bi-allelic RNU6ATAC/RNU4ATAC (Tables 1–2) and HGG Adv 2026 RNU6ATAC (Table 1): **new gene RNU6ATAC** (9 variants).
- iScience 2026 Table S1: 452 patient-variant reports compiled from 65 publications (separate "compiled" tier; no per-variant
  classification, so v4s additionally requires ≥2 independent publications).
- Deliberately NOT used: Nat Genet 2026 RNU2-2 ST5 (prioritised candidate genotypes in unsolved cases, 4 alleles are
  gnomAD homozygotes); LOVD (blocks automated access; respected).
**Bug found and fixed:** literature numbers RNU4ATAC on NR_023343.1, which starts 1 nt upstream of the GENCODE model.
Literature n. positions are shifted accordingly; validated on 7 classic variants that reproduce their ClinVar P/LP records.
(ClinVar-derived labels were never affected: they use genomic coordinates.)
**Evaluation rule added:** a held-out unit needs ≥3 pathogenic AND ≥3 controls (RNU6ATAC has 1 control with a structure
node, so it is used for training only; an earlier AUROC 1.0 on it was meaningless and is not reported).

| label set | pathogenic (distinct) | controls |
|---|---|---|
| v2 | 119 (with structure node) | 1,141 |
| v3 | 146 | 4,448 rows |
| v4 | 266 | 4,418 rows |
| v4s (compiled ≥2 publications) | 161 | 4,441 rows |

### v2 + nested selection on labels v4 (held-out units; 243 pathogenic with a structure node)
| context | snRNA-VEP (nested) | v2 | contacts+phyloP | contacts | CADD | phyloP447 |
|---|---|---|---|---|---|---|
| dominant (4 genes) | **0.807** (v2 chosen in 4/4) | 0.807 | 0.732 | 0.716 | 0.734 | 0.635 |
| recessive (4 genes) | **0.781** (combo chosen in 4/4) | 0.730 | 0.781 | 0.732 | 0.696 | 0.717 |

With twice the labels the structure model stays clearly best for dominant variants (+0.07 over CADD, +0.17 over phyloP), and
for recessive variants the nested procedure consistently selects the untrained contact+conservation combination.
The final system is therefore **mechanism-aware**: a trained multi-state structure model for dominant variants and
structure+conservation evidence for recessive variants, with the choice made by nested cross-validation, not by us.

**Sensitivity (labels v4s: compiled variants need ≥2 independent publications; 161 pathogenic):** dominant nested
**0.826** (v2 chosen 4/4; CADD 0.731, phyloP 0.643) — robust. Recessive nested 0.722: with fewer recessive labels the
inner selection alternates between contacts alone and contacts+phyloP (the latter scores 0.785 if always used), so the
recessive choice is label-volume dependent; reported as such.

## v3 (RNA language model + per-mechanism models; scripts/13_v3_model.py, labels v4) — negative result, reported
- RNA-FM zero-shot (masked marginal for SNVs, Δlog-likelihood otherwise; weights verified: 93% masked accuracy on U6):
  dominant 0.471 (inverted on RNU6: 0.085), recessive 0.724 (≈ phyloP), SGE ρ 0.01–0.15. **A sequence language model
  captures conservation-like signal but not the structural mechanism of dominant snRNA variants.**
- v3 vs v2: dominant 0.814 vs 0.807 (n.s.), recessive worse (0.657); gene-balanced v3 recessive 0.749 (< combo 0.781).
- Adding v3 candidates to nested selection lowered dominant to 0.776 (more candidates, noisier selection).
- **Final system stays v2-nested with its original 4 candidates** (specified before v3 existed).

## Prospective ClinVar test (scripts/17_prospective.py): time split AND gene split
Train only on knowledge public before the ClinVar 2025-05-04 release (54 pathogenic), only on other genes; test on the
41 variants that became P/LP afterwards (ClinVar 2026-10-04) vs the gene's controls.

| new P/LP after 2025-05 | v2 | contacts+phyloP | contacts | phyloP447 | CADD |
|---|---|---|---|---|---|
| RNU4-2 dominant (8 vs 298) | **0.716** | 0.570 | 0.731 | 0.382 | 0.522 |
| RNU2-2 dominant* (20 vs 15) | 0.800 | 0.907 | 0.782 | 0.810 | 0.665 |
| RNU4ATAC recessive (7 vs 12) | 0.821 | 0.857 | 0.738 | 0.881 | 0.875 |
| RNU12 recessive (6 vs 7) | 0.714 | 0.833 | 0.917 | 0.679 | 0.514 |

The final system (v2 for dominant, contacts+phyloP for recessive) gives dominant mean 0.758 (CADD 0.594, phyloP 0.596) and
recessive mean 0.845 (CADD 0.695, phyloP 0.780). Rank-pooled v2 − CADD +0.111 [−0.018, +0.235], p=0.09: consistent
direction, underpowered (41 new variants). *Mode from ClinVar condition text; some RNU2-2 entries may be recessive.

## Data sources tried and not usable (documented for transparency)
LOVD (blocks automated access), UK Biobank Allele Frequency Browser (robots.txt disallows /api; Cloudflare challenge),
All of Us data browser (counts rounded to 20, rare homozygotes hidden), jMorp (data-use agreement), medRxiv PDFs (403).
Manual option for the authors: query UKB AFB in a browser for homozygote counts in RNU2-2, RNU4ATAC, RNU12, RNU6ATAC
to enlarge recessive controls.

## Labels v5: UK Biobank WGS homozygote controls (scripts/05h; exports made by a person from afb.ukbiobank.ac.uk)
~490k genomes. Reference-checked; QC: allele number ≥ 50% of the gene median (1 RNU4ATAC site removed). Homozygous
variants in RNU2-2 80, RNU4ATAC 5, RNU12 16, RNU6ATAC 3. Recessive controls with a structure node: RNU2-2 19→71,
RNU12 7→16, RNU4ATAC 12→14, RNU6ATAC 1→2 (still training-only).
- **Constraint:** no UKB homozygotes in RNU2-2 n.1–69 (dominant hotspot half) or RNU12 n.1–41 (branch-point pairing end).
- **16 conflicts:** variants labelled pathogenic (recessive) carried homozygously by ≥1 UKB adult: RNU2-2 15 (11 from the
  compiled tier, 4 curated/ClinVar), RNU4ATAC 1 (n.8C>T, LP in the GIM 2026 cohort; AF 0.06%). Likely hypomorphic alleles
  or label errors; kept pathogenic in v5, removed in v5c.

### v2 / nested on labels v5 (held-out units)
| | v2 | contacts+phyloP | contacts | CADD | phyloP447 |
|---|---|---|---|---|---|
| dominant (4 genes) | **0.804** | 0.700 | 0.688 | 0.659 | 0.600 |
| recessive (4 genes) | 0.632 | **0.690** | 0.655 | 0.627 | 0.637 |
Nested: dominant 0.804 (v2 chosen 4/4); recessive 0.655.
**Stringent, realistic recessive controls (rare variants seen homozygous in healthy adults, in the same 3′ regions as
recessive pathogenic variants) make recessive prediction hard for every method (0.63–0.69); the dominant result is
robust and its margin grows** (RNU2-2 dominant: v2 0.784 vs CADD 0.473, phyloP 0.702).

## Mechanism (scripts/18a, 18b; labels v5): which proteins and states
Mantel–Haenszel odds ratios (strata = genes), position level, BH-FDR:
- **Dominant** variants are enriched at contacts with the activation/catalytic machinery: PRPF8 OR 6.4 (FDR 5e-4),
  SART1 5.8 (0.003), **RBM42 12.0 (0.003; independently recovers the 'RBM42-interacting region' of severe ReNU, Nava et
  al. 2025)**, SNRNP200/Brr2 4.4 (0.006), SNRNP27 4.9 (0.009), SF3B1 13.6 (0.006), PPIL2 17.1; Sm ring not enriched (0.8).
- **Recessive** variants are enriched at snRNP-assembly contacts: SNU13 (k-turn) OR 4.4 (FDR 0.036), SF3A 5.5 (0.07),
  Sm ring 1.7 (trend).
- States: protein contact in the **B complex** (activation) separates dominant pathogenic from controls best for RNU4-2
  (0.64) and RNU6 (0.84).
**Interpretation: dominant snRNA disease = disrupted contacts of the activation/catalytic machinery; recessive = disrupted
snRNP biogenesis contacts.** Derived systematically from 11 cryo-EM states, consistent with independent literature.

## FINAL (labels v5) — paired tests, prospective, time split, SGE (scripts/16, 17, 12)
Final system = v2-nested (v2 for dominant, chosen 4/4; contacts+phyloP / contacts for recessive). Rank-pooled within unit.

| comparison (held-out genes) | v5 ΔAUROC [95% CI], p | v5c (conflicts removed) |
|---|---|---|
| dominant: final − CADD | **+0.192 [0.096, 0.283], p<0.001** | +0.205 [0.106, 0.299], p<0.001 |
| dominant: final − phyloP447 | **+0.181 [0.100, 0.267], p<0.001** | +0.198, p<0.001 |
| dominant: final − contacts alone | **+0.100 [0.048, 0.154], p=0.001** | +0.106, p<0.001 |
| recessive: final − CADD / phyloP | +0.014 (p=0.69) / +0.003 (p=0.95) | +0.023 / +0.010 (n.s.) |
| both: final − CADD | +0.064 [0.007, 0.116], p=0.021 | +0.074, p=0.013 |
| both: final − phyloP447 | +0.057 [0.012, 0.100], p=0.015 | +0.068, p=0.005 |

**Prospective ClinVar (train on knowledge before 2025-05-04, other genes only; 41 later P/LP):**
dominant v2 RNU2-2 0.797, RNU4-2 0.716 (mean 0.757) vs CADD 0.458/0.522, phyloP 0.670/0.382; recessive (system uses
contacts+phyloP) RNU4ATAC 0.806, RNU12 0.672 vs CADD 0.800/0.512, phyloP 0.832/0.536. Pooled v2 − CADD **+0.173
[0.039, 0.304], p=0.009**; v2 − phyloP +0.111 (p=0.11).
**Time split (train on RNU4-2 only):** RNU2-2 0.719 (phyloP 0.702, CADD 0.473), RNU5B-1 0.672 (CADD 0.920), RNU6 0.928.
**SGE zero-shot (no RNU4-2 labels):** v2 ρ 0.402 (no U4 family 0.404); graph-smoothed contacts 0.437; contacts 0.347;
CADD 0.072; phyloP −0.004.

## Resource (scripts/19): results/snrnavep_scores.tsv.gz
13,280 variants in 17 spliceosomal snRNA genes (disease genes + paralog/candidate genes RNU5A-1, RNU4-1, RNU6-7,
RNU11, RNU5D/E/F-1): dominant and recessive scores, within-gene percentiles, features, CADD, training label.
Candidate gene note: RNU5A-1 de novo patient variants (Nava et al. 2025; VUS) rank around the 59–61st within-gene
percentile (loop I), n.62G>A 51st, n.73C>T 28th; RNU5B-1 pathogenic median 74th. Weak signal; descriptive only.

## Limitations to state
GENCODE v27 RNU4-2 model is 141 nt (current annotation 145 nt): 3′-terminal 4 nt not scored (no known patient
variants). RNU5B-1 (7) and RNU6 (4) have few positives. Recessive prediction from position is limited with stringent
controls. Not calibrated to PP3 (BP4_supporting only). Literature-compiled tier lacks per-variant classification
(v4s/v5c sensitivity analyses).
