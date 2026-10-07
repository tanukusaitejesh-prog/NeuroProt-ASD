# Pre-specified analysis plan (frozen before any new analysis)

Written 2026-10-07T19:11Z (00:41 IST 8 Oct), committed before scripts 22+ were written or run. The git commit
timestamp of this file is the registration time. Anything not specified here and run later is labelled
**exploratory** in RESULTS and in the manuscript.

## 1. Frozen elements (no changes after this commit)

| element | frozen value |
|---|---|
| primary label set | `labels_v5c.tsv` (v5 with the 16 UKB-homozygote conflicts removed). `labels_v5.tsv` = sensitivity analysis only |
| structural features | within-gene percentile ranks of: fraction of splicing states with protein contact; max protein residues; fraction of states with snRNA contact; 1-hop graph-smoothed contact fraction (row-normalised adjacency incl. self loop, RNA-RNA contact graph); phyloP447 |
| model (v2) | StandardScaler + L2 logistic regression, C = 0.3, class_weight = balanced; design = [F, F x recessive, recessive] (11 coefficients) |
| final system | v2-nested: scorer chosen per held-out unit and context from {v2, contacts+phyloP mean rank, contacts, phyloP} by inner leave-one-unit-out on training units |
| structures | 11 cryo-EM structures / 10 states as in `data/processed/contacts_per_structure.tsv` |
| evaluation | held-out unit (gene; RNU6-1/2/8/9 as one unit); a unit/context is evaluated only with >= 3 pathogenic and >= 3 controls; primary metric = mean per-unit AUROC, dominant context (4 units: RNU4-2, RNU2-2, RNU5B-1, RNU6) |

No hyperparameter, feature, or label change is allowed after this point. A failing result is reported, not repaired.

## 2. Gate experiments (decide the paper's framing)

### G1: multi-state vs single-state (script 22)
Comparisons use identical features recomputed from the selected structures (script 20 definitions; ALL must
reproduce the main features).
- **Single-state** = all structures of one state (tri-snRNP = 3JCR + 6QW6; otherwise one structure per state).
- **NESTED single state** (primary comparator): for each held-out unit, the state is chosen by mean AUROC over the
  other units only.
- **ORACLE single state** (descriptive, optimistic): the best state chosen on the test units.
- **Leave-one-state-out (LOSO)**: ALL minus each state.
- **Degree-preserving shuffled-contact null** (200 replicates): within each structure and paralog family, the
  per-position contact profiles (protein residues, snRNA contacts) are permuted among resolved positions, which
  keeps each structure's degree sequence. The RNA-RNA contact graph is rewired by double-edge swaps (10 x |E|
  swaps per structure), which keeps every node's degree. The full model is refitted on every replicate.

**G1 passes** if all three hold (dominant context, labels v5c):
1. mean AUROC(ALL) > mean AUROC(NESTED single state), and the 95% position-block bootstrap CI of the mean per-unit
   difference excludes 0;
2. ALL > NESTED single state in >= 3 of 4 dominant units;
3. ALL exceeds the 95th percentile of the shuffled-contact null (empirical one-sided p < 0.05).
ALL vs ORACLE: non-inferiority with margin 0.02 is reported but is not part of the gate (the oracle uses test labels).

### G2: stronger baselines (script 23), identical held-out protocol
1. ViennaRNA base-pair disruption (ddG of folding, change in pairing, duplex ddG; existing features) and a
   trained secondary-structure model (v2 design, structural features replaced).
2. Annotated structural elements: element membership (Sm site, k-turn, stems, loops, single-stranded functional
   regions) from Rfam / literature secondary structures; score = pathogenic rate of the element class in training
   units (leave-unit-out), fallback 0.5.
3. Published RNU4ATAC RNA-structure score, if it can be obtained; otherwise reported as not available.
4. Per-position population variant density: number of OTHER gnomAD v4.1 PASS alleles at the same position
   (the variant's own allele excluded) and in a +/-2 nt window. Circularity caveat: for RNU5B-1 and RNU6 the
   controls themselves are gnomAD variants, so this baseline is advantaged there.
5. MitoTIP-style paralog density: number of training-unit pathogenic variants at the aligned family position
   (+/-1 nt), leave-unit-out.
6. Distance to the critical region: minus the distance (family coordinates) to the nearest training-unit
   pathogenic position, leave-unit-out.

**G2 fails** if population density (4) reaches mean dominant AUROC >= final system − 0.02, or the
position-block CI of (final − density) includes 0. Other baselines matching the final system (difference <= 0.02)
are reported and the claims are narrowed accordingly.

**Decision rule:** if G1 or G2 fails, the paper is reframed as a resource and benchmark (scores plus negative
results), and the target tier is lowered.

## 3. Rigor analyses (script 24)
- **Uncertainty**: hierarchical bootstrap (resample units, then positions within units, 2,000 replicates) and the
  position-block bootstrap, for the mean per-unit ΔAUROC; per-unit AUROCs with position-block 95% CIs. A logistic
  GLMM with a gene random intercept (score effect) is a secondary analysis. **No p-values in the abstract.**
- **Allele-specific features** (script 25): per state, whether the position is base-paired (geometric detection in
  the cryo-EM models) and whether the alternative allele keeps a canonical or wobble pair; transition vs
  transversion; base vs backbone protein contact. Validation = within-position allele ranking in the RNU4-2 SGE
  (Spearman within positions that have all 3 alleles measured, averaged; and the fraction of positions where the
  most damaging allele is ranked first). Indels: scored at the affected node (the first deleted or the anchor+1 base),
  with allele features set to "pair broken" for deletions of paired bases. This is stated in Methods.
- **Dominant vs recessive interaction** within RNU4-2 and RNU2-2: stacked logistic model
  y ~ (contact + graph-contact + phyloP) x context + gene, rows = dominant cases + controls and recessive cases +
  controls; position-clustered SEs; joint Wald test of the interaction terms. The mechanism claim is supported if
  the joint test gives p < 0.05 with contact x recessive < 0 and phyloP x recessive > 0.
- **Family-projection confound**: rerun held-out evaluation with the minor-spliceosome genes (RNU4ATAC, RNU12,
  RNU6ATAC) excluded from training.
- **Gradient**: Spearman of scores with the continuous RNU4-2 SGE function score, plus the SGE functional classes;
  severity / splicing-extent data are used if public.

## 4. Novelty additions (computational parts)
- Timestamped predictions for every possible SNV in RNU2-2 and RNU5B-1 (SHA-256 recorded and committed; Zenodo
  deposit by the authors).
- Calibration against RNU4-2 SGE classes (proof of concept). BP4_supporting is dropped unless its lower bound
  survives the hierarchical bootstrap.
- Scale-up to other snRNA genes that can be aligned to a resolved family.
