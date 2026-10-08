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
