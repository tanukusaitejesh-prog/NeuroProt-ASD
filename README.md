# snRNA-VEP

Variant-effect prediction for spliceosomal snRNA genes in neurodevelopmental disorders
(*RNU4-2* / ReNU, *RNU2-2*, *RNU5B-1*, *RNU5A-1*, *RNU4ATAC*, …). See
[`PROPOSAL_snRNA-VEP.md`](PROPOSAL_snRNA-VEP.md) for the motivation, novelty audit and evaluation plan.

## Pipeline

| Step | Script | Needs | Output |
|---|---|---|---|
| 0 | `scripts/00_download_local.py` | internet (rcsb.org) | `data/raw/pdb/*.cif.gz` spliceosome structures |
| 1 | `scripts/01_build_catalogue.py` | `gsutil` (public GCS) | `data/processed/small_rna_catalogue.tsv` (GENCODE v27 small RNAs) |
| 2 | `scripts/02_fetch_gnomad.py` | public GCS | gnomAD v4.1 genome variants + reference for core genes ±2 kb |
| 3 | `scripts/03_depletion_scan.py` | step 2 | context-aware, within-gene O/E depletion (`results/depletion_*`, `figures/depletion_profiles.png`) |
| 4 | `scripts/04_build_features.py` | steps 2–3 (+0, +4b if available) | `data/processed/variant_features.tsv.gz`: every SNV and 1-nt indel with features |
| 4c | `scripts/04c_add_contacts.py` | `data/raw/pdb/` | refresh cryo-EM contact features without re-folding |
| 4d | `scripts/04d_add_patient_variants.py <patient tables>` | patient table | featurize patient variants outside the enumerated set (multi-nt indels) |
| 4e | `scripts/04e_position_conservation.py` | step 2 | interim position-level phyloP from gnomAD annotations (gaps interpolated, no missingness flag) |
| 4b | `scripts/04b_external_scores.py cadd phylop rnalm evo2 alphagenome=PATH` | internet / GPU | `data/external/*.tsv`; rerun step 4 to merge |
| 5a | `scripts/05a_clinvar_labels.py` | public GCS | interim patient labels from ClinVar 2025-05 (P/LP) |
| 5c | `scripts/05c_zero_shot_dominant.py --gene RNU4-2` | 5a | zero-shot feature check on a dominant gene + cross-gene transfer model |
| 5 | `scripts/05_train_eval.py --task dominant` | `data/curation/patient_variants.tsv` | leave-one-gene-out benchmark, paired bootstrap vs baselines, ablations, time split, SGE zero-shot |

```bash
pip install -r requirements.txt
python scripts/01_build_catalogue.py && python scripts/02_fetch_gnomad.py
python scripts/03_depletion_scan.py
python scripts/00_download_local.py            # structures
python scripts/04_build_features.py
python scripts/04b_external_scores.py cadd phylop rnalm   # evo2 on a GPU
python scripts/04_build_features.py            # merge external scores
python scripts/05_train_eval.py --task dominant
python -m pytest -q tests/
```

## Package

- `snrna_vep/coords.py`: HGVS n. ↔ GRCh38, left-normalized indels (minus-strand aware)
- `snrna_vep/variants.py`: enumerate all SNVs / 1-nt indels per gene with join keys
- `snrna_vep/mutation_model.py`: trinucleotide-context mutability from flanks; O/E depletion
- `snrna_vep/features_rna.py`: ViennaRNA ΔΔG (fold and U4/U6-type duplex), unpaired probability
- `snrna_vep/align.py`: paralog/analog alignment onto family reference coordinates
- `snrna_vep/contacts.py`: per-nucleotide protein/RNA contacts across cryo-EM states (gemmi);
  snRNA chains found by sequence, not hand-entered chain IDs
- `snrna_vep/external.py`: CADD, phyloP, AlphaGenome Atlas, RNA-LM and Evo 2 scores
- `snrna_vep/labels.py`: patient table schema, gnomAD control definitions (dominant vs recessive), SGE loader
- `snrna_vep/model.py`: gene-agnostic GBM / logistic model, in-fold paralog transfer, LOGO and time split
- `snrna_vep/evaluate.py`: AUROC/AUPRC/sens@95%spec, bootstrap and paired-bootstrap tests

## Evaluation rules

- gnomAD-derived features (`population` group) are off by default because controls come from gnomAD.
- Each unique variant counts once; patients are collapsed across papers (`collapse_patients`).
- Recessive task uses gnomAD homozygotes as controls, never heterozygous carriers.
- Uncallable multicopy loci (RNU1-1..4, RNU2-1) are excluded.

## First results (gnomAD only, no patient labels)

`scripts/03_depletion_scan.py`: snRNA genes carry 2–9× more gnomAD SNVs than a flank-trained context
model predicts (consistent with reported snRNA hypermutability), so depletion is measured relative to
each gene's own rate. Ranked by their most depleted 10-nt window, the four known dominant NDD snRNA
genes come 1st, 2nd, 4th and 5th of 23 callable genes (RNU2-2, RNU5B-1, RNU4-2, RNU5A-1; exact
rank-sum permutation p = 4/8855 ≈ 5×10⁻⁴). RN7SK, with no known disease link, ranks 3rd. In RNU4-2 the
minimum is at n.67, inside the ReNU critical region. These genes were partly discovered from this
signal, so this checks the pipeline rather than validating the model.

`scripts/06_label_free_checks.py`: ViennaRNA single-molecule features (|ΔΔG fold|, pairing probability)
do **not** track gnomAD depletion in RNU4-2, RNU2-2, RNU5A-1 or RNU5B-1 (|Spearman| < 0.1, n.s.); only
RNU4ATAC shows a weak signal (paired positions more depleted, ρ = −0.24, p = 0.006). Whole-molecule
RNAcofold U4·U6 ΔΔG even trends the wrong way, likely because the MFE dimer does not reproduce the
native U4/U6 pairing. This matches the clinical guidance that isolated-RNA folding is not informative
here and puts the weight on cryo-EM contact profiles, conservation and language-model features.

**Cryo-EM contact profiles** (`scripts/04c_add_contacts.py`; 11 human spliceosome structures from the
AWS PDB snapshot covering tri-snRNP, pre-B, B, Bact, C*, P, 17S U2 and the minor pre-B/Bact). snRNA chains
are matched by sequence (U1, U2, U4, U5, U6, U11, U12, U4atac, U6atac all found). Against within-gene
depletion, with a circular-shift permutation that respects autocorrelation along the RNA:

| Gene | snRNA–snRNA contact frequency across states: ρ (p_shift) | max protein contacts: ρ (p_shift) |
|---|---|---|
| RNU2-2 | −0.62 (0.006) | −0.21 (0.28) |
| RNU5A-1 | −0.30 (0.011) | −0.23 (0.09) |
| RNU5B-1 | −0.30 (0.011) | −0.20 (0.24) |
| RNU4-2 | −0.20 (0.24) | −0.19 (0.10) |

So RNA–RNA interactions across the splicing cycle carry constraint signal that isolated-RNA folding does
not, while protein-contact counts are suggestive but not significant after the permutation. Naive
Spearman p-values (`label_free_feature_vs_depletion.tsv`) overstate significance and should not be quoted.
This remains a label-free proxy; the decisive test is the patient-vs-population benchmark.

### First labelled results (ClinVar 2025-05 P/LP; interim, small)

**RNU4-2 zero-shot** (`results/zero_shot_RNU4-2_clinvar.tsv`; 18 ClinVar P/LP vs 215 gnomAD variants, no
training on RNU4-2):

| Score | AUROC [95% CI] |
|---|---|
| max protein contacts across cryo-EM states | **0.73** [0.64, 0.82] |
| fraction of states with protein contact | **0.73** [0.62, 0.82] |
| snRNA–snRNA contact frequency | **0.69** [0.58, 0.79] |
| logistic model trained on the other genes (transfer) | 0.61 [0.50, 0.72] |
| phyloP (241 mammals, position-level) | 0.34 (inverted) |
| ViennaRNA \|ΔΔG\| / duplex ΔΔG / pairing | 0.34–0.38 (inverted) |
| within-gene depletion | 0.98 (**circular**: controls are gnomAD variants) |

Conservation and isolated-RNA folding fail on RNU4-2, as the SGE paper reported for CADD (AUC 0.65 on its
own variant set); cryo-EM contact profiles are the only non-circular features above chance. Depletion
needs controls independent of gnomAD (UKB / All of Us counts from the papers) before it can be used.

**Recessive genes, leave-one-gene-out** (`results/benchmark_logo_recessive_clinvar.tsv`; RMRP, RNU4ATAC,
RNU12, RNU7-1; 105 P/LP vs 31 gnomAD homozygote controls): the trained model is near chance (GBM 0.54,
logistic 0.61); RMRP dominates the positives, and its pathogenic alleles are mostly insertions, so
|ΔΔG| alone reaches 0.68. This interim set is too small and too RMRP-heavy to judge the model; the
dominant spliceosomal benchmark needs the curated supplement table.

CADD could not be computed in the cloud session (no CADD host access; Hail does not build on this Python),
so the CADD comparison runs locally via `04b_external_scores.py cadd`.
