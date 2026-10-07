# Handoff: cloud session → local Claude Code

Read this first, then `PROPOSAL_snRNA-VEP.md` and `README.md`.

## Project in one paragraph
snRNA-VEP is a variant-effect predictor for spliceosomal snRNA genes that cause neurodevelopmental
disorders (*RNU4-2*/ReNU, *RNU2-2*, *RNU5B-1*, *RNU5A-1*, plus recessive *RNU4ATAC*, *RNU12*, *RMRP*, *RNU7-1*).
Generic predictors fail here (CADD AUC 0.65 vs saturation genome editing 0.95 on RNU4-2), and 2026 clinical
guidance says no in-silico tool is usable for these genes. Main idea: per-nucleotide contact profiles across
human spliceosome cryo-EM states, paralog transfer, and gene-held-out evaluation.

## State when the cloud session ended (October 2026)
- Pipeline in `snrna_vep/` and `scripts/` (steps 00–06; see README table). 7 offline tests pass.
- `data/processed/variant_features.tsv.gz`: 25,009 variants (all SNVs and 1-nt indels in 23 callable core
  genes, plus ClinVar multi-nt indels) with position, ViennaRNA, paralog, depletion, cryo-EM contact
  (11 structures) and interim phyloP (`phylop_gnomadpos`) features.
- Interim labels: `data/curation/patient_variants_clinvar.tsv` (ClinVar 2025-05 P/LP).
- The cloud session could not reach nature.com, medRxiv, NCBI, rcsb.org, HuggingFace or CADD hosts.
  Everything below needs normal internet.

## Results so far (details and caveats in README)
- Label-free: within-gene gnomAD depletion ranks the 4 dominant NDD snRNA genes 1st, 2nd, 4th, 5th of 23.
- snRNA–snRNA contact frequency tracks depletion in RNU2-2, RNU5A-1, RNU5B-1 (circular-shift p ≤ 0.011).
- RNU4-2 zero-shot (18 ClinVar P/LP vs 215 gnomAD): protein contacts AUROC 0.73 [0.64–0.82];
  phyloP 0.34 and ViennaRNA folding 0.34–0.38 (inverted); depletion 0.98 but circular.
- Recessive ClinVar leave-one-gene-out: model near chance (small, RMRP-heavy set; not the target task).

## To-do, in order
1. **Reproduce**: `python -m pytest -q tests/` then `python scripts/05c_zero_shot_dominant.py --gene RNU4-2`
   (expect max_protein_res ≈ 0.73).
2. **Decisive check**: `python scripts/04b_external_scores.py cadd phylop` → `python scripts/04f_merge_external.py`
   → rerun `05c`. Question: do contact features beat `cadd_phred` on RNU4-2? (CADD scores SNVs only.)
3. **Novelty re-check** (full text): AlphaGenome Atlas medRxiv 10.64898/2026.09.16.26363192; snRNA guidance /
   RNUdb medRxiv 10.64898/2026.08.03.26359558; RNU4-2 SGE Nature s41586-026-10334-9 (which tools compared);
   papers citing the SGE paper; bioRxiv/medRxiv last 6 months for "snRNA variant effect prediction".
4. **Curate** `data/curation/patient_variants.tsv` from the supplements of the papers listed in the proposal
   (template: `patient_variants.template.tsv`). Use GENCODE names (RNU2-2 = RNU2-2P), GRCh38, one row per
   variant per paper, first-publication month, cohort (to drop duplicate patients across papers).
   Also pull the RNU4-2 SGE scores into `data/external/sge_rnu4-2.tsv` (columns `key`, `sge_score`), and
   UK Biobank / All of Us counts as gnomAD-independent controls if the papers report them.
5. `python scripts/04d_add_patient_variants.py data/curation/patient_variants.tsv` →
   `python scripts/04e_position_conservation.py` → `python scripts/05_train_eval.py --task dominant`.
6. Optional: RNA language model (`pip install multimolecule torch`; `04b ... rnalm`), Evo 2 on a GPU.

## Rules to keep
- Population-derived features (`rel_oe_w10`, depletion) are circular whenever controls come from gnomAD.
- Count each unique variant once; never random-split positions; report per gene and leave-one-gene-out.
- Commit to branch `claude/neuro-disorders-aiml-idea-6ce9lt`.
