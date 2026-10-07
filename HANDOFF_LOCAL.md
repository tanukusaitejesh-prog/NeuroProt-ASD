# Handoff for local Claude Code (7 Oct 2026, 19:30 UTC)

Read first: RESULTS_TONIGHT.md (all results, in order), MANUSCRIPT_DRAFT.md (draft v0.1), PROPOSAL_snRNA-VEP.md.

## State
- Final label set: data/curation/labels_v5.tsv (build: `python scripts/05e_labels_v2.py --v5`; v5c = conflicts removed).
- Model scripts read the label set from the env var SNRNA_LABELS (default labels_v2.tsv); results get a tag suffix
  (e.g. results/v2_heldout_labels_v5.tsv). Always run with `SNRNA_LABELS=labels_v5.tsv`.
- Final system = v2-nested (scripts/11_v2_model.py): v2 for dominant, contacts+phyloP / contacts for recessive.
- Done: held-out genes, SGE zero-shot (12), prospective ClinVar (17), paired tests (16), mechanism (18a/18b),
  score resource (19), figures (15), ablation multi- vs single-structure (20: ALL 0.806 vs oracle single 0.795 vs
  nested single 0.597 for dominant).
- Raw data (data/raw/: PDB mmCIFs, ClinVar VCF) is gitignored; scripts 01-04 re-download it. UKB AFB raw exports
  (data/external/ukb_afb/*.csv) are not in git; the derived table data/curation/ukb_afb_homozygotes.tsv is.
- RNA-FM scoring used a separate venv (multimolecule 0.0.6, transformers 4.48.3; weights remapped by hand);
  its output data/external/rnafm.tsv is committed.

## Next steps (reviewer-driven; aim: novelty ~4)
1. Run `SNRNA_LABELS=labels_v5.tsv python scripts/21_robustness.py` (position-block + gene-level bootstrap,
   ViennaRNA/constraint baselines + trained secondary-structure model, within-gene dominant vs recessive test).
   Record in RESULTS_TONIGHT.md; commit.
2. Allele-specific structural scoring: distinguish alleles at a position (base-pair preservation incl. wobble in each
   state; base vs backbone protein contacts); validate on within-position allele ranking with RNU4-2 SGE (all 3 alleles
   measured). Compare with CADD.
3. Generalisation to another RNA machine: RMRP (RNase MRP; cartilage-hair hypoplasia; ClinVar P/LP) with human RNase
   MRP cryo-EM structures, zero-shot from the spliceosome-trained model.
4. Manuscript fixes: soften abstract (no p<0.001 headline; use clustered CIs), drop BP4 claim, cite MitoTIP
   (paralog/position transfer precedent), Nava et al. 2025 (RBM42 region n.68-70; Snu66/SNRP27K/RBM42 stabilise the
   ACAGAGA loop; their SVM is a methylation episignature classifier, not a variant predictor), SGE / recessive RNU4-2 /
   RNU2-2 papers (qualitative dominant/recessive mechanisms). Add ablation + robustness results.
5. Post a medRxiv preprint soon (scoop risk: Whiffin-lab guidance preprint, 4 Aug 2026, calls for snRNA tools).
