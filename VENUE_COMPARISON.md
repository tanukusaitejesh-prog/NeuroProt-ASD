# MuST-VEP against the two comparable papers in the target venue

Both comparators are *Briefings in Bioinformatics* 2026, Volume 27, Issue 3, "Problem Solving Protocol" — the exact
slot this manuscript targets. Both are from the same group (Drugparadigm Research Laboratory / Keshav Memorial
Engineering College, Hyderabad), both unfunded, both ~6 months from submission to acceptance.

- **InversePep** — Chilakamarri *et al.*, *Brief Bioinform* 2026;27(3):bbag277. doi:10.1093/bib/bbag277
- **SE(3)-PROTACs** — Kothakapu *et al.*, *Brief Bioinform* 2026;27(3):bbag228. doi:10.1093/bib/bbag228

---

## Head to head

| | InversePep (bbag277) | SE(3)-PROTACs (bbag228) | **MuST-VEP (this work)** |
|---|---|---|---|
| **Problem** | peptide inverse folding | PROTAC degradation, binary | snRNA variant effect, NDD genes |
| **Representation** | GVP-GNN + diffusion + Transformer | SE(3)-Transformer + ESM-2 8M | multi-state contact profiling, 11 cryo-EM structures |
| **Dataset size** | 38,280 train / 552 test | 1,979 | 25,201 scored / 4,789 labelled |
| **Baselines** | **2** (ProteinMPNN, ESM-IF1) | **7** | **13** |
| **Ablation rows** | **4** | **8** (two studies) | **17** |
| **Evaluation splits** | random only | random / cluster / temporal | **held-out-gene / position-block / temporal / prospective** |
| **Headline metric** | TM-score 0.51 mean, **0.48 median** | random AUROC 0.871 | dominant AUROC **0.810** |
| **Realistic-split metric** | *not reported* | **cluster 0.650, temporal 0.701** | **held-out-gene 0.766 [0.713, 0.811]** |
| **Statistical null model** | none | none | **degree-preserving shuffle, 200×, p = 0.005** |
| **Confidence intervals** | none | none | **position-block + hierarchical bootstrap on every comparison** |
| **Pre-registration** | none | none | **timestamped before any analysis (commit 350165e)** |
| **Orthogonal ground truth** | none ("future work") | none | **saturation genome editing, 485 measurements, zero-shot** |
| **Clinical calibration** | none | none | **ClinGen Bayesian → PP3_moderate on the conservative bound** |
| **Negative results reported** | none | ESM-3 ≯ ESM-2 | **4** (G2 gate, mechanism split, trajectory, partner identity) |
| **Failed pre-specified gate reported as failed** | n/a | n/a | **yes — and it became the main finding** |
| **Independent machine replication** | none | none | minor spliceosome (*RNU4ATAC*, *RNU12*, *RNU6ATAC*) |
| **Precomputed resource released** | model weights | model weights | **13,280 variants, 17 genes, with callability annotation** |

---

## What the comparators did that we must match

These are the format requirements, and they are the reason both papers are published:

1. **A named method with an acronym.** Done — MuST-VEP.
2. **An architecture/workflow figure as Figure 1.** Done — `figures/fig1_workflow.png`, and unlike both
   comparators two of its three panels are data rather than schematic.
3. **Tables of numbers where the method's row is bold and biggest in most columns.** Done — Tables 1–4.
4. **A "Key points" box.** Done.
5. **Methods before results.** Done.
6. **Code on GitHub.** Repository exists; needs a README with install and a worked example.

## What we have that neither comparator has

- A **pre-registration** committed before any evaluation, with gate criteria that could fail — and one that did.
- **Uncertainty quantification on every number.** Neither comparator reports a single confidence interval.
- A **null model.** Neither comparator has one, so neither can exclude that its geometry features work through
  degree or size alone.
- **Orthogonal ground truth.** Both comparators validate against a proxy (refold TM-score; a binary DC50 cutoff)
  and both defer wet-lab validation to future work. MuST-VEP transfers zero-shot to 485 real functional
  measurements that no model saw.
- **Four reported negative results**, including two extensions killed by their own pre-declared criteria.
- A **clinical calibration** in the framework practitioners actually use.

## Honest assessment of where we are weaker

- **Effective sample size.** Two of four dominant genes contribute 7 and 4 pathogenic variants. The hierarchical
  intervals in Table 1 are wide and several cross zero. SE(3)-PROTACs has 1,979 independent samples; we have
  4,789 rows but only 5 labelled genes, and gene is the unit of inference.
- **One functional assay, one gene.** *RNU4-2* is the only snRNA with a saturation screen, so the headline result
  cannot yet be replicated in a second gene.
- **No learned geometric encoder.** Both comparators' novelty is architectural. Ours is not: the features are
  hand-counted contacts. Tested honestly, the learned alternatives we tried did not help (Table 2), and we say so —
  but a reviewer primed on geometric deep learning may read this as a gap rather than as a finding.
- **Recessive prediction does not work**, for any method tested, and we report that rather than dropping the arm.

## Recommended framing

Lead with the finding, not the method. The paper is not "we built another predictor"; it is **"population
constraint and molecular function are not the same thing, and benchmarks built on curated labels cannot tell them
apart"** — demonstrated in the one non-coding gene class where both readouts exist for the same variants. The
predictor is the instrument that makes that measurement possible.

That framing also protects the weakest flank: if a reviewer discounts the AUROCs on sample-size grounds, the
label-versus-function result stands independently, because it is a within-variant comparison of two score types on
485 measurements, not a cross-gene classification claim.
