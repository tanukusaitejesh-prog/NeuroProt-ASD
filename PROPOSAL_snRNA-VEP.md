# Proposal: snRNA-VEP — a variant-effect predictor for spliceosomal snRNA genes in neurodevelopmental disorders

**Status:** idea + novelty audit (October 2026). Data-pipeline scaffold in `snrna_vep/`.
**Self-rated novelty:** 4 / 5 · **Feasibility:** high (public data, one GPU, ~3–4 months for one person)

---

## 1. One-line pitch

Since 2024, small non-coding spliceosomal RNA genes have become one of the largest single causes of
NDDs: *RNU4-2* (ReNU syndrome), *RNU2-2*, *RNU5B-1*, *RNU5A-1*, with dominant and recessive forms.
Clinical labs still can't score a new variant in these genes computationally: CADD and the other
existing predictors fail, and the 2026 clinical guidance tells labs not to rely on them. We build
the first predictor made for these genes. It reads structure across the whole splicing cycle, uses
evidence from paralogs, and tells dominant variants from recessive ones. We test it with
leave-one-gene-out and time splits.

## 2. Why the gap is real (evidence)

| Evidence | Source |
|---|---|
| De novo variants in an 18-bp region of *RNU4-2* explain ~0.4% of NDD; one recurrent insertion (n.64_65insT) | Chen et al., *Nature* 2024 — https://www.nature.com/articles/s41586-024-07773-7 |
| U4 and U5 snRNA genes (*RNU5B-1*, *RNU5A-1*) cause NDD through splicing disruption | *Nat Genet* 2025 — https://www.nature.com/articles/s41588-025-02184-4 |
| *RNU2-2* and *RNU5B-1* are NDD genes | *Nat Genet* 2025 — https://www.nature.com/articles/s41588-025-02209-y ; https://www.nature.com/articles/s41588-025-02159-5 |
| Recessive *RNU2-2* is the most common known recessive NDD; recessive *RNU4-2* is a separate syndrome | *Nat Genet* 2026 — https://www.nature.com/articles/s41588-026-02539-5 ; https://www.nature.com/articles/s41588-026-02554-6 ; https://www.nature.com/articles/s41588-026-02547-5 |
| Saturation genome editing (SGE) of *RNU4-2*: SGE AUC **0.95** vs CADD AUC **0.65** for ReNU vs population variants; a CADD cutoff catching all ReNU SNVs also flags 56.4% of non-pathogenic ones | *Nature* 2026 — https://www.nature.com/articles/s41586-026-10334-9 |
| 2026 clinical guidance: CADD thresholds that catch every ReNU/RNU4ATAC SNV also flag 55–62% of gnomAD SNVs; structure tools including AF3 can't model the dynamic spliceosome context | medRxiv 2026.08.03.26359558 (Guidance for clinical variant classification in snRNA genes; introduces the RNUdb *database*, not a predictor) |
| Dominant U4/U6 variants also cause retinitis pigmentosa, so one gene can give different phenotypes depending on where the variant sits | *Nat Genet* 2025 — https://www.nature.com/articles/s41588-025-02451-4 |

SGE exists only for *RNU4-2*. Every other snRNA gene, and every new candidate, has no usable
computational evidence (ACMG PP3/BP4).

## 3. Novelty audit (what exists vs. what this adds)

Searched up to October 2026 by web search only (Nature/medRxiv/bioRxiv full texts were not reachable
from the cloud session, so these conclusions rest on abstracts and search snippets).

| Exists already | Does it cover snRNA VEP? |
|---|---|
| Generic noncoding VEPs: CADD, Evo 2, AlphaGenome, GPN-MSA, Gnocchi constraint | Not snRNA-specific. CADD is documented to fail. No published snRNA evaluation of Evo 2 or AlphaGenome found |
| AlphaGenome Atlas (medRxiv Sept 2026, 10.64898/2026.09.16.26363192): scores for all 9B SNVs | Generic and regulatory; abstract does not mention snRNAs (full text not yet read). **Strongest baseline to beat** |
| Minor-spliceosome review (Frilander, *RNA* 2025): maps RNU4ATAC/RNU12 variants onto cryo-EM structures | Descriptive mapping only, no predictive model. Prior art for the structure idea, cite it |
| RNU4-2 SGE functional map | Covers one gene, from a wet-lab assay. No model that transfers to other genes |
| RNUdb (2026) | A database of patient and population counts plus paralog-equivalent positions. No predictive model |
| RNA language models (RiNALMo, RNA-FM), RNA structure tools | Not applied to snRNA pathogenicity |
| snoRNA–disease association GNNs (GL4SDA, SPGA) | Gene–disease link prediction, not variant effect |

**What's new here (each piece unpublished as far as we found):**
1. **The task itself.** First variant-level pathogenicity model for spliceosomal snRNA genes, scoring SNVs *and indels*.
2. **Structure across the whole splicing cycle.** Each nucleotide gets a profile of RNA–RNA and RNA–protein
   contacts across cryo-EM states (tri-snRNP → pre-B → B → Bact → C → P), instead of one fold
   prediction, which the guidance says fails.
3. **Paralog-aware transfer.** Evidence moves between U4/U5/U6/U2 paralogs and the major/minor
   spliceosome analogs (U4↔U4atac, U6↔U6atac), weighted by each copy's expression.
4. **Mechanism and phenotype head.** Separate outputs for dominant NDD, recessive NDD, retinitis
   pigmentosa, and benign. Those regions are now known to differ.
5. **Time-split, leave-one-gene-out evaluation.** "Could RNU2-2 have been predicted from RNU4-2?"
   Train on variants published before 2025, test on discoveries from 2025–26.

### Still to check before investing (needs full-text access)
1. AlphaGenome Atlas full text: any snRNA / RNU4-2 / RNU2-2 evaluation?
2. snRNA guidance / RNUdb full text: any predictor or calibrated PP3/BP4 score?
3. RNU4-2 SGE paper: which in-silico tools were compared (CADD confirmed; others?)
4. bioRxiv/medRxiv/Scholar, last 6 months, plus papers citing the SGE paper: "snRNA variant effect
   prediction", "RNU4-2 machine learning", "RNU2-2 pathogenicity prediction", "ReNU classifier".

Decision rule: a multi-gene snRNA ML predictor already exists → pivot to the mechanism head and the
time-split test. Only single-gene models or benchmarks exist → cite them and beat them. Nothing → go.

## 4. Data (all public)

- **Patient variants:** supplements of the papers in §2, plus LOVD (lots of RNU4-2 entries), ClinVar,
  and RNUdb. Older genes: *RNU4ATAC* (MOPD1/Roifman/Lowry-Wood), *RNU12* (ataxia), *RNU7-1* (AGS),
  *RMRP* (CHH).
- **Functional ground truth:** RNU4-2 SGE function scores (Nature 2026 supplement / MaveDB).
- **Population controls:** gnomAD v4.1 genomes, All of Us allele-frequency browser, UK Biobank counts
  as reported in the papers/RNUdb. **Apply a mappability filter** because these genes are multicopy.
- **Structure:** PDB cryo-EM spliceosome states (human tri-snRNP, pre-B, B, Bact, C, P, U2 snRNP,
  minor spliceosome). Compute nucleotide contacts to proteins (PRPF8, PRPF31, SNU13, SF3B…) and to
  pre-mRNA and other snRNAs.
- **Sequence and evolution:** RNA LMs (RiNALMo, RNA-FM), Evo 2 (1B/7B) log-likelihood ratios,
  phyloP/Zoonomia, Gnocchi, plus covariation from Rfam alignments.
- **Paralog expression:** small-RNA-seq / ENCODE to weight which copies are actually expressed (for
  example, RNU2-2 was long called a pseudogene).

## 5. Model

Inputs for each candidate variant (gene, position, alt, including indels):
`[RNA-LM ΔLL, Evo2 ΔLL, conservation, covariation, multi-state contact profile, ΔMFE in each
state's local structure, paralog-aligned position embedding, paralog expression weight]`
→ small gated MLP or GNN over the snRNA contact graph (nodes are nucleotides and protein partners,
edges are contacts in each state)
→ heads: **(a)** regression on SGE score (auxiliary, RNU4-2 only), **(b)** multi-class
{dominant-NDD, recessive-NDD, RP, benign}.

Rules that keep the evaluation honest (lessons from NeuroProt-ASD):
- **Don't use AF/population counts as features** when gnomAD variants serve as the "benign" labels.
- Count each *unique* variant once. Recurrent n.64_65insT counts once, not 100+ times.
- Report everything per gene and leave-one-gene-out. No random splits across positions.

## 6. Evaluation and what "SOTA" means here

| Test | Metric | Baselines |
|---|---|---|
| Zero-shot on RNU4-2 SGE (model never sees SGE) | Spearman, AUROC for SGE-depleted variants | CADD 1.7, phyloP, Gnocchi, Evo 2, AlphaGenome / AlphaGenome Atlas AVI, RNA-FM, ΔMFE |
| Leave-one-gene-out: patient vs population (RNU2-2, RNU5B-1, RNU5A-1, RNU4ATAC, RNU12 …) | AUROC, AUPRC, sensitivity at 95% specificity (that's what PP3 needs) | same |
| **Time split**: train ≤2024, test on variants first published 2025–26 | same | same |
| Dominant vs recessive vs RP separation | macro-F1 | position-only and region-rule baselines |
| ACMG calibration (Pejaver-style) | thresholds giving supporting/moderate/strong PP3/BP4 | CADD (expected to fail) |

**Success = SOTA claim:** beat every generic VEP on leave-one-gene-out and time-split AUROC, and reach
calibrated PP3 "moderate" evidence for at least one gene without SGE.

**Discovery output:** score all ~2,000 snRNA/snoRNA/scaRNA loci in RNAcentral. Nominate candidate
"critical regions" in genes not yet linked to disease (for example the *RNU6* copies, *RNU1*,
*RNU11*, *RNU6ATAC*, *RN7SK*). Check them against patient variants published after the cutoff and
against de novo variants in public genome cohorts where accessible.

## 7. Plan (~14 weeks)

1. Wk 1–3: curate the patient and population variant table (largest manual step), with date stamps per variant.
2. Wk 3–5: features (contacts, LMs, conservation, paralog alignment); run the baseline VEPs.
3. Wk 6–9: model, LOGO and time-split evaluation, ablations (drop structure / paralog / LM).
4. Wk 10–12: calibration, genome-wide scoring, candidate regions.
5. Wk 13–14: write-up. Release the score table + web lookup (preprint early; the field moves fast).

## 8. Risks and stop criteria

- **Small labels:** few hundred unique pathogenic variants in ~10 genes. Mitigation: features that
  don't depend on the gene, auxiliary SGE head, low-capacity model. *Stop if* LOGO AUROC is no
  better than conservation alone (bootstrap CI overlaps).
- **Ascertainment bias:** patient variants cluster where people looked. Mitigation: time split, SGE
  as an unbiased test set.
- **Multicopy mapping:** gnomAD may miss or mis-map variants. Mitigation: mappability filter;
  sensitivity analysis on high-confidence sites only.
- **Getting scooped:** active labs (SGE / RNUdb groups) could publish a model. Mitigation: post a
  preprint early. The cross-gene, mechanism-aware, time-split angle stays distinct even if an
  RNU4-2-only SGE model appears.

## 9. Ideas checked and rejected (novelty too low)

| Idea | Why rejected |
|---|---|
| Better ASD gene prioritization (GNN / PLM / scRNA features) | Crowded: forecASD, DAMAGES, A-risk, DeepND, HyperAD, Baylor AJHG 2025 single-cell model |
| ML prior inside de novo likelihood (TADA-style) | Done by VBASS (Commun Biol 2023) and GeneBayes (2024) |
| Multi-disorder joint de novo model | mTADA (Nat Commun 2020) + NAR GB 2025 joint analysis |
| Leakage/time-split benchmark of ASD gene predictors | Useful but likely a ~3 novelty "benchmark" paper again, similar to NeuroProt-ASD's outcome |
| ABIDE fMRI / home-video / EEG ASD classifiers | Heavily saturated; site leakage already well documented |
