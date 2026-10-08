# Judging the two published papers, and the idea that comes out of them

*Written 2026-10-09 overnight. Both PDFs read in full.*

---

## Part 1 — What these two papers actually are

Both are from the **same lab**: Vani Kondaparthi, Drugparadigm Research Laboratory / Keshav Memorial
Engineering College, Hyderabad. Both are in **Briefings in Bioinformatics 2026, Vol 27 Issue 3** —
article numbers bbag277 and bbag228, i.e. back to back. Both are in the **"Problem Solving Protocol"**
category. Both state **"The research did not receive any specific grant."** Both have undergraduate CS
students as first authors.

Timeline: InversePep received 8 Nov 2025, accepted 10 May 2026. SE(3)-PROTACs received 13 Oct 2025,
accepted 15 Apr 2026. **Six months from submission to acceptance, for both.**

### InversePep — honest read

**What it is.** Structure-conditioned inverse folding for peptides: given a 3D backbone, generate an
amino-acid sequence that folds into it. GVP-GNN encoder → diffusion over sequence space → Transformer
decoder with AdaLN conditioning.

**Novelty, judged as of May 2026: 2 / 5.**
Every component is off the shelf and was off the shelf well before submission:
- GVP-GNN — Jing et al. 2021
- diffusion with self-conditioning — Chen, Zhang & Hinton 2022 (their ref 34)
- AdaLN conditioning — DiT, Peebles & Xie 2023 (ref 40)
- bucket batching — ref 24
- the task itself — ProteinMPNN 2022, ESM-IF1 2022

The assembled novelty is "apply sequence diffusion inverse folding to peptides rather than proteins, and
score with refold-TM instead of sequence recovery." That reframing of the metric is the one genuinely
sensible idea in the paper.

**Metrics, honestly: 2 / 5.**
- Mean TM-score **0.51**, median **0.48**. TM-score 0.5 is the conventional "same fold" threshold. At
  the median they are *below* it. The paper's own text concedes this: "TM-scores above 0.5 suggest that
  the overall fold and topology are correctly captured."
- Sequence recovery in their own Figure 4: **0.11, 0.04, 0.20**. That is essentially no recovery. They
  reframe it as diversity.
- Baselines are ProteinMPNN and ESM-IF1 — **both trained on proteins, not peptides**, while InversePep
  was trained on peptides. The comparison is structurally unfair in their favour and they do not say so.
- Physicochemical validation: **n = 5** peptides (Table 3).
- No wet lab, no MD. Explicitly deferred: "Future work will focus on ... validating findings in wet-lab
  experiments."
- Random split only. No cluster split, no temporal split.

**Potential / impact: 2 / 5.** Code is released. It will be cited as "another inverse folding model."

**What got it published:** a clean single problem, a named method, two baselines, a four-variant ablation
study, an architecture figure, three results tables, and a GitHub link.

### SE(3)-PROTACs — honest read

**What it is.** Binary prediction of whether a PROTAC degrades its target. SE(3)-equivariant transformer
over the PROTAC molecular graph + ESM-2 embeddings of target and E3 ligase + a "pairwise interaction
mechanism" (all-pairs target×ligase scoring, mean-pooled, sigmoid-gated) → MLP.

**Novelty, judged as of Apr 2026: 2.5 / 5.** The SE(3)-Transformer and ESM-2 are off the shelf. The
pairwise interaction module is the one original contribution and it is about six lines of algebra. But it
is a *real* contribution, and the ablation (Table 5) does show it beats multi-head attention, gated
fusion, tripartite attention and TAN.

**Metrics, honestly: 3 / 5 — and this is the interesting part.**
- Random split: accuracy 0.808, AUROC 0.871. Looks good.
- **Cluster split (cold target, n = 128): accuracy 0.656, AUROC 0.650.**
- **Temporal split (post-March 2023): accuracy 0.641, AUROC 0.701.**
- On the cluster split, **RF + Morgan fingerprints beats them on AUPR** (0.6366 vs 0.5886). They argue
  around it via MCC.
- Dataset: **1,979 PROTACs.** Small.
- ESM-2 **8M** — the smallest ESM-2 there is, 6 layers. Justified by dataset size; the ablation shows
  ESM-3 does not help.
- No wet lab.

So on the only two realistic evaluations, AUROC is 0.65 and 0.70 — barely better than useless. **But they
ran those splits and published the bad numbers.** Seven baselines. Two separate ablation studies. Six
metrics. That honesty is exactly why it is publishable, and it is the single best thing in either paper.

**Potential / impact: 2.5 / 5.** Positioned as a "computational pre-filter", which is the right, modest
claim.

---

## Part 2 — The calibration you actually need

Read that again and hold it next to your own work.

Two papers got into a **Q1 journal** (BiB, impact factor ~6.8, JCR Q1 in Mathematical & Computational
Biology) with: zero new biology, zero experimental validation, off-the-shelf architectures recombined,
datasets of 552 and 1,979, performance that is weak under realistic evaluation, and no funding.

What they *did* have, in both cases, was a **format**:

1. one clean nameable problem
2. a named method with an acronym
3. baselines in a table (2 and 7)
4. an ablation table (4 and 8 variants)
5. an architecture figure
6. multiple evaluation splits — SE(3)-PROTACs has random / cluster / temporal
7. code on GitHub

That is the whole recipe. It is a format, not a discovery.

**Your snRNA-VEP work already beats both of these on substance and loses to them on format.** You have a
timestamped pre-registration, held-out-unit evaluation, degree-preserving shuffled nulls, an
ascertainment-free orthogonal functional readout (saturation genome editing), ClinGen Bayesian ACMG
calibration, and honestly-reported negative results. Neither of these papers has *any* of those things.
What you don't have is an acronym, a seven-row baseline table, an ablation table, and an architecture
figure.

So there are two separate jobs, and they are not the same job:

- **Job A (cheap, certain):** put the existing snRNA work into this format. That is mostly writing.
- **Job B (the real one):** find the discovery that makes it more than a format exercise.

Part 3 is Job B.

---

## Part 3 — The idea: conformational hinges

### The observation these papers point to

The engine of both papers is the same: **encode 3D geometry equivariantly, bring in a pretrained language
model for the sequence partner, and add one named module that couples two things nobody has coupled.**
SE(3)-PROTACs couples the *two protein interfaces* via the PROTAC scaffold. That coupling module is where
all its performance comes from.

Ask the obvious question: for the spliceosome, **what is the axis nobody has coupled?**

It is not space. It is **time**.

### The gap

Gate G1 established that pooling 11 cryo-EM states beats any single state (0.815 vs 0.631, beats 200/200
degree-preserving nulls, p = 0.005). But it pooled them as an **unordered bag**: the feature was "fraction
of states in which this nucleotide touches a protein."

That throws away the one thing the spliceosome uniquely provides among all macromolecular machines. The
11 structures are not independent snapshots. **They are an ordered assembly pathway**, and every
nucleotide has a *trajectory* through it:

```
U2-snRNP -> tri-snRNP -> pre-B -> B -> Bact -> Bact-mature -> C* -> P      (major)
                      minor-preB -> minor-Bact                             (minor)
```

U4 is present in tri-snRNP, pre-B and B — and **gone from Bact onward**, because U4/U6 unwinding by BRR2
*is* the activation step. U5 loop I holds the exons through both transesterifications. U2 recognises the
branch point early. These are different jobs at different times, and no variant effect predictor in
existence — CADD, SpliceAI, phyloP, gnomAD constraint, RNA-FM, or your own v2-nested model — knows any of
this.

### The hypothesis

> **Pathogenic snRNA variants concentrate not at the most buried nucleotides, but at nucleotides whose
> protein partners turn over between consecutive states — the positions that must release one partner and
> acquire another for the cycle to proceed.**
>
> A permanently buried nucleotide is a constitutive scaffold: structurally important, but either
> mutationally lethal-early or compensated. A **rewired** nucleotide is a kinetic checkpoint, and that is
> where disease lives.

This is mechanistic, falsifiable, and NDD-specific. It predicts, without being told, that the RNU4-2
critical region (n.64–90 — the ReNU syndrome hotspot) should sit in the part of U4 that is remodelled in
preparation for ejection, not in the part that is most deeply buried while U4 is bound.

It also gives the recurrent `n.64_65insT` allele a mechanism it currently lacks: a T4 homopolymer sitting
in the stretch that has to be unwound.

### Why this is a discovery and not a benchmark

The headline is biology: *disease variants in spliceosomal snRNAs sit at conformational hinges, not at
static cores.* The model is the instrument, not the point. That is the framing you have been asking for,
and it is the thing neither published paper has at all.

It also comes with a **built-in replication across an independent machine.** The minor spliceosome
(U11 / U12 / U4atac / U6atac) is a separate evolutionary realisation of the same mechanism and causes a
*different* NDD spectrum — RNU4ATAC gives Taybi-Linder / Roifman / Lowry-Wood, RNU12 gives early-onset
cerebellar ataxia. You have 8Y6O (minor pre-B) and 7DVQ (minor Bact), and 83 pathogenic variants across
RNU4ATAC, RNU6ATAC, RNU11 and RNU12. If the hinge rule holds in both machines, that is a replication no
reviewer can wave away.

### The gate, declared before running

`scripts/51_trajectory_gate.py`. Per `(family, ref_pos)`:

| feature | meaning |
|---|---|
| `burial` | fraction of its resolved states with >=1 protein partner — **this is the G1 feature, the thing to beat** |
| `mean_degree` | mean number of distinct protein partners |
| `rewiring` | mean Jaccard **distance** between partner sets at consecutive states |
| `max_turnover` | largest single-transition Jaccard distance |
| `n_lost_max` | most partners released at any one transition (directional: preparing for ejection) |
| `core_frac` | fraction of partners present in *every* state (constitutive core) |

Three tests:

- **T1a** pathogenicity, leave-one-gene-out: `burial + mean_degree` vs the same **+ trajectory**.
- **T1b** measured function, RNU4-2 SGE, 485 variants, zero-shot: partial Spearman of
  `rewiring | burial` and the reverse. **This is the test that matters** — it is the only
  ascertainment-free readout, and it is the one that killed the constraint baseline before.
- **T1c** null: permute partner *identities* among positions within each structure, **matched on partner
  count**, so that `burial` and `mean_degree` are bit-identical and only cross-state coherence of partner
  identity is destroyed. 200 draws.

**Kill criterion, stated in the script before it ran:** if rewiring adds no held-out AUROC in T1a *and*
its partial correlation with SGE function given burial is indistinguishable from 0 in T1b, the hinge
hypothesis is dead and this line stops. No rescuing.

### RESULT: the gate FAILED. The hypothesis is dead.

Run 2026-10-09. `scripts/51_trajectory_gate.py`, `results/trajectory_gate_*.tsv`.

**T1b — measured function (RNU4-2 SGE, 485 variants, zero-shot):**

| feature | rho vs function | partial given `burial` | partial given `mean_degree` |
|---|---|---|---|
| `burial` | **+0.245** | — | — |
| `mean_degree` | **+0.298** | — | — |
| `rewiring` | +0.150 | **+0.092** | +0.088 |
| `max_turnover` | +0.154 | +0.094 | +0.091 |
| `n_lost_max` | +0.182 | +0.088 | +0.091 |
| `core_frac` | +0.134 | −0.092 | −0.140 |

Trajectory features correlate with function on their own, but almost all of it is burial. Partial correlation given
burial is **+0.092**, below the +0.10 kill line declared in the script before it ran. Meanwhile `mean_degree` given
`rewiring` is **+0.299** — the count explains the trajectory, not the other way round.

**T1c — the null is decisive and it is worse than "no effect":**

| null | observed | null mean ± sd | p |
|---|---|---|---|
| degree-matched partner-identity shuffle | +0.150 | **+0.281 ± 0.031** | **1.000** |
| free partner-set shuffle | +0.150 | +0.182 ± 0.039 | 0.801 |

Randomly permuting which position holds which partner set produces a **stronger** rewiring–function correlation than
the real data. The real trajectory is not merely uninformative, it is *less* coherent than chance under this measure.
There is no signal to rescue.

**T1a — pathogenicity, leave-one-gene-out (n = 813, 212 pathogenic, 5 genes):** `burial + mean_degree` 0.574 →
`+ trajectory` 0.603. The +0.030 looks like a gain until it is broken out by gene: RNU5B-1 0.526 → **0.960** on
*seven* positives, while RNU2-2P loses 0.101 and RNU12 loses 0.188. That is variance on a tiny positive count, not
an effect.

**Verdict: dead, by the criterion written down before the run. Not revisited.**

### Second gate, same night: does partner IDENTITY transfer? Also dead.

`scripts/52_partner_identity_gate.py`. The existing per-protein Mantel–Haenszel ORs are large — PPIL2 31.8,
SF3B1 22.2, SRRM2 17.5, PRPF8 7.7, SART1 6.7, SNRNP200 4.7 — and they are the obvious hook for the SE(3)-PROTACs
move of embedding each partner protein with ESM-2. But they are stratified by gene only, and pathogenic variants
cluster *within* a gene, so any protein whose footprint covers the cluster gets a large OR for free.

The test that removes the confound is whether identity transfers to a **held-out gene**:

| model | LOGO AUROC |
|---|---|
| `burial + mean_degree` | **0.574** |
| identity only | **0.437** (below chance) |
| `burial + mean_degree` + identity | **0.548** (worse than base) |

Only **9** protein partners appear at labelled positions in ≥3 of the 5 labelled genes, and **8 of those 9 are the
Sm ring** (SNRPB, SNRPD1/D2/D3, SNRPE, SNRPF, SNRPG) — which binds every snRNA by construction and so says nothing.
On RNU4-2 SGE, only PRPF8 occupancy survives conditioning on partner count (+0.112); everything else is ~0.

**Conclusion: the large per-protein ORs are positional clustering. Partner identity carries no transferable
information, and an ESM-2 partner-embedding layer is not worth building.** That is a useful negative — it says the
most load-bearing component of SE(3)-PROTACs does not port to this problem.

### What this means for the model

Two independent attempts to find structure beyond "how many proteins touch this nucleotide, in how many states"
both failed against their own pre-declared controls. Combined with G2 (constraint ties on clinical labels) and the
mechanism split (wrong sign), that is now **four** interpretive layers that have died while every *measurement* has
survived.

The honest reading is not that the next idea will work. It is that this dataset — 5 labelled genes, 1 gene with
functional data, 11 structures — has given what it has to give. The surviving asset is real and it is enough.

### The model, built their way but better

Both candidate "new modules" are now empirically excluded, which actually **simplifies** the build and makes it more
honest than either published paper. The model to ship is the one that survived every control:

- **multi-state contact profiling** across the 11 structures (G1: 0.815 vs 0.631 single-state, beats 200/200
  degree-preserving nulls, p = 0.005)
- **allele-specific pairing geometry** (within-position rho +0.346 vs CADD −0.077)
- **zero-shot transfer to measured function** (SGE rho 0.437 vs constraint 0.206, CADD 0.072, phyloP −0.004;
  partial structure|constraint +0.394 against constraint|structure +0.025)

and the two dead layers become the **ablation section** — the part SE(3)-PROTACs got credit for. Reported as
negatives:

- trajectory/hinge features add nothing beyond multi-state counts, and fail a degree-matched null (p = 1.00)
- partner identity does not transfer across genes (held-out AUROC 0.574 → 0.548); the large per-protein ORs are
  positional clustering
- an SE(3)-equivariant learned encoder, if built, should be reported against the hand-crafted counts whatever the
  outcome — "the learned geometry does not beat six hand-counted features" is a *better* ablation result than
  InversePep's, and it is true

That is a stronger ablation table than either paper has, because every entry is a real control rather than a
component swap.

### The comparison table that writes itself

| | InversePep | SE(3)-PROTACs | this |
|---|---|---|---|
| geometric encoder | GVP-GNN | SE(3)-Transformer | multi-state contact profiling over 11 cryo-EM structures |
| sequence partner model | — | ESM-2 8M | tested and **excluded** (identity does not transfer) |
| coupling module | — | pairwise (spatial) | tested and **excluded** (trajectory fails its null) |
| ground truth | refold TM-score (proxy) | DC50 binary (proxy) | **saturation genome editing, 485 real measurements** |
| evaluation splits | random | random / cluster / temporal | **random / held-out-gene / position-block / prospective** |
| null models | none | none | **degree-matched partner-identity shuffle, 200x** |
| pre-registration | none | none | **timestamped, SHA-256** |
| orthogonal validation | none ("future work") | none | **zero-shot transfer to a wet-lab assay** |
| clinical calibration | none | none | **ClinGen Bayesian ACMG -> PP3 strength** |
| independent replication | none | none | **minor spliceosome, separate machine, separate disease** |

Every row is a win. That table is the paper.

---

## Honest risk register — resolved

| risk | called at | outcome |
|---|---|---|
| rewiring adds nothing beyond burial | ~35% | **realised.** Partial +0.092, below the line; fails its null at p = 1.00. |
| "hinge" is a re-parameterisation of G1 | the real intellectual risk | **realised.** The degree-matched null settled it, exactly as designed. |
| partner identity is positional clustering | not called in advance | **realised.** Held-out-gene AUROC 0.574 → 0.548. |
| 4,789 labels over ~6 genes is tiny for a deep model | high, but normal here | stands. SE(3)-PROTACs trained on 1,583. |
| SGE exists for one gene only (RNU4-2) | known, unavoidable | stands. Still the only ascertainment-free readout anyone has. |

Total cost of both kills: about two hours, nothing built on top. That is what the gates were for, and they worked.

---

## Part 4 — The recommendation

Fourteen-odd idea attempts into this project, the pattern is consistent and worth saying plainly: **every
interpretive claim has died on its controls, and every measurement has survived.** Tonight added two more deaths to
that ledger. The dataset — 5 labelled genes, 1 gene with functional data, 11 structures — is not withholding a
discovery that one more clever framing will unlock. It has given what it has.

That is not the bad news it sounds like, because of Part 2. **The two published papers prove the bar is a format,
not a discovery.** They are in a Q1 journal with no new biology, no validation, off-the-shelf architectures, 552 and
1,979 samples, and AUROC 0.65 on their realistic splits.

So the recommendation is: **stop looking for idea #15. Do the format port on the work that survived.**

Concretely, for `MANUSCRIPT_v2.md`:

1. **Name the method and give it an acronym.** Non-negotiable; both papers have one.
2. **Build the baseline table to 7+ rows.** You already have CADD, phyloP447, gnomAD constraint, gnomAD density,
   ViennaRNA ddG, Rfam element, paralog density, distance-to-critical-region, RNA-FM (`data/external/rnafm.tsv`),
   and the published Benoit-Pilven 2020 structure score. That is ten. SE(3)-PROTACs has seven.
3. **Build the ablation table.** Tonight gave you three honest rows for free; add leave-one-state-out, drop
   allele features, drop the recessive interaction.
4. **Report four splits**, where SE(3)-PROTACs reports three: random, **held-out-gene**, position-block, and the
   **prospective timestamped** set (2,777 predictions, SHA-256 `b5b2249a...` — neither paper has anything like it).
5. **Draw the architecture figure.** Both papers have one and it does a lot of work.
6. **Lead with the result that survived everything:** structure predicts *measured* function while population
   constraint does not (partial +0.394 vs +0.025). Neither published paper has any orthogonal ground truth at all.

The comparison table in Part 3 is the honest case, and every row of it is already true.

**Target:** Briefings in Bioinformatics, "Problem Solving Protocol" — the exact slot both of these occupy, ~6 months
to acceptance. On substance you are ahead of both. The gap was never the science.

**Second paper, unchanged:** `MANUSCRIPT_RECURRENCE.md` stands on its own — the 20× residual after long-read
rate correction is a measurement, and measurements have been the thing that survives here.

## What is explicitly NOT being re-proposed

Ruled out from memory and prior work, do not revisit:
- benchmark/evaluation-only framings (rejected twice)
- control-set circularity for snRNA — tested; constraint survives UKB controls (0.743 vs 0.785 on
  RNU4-2), the circularity defence does **not** rescue G2
- PP3 / PM1 / PP2 dependence — that is the *other* paper, 20 drafts deep in `ndd_pp3_gate`
- dominant vs recessive mechanism split — pre-specified criterion not met, wrong sign
- de novo-trained classifiers of any kind — fail to transfer, three separate times
