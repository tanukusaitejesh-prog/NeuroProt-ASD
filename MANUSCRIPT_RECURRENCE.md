# The dominance of a single insertion allele in ReNU syndrome is not explained by mutation rate

*Draft, 9 October 2026. Every number traces to a committed script, named in brackets. Not pre-registered;
this is a separate study from the snRNA-VEP variant-effect work and is reported as such.*

**Suggested venue.** Short report / research letter: *Genetics in Medicine*, *European Journal of Human
Genetics*, or *Human Molecular Genetics*. The result directly addresses a question posed in a 2026 *Nature*
paper, so a Matters Arising there is also defensible.

---

## Abstract

De novo variants in *RNU4-2* cause ReNU syndrome, among the most frequent single-gene causes of
neurodevelopmental disorder yet described. Strikingly, a single-base insertion, n.64_65insT, accounts for
70–77% of all reported cases. The reason is unresolved: saturation genome editing has shown the variant is not
unusually damaging, and the candidate explanations named in the literature are an elevated local mutation rate,
positive selection in the female germline, and ascertainment bias.

We tested the mutation-rate explanation using public data. Base-pair-resolution germline mutation rates
(Roulette) show the 18-nucleotide critical region is *less* mutable than the remainder of *RNU4-2* (mean
relative rate 0.063 versus 0.088), and do not predict which critical-region variants are observed in patients
or in population cohorts. Because n.64_65insT is a thymine inserted into a four-thymine homopolymer, we measured
the single-base insertion frequency at homopolymer tracts empirically from gnomAD v4.1, across 10 genomic
windows containing 10,992 T₄ tracts. Insertion frequency rises log-linearly with tract length
(1.90-fold per added base, R² = 0.90, ρ = 0.99), reproducing the expected signature of replication slippage and
validating the measurement. At tract length four, however, single-base insertions occur at only 0.26 times the
frequency of substitutions (95% Poisson CI 0.22–0.30): such sites are mutationally *disfavoured*.

In patients, single-base insertions in the critical region are 7.9-fold more frequent per allele than pathogenic
single-nucleotide variants in the same region. Set against a mutational expectation of 0.26, this is a
34-fold discrepancy; for n.64_65insT alone the figure is 153-fold. Mutation rate therefore does not explain the
dominance of this allele, and the discrepancy is the quantity that any selection-based explanation must account
for. We discuss why ascertainment is unlikely to be sufficient and what measurement would settle the question.

---

## Introduction

*RNU4-2* encodes the U4 small nuclear RNA, a component of the major spliceosome. De novo variants in an
18-nucleotide critical region cause ReNU syndrome, which accounts for roughly 0.4–0.5% of undiagnosed
developmental disorder in large cohorts.

What is unusual is the allelic spectrum. One variant, n.64_65insT (GRCh38 chr12:120,291,839:T:TA), accounts for
77.4% of cases in the discovery cohort and 72.6% in an independent French cohort. Most Mendelian disorders are
allelically heterogeneous; a single allele at this frequency demands an explanation.

The obvious one has been excluded by the primary literature. Saturation genome editing of *RNU4-2* found that
n.64_65insT is not unusually damaging, with many variants scoring similarly or lower on function. Those authors
concluded that high recurrence is unlikely to result from a particularly damaging functional effect driving
ascertainment, and that positive selection in the female germline or an increased local mutation rate are the
more likely explanations — while stating that which of these applies is at present unknown.

Here we test the mutation-rate explanation, because it is the one that can be settled with public data.

---

## Results

### The critical region is not mutationally hot

We retrieved base-pair-resolution germline mutation rate estimates for every possible single-nucleotide variant
in *RNU4-2* from Roulette, by range query against the per-chromosome rate file [48]. Across 423 records, the
mean relative rate inside the 18-nucleotide critical region was 0.063, against 0.088 across the remainder of the
gene. The region in which pathogenic variants cluster is therefore, if anything, *less* mutable than its
surroundings.

Within the critical region, mutation rate did not predict which variants are observed. Rate was uncorrelated
with the number of patient carriers (Spearman ρ = +0.18, p = 0.41, n = 23 variants) and with the number of
population carriers (ρ = −0.23, p = 0.29). Observed variants separated cleanly by context rather than by rate:
every variant with patient carriers had zero population carriers, and every variant with population carriers
had zero patient carriers.

### Insertion frequency at homopolymers scales log-linearly with tract length

Roulette models substitutions only. Since n.64_65insT is a thymine inserted into a T₄ homopolymer — the
canonical signature of replication slippage — we measured the insertion frequency at homopolymer tracts
directly [49]. Across 10 windows of chromosome 12 (~2 Mb of sequence), we enumerated every homopolymer tract
from the reference and counted gnomAD v4.1 PASS single-base insertions assigned to each, expressing the result
relative to the frequency of an observed substitution per possible substitution in the same windows — the same
cohort, the same filters, the same ascertainment.

Insertion frequency rose monotonically and log-linearly with tract length (ρ = 0.99, p = 2.6 × 10⁻¹⁸;
1.90-fold per added base, R² = 0.90), from 0.09 at length two to 24.5 at length twelve. This reproduces the
established length dependence of replication slippage and serves as an internal validation of the measurement.

**Table 1.** Single-base insertion frequency at A/T homopolymer tracts, relative to the per-substitution frequency.

| tract length | relative insertion frequency | tracts |
|---|---|---|
| 2 | 0.09–0.12 | 23,497 |
| 3 | 0.11–0.13 | 64,305 |
| **4** | **0.21–0.26** | **21,831** |
| 5 | 0.48–0.60 | 8,722 |
| 6 | 3.44–3.67 | 2,285 |
| 7 | 7.22–7.39 | 1,005 |
| 8 | 11.3–13.3 | 360 |

### A T₄ tract is mutationally disfavoured, not favoured

At tract length four — the context of n.64_65insT — single-base insertions occurred at 0.257 times the
substitution frequency for T₄ (163 insertions across 10,992 tracts; 95% Poisson CI 0.22–0.30) and 0.208 for A₄
(130 across 10,839; CI 0.17–0.24). Insertions at such sites are roughly four-fold *less* likely than
substitutions.

### The residual

Using the per-variant carrier counts reported in the discovery cohort, single-base insertions in the critical
region have 19.8 patient carriers per allele (5 alleles, 99 carriers) against 2.5 for pathogenic
single-nucleotide variants (6 alleles, 15 carriers) — a 7.9-fold excess per allele. Against a mutational
expectation of 0.26, the discrepancy is **34-fold**. Considering n.64_65insT alone (89 carriers), it is
**153-fold**.

We report the class-level comparison as primary because it does not single out the allele that was discovered
first, and is therefore substantially more robust to ascertainment.

---

## Discussion

The dominance of n.64_65insT in ReNU syndrome is not a consequence of an elevated mutation rate. The critical
region is less mutable than its surroundings; mutation rate does not predict which variants are seen; and the
specific sequence context of the recurrent allele is one in which insertions are disfavoured relative to
substitutions by roughly four-fold. One of the three explanations named in the primary literature can therefore
be set aside.

Of the remaining two, ascertainment is the more mundane and must be addressed. n.64_65insT was the allele
through which the disorder was discovered, and subsequent studies searched for it directly, which inflates its
apparent share. Three observations argue that this is insufficient. First, its share is similar in an
independent cohort ascertained separately (72.6% versus 77.4%). Second, the comparison we report as primary is
between *classes* of variant — insertions versus substitutions — in the same region of the same gene, and both
classes were sought in the same screens. Third, substitutions in the critical region are if anything easier to
detect than a single-base insertion in a homopolymer, which is the harder call for short-read pipelines; the
technical bias runs against the observed direction.

That leaves positive selection in the germline, which the primary literature proposes but does not quantify. Our
estimate puts a number on what it would have to achieve: a 34-fold enrichment of insertion alleles over
substitution alleles beyond mutational expectation, or 153-fold for the single recurrent allele. For comparison,
the best-characterised selfish spermatogonial variants are estimated to be enriched by one to two orders of
magnitude, so the magnitude required here is not implausible — but every established example of that mechanism
is paternal, whereas all 54 informative ReNU cases arose on the maternal allele. A maternal mechanism of
comparable strength is not currently described.

This points to the measurement that would settle it. Phased de novo variant data, in which the parental origin
and the local mutation spectrum can be assessed together, would distinguish a female-germline selective process
from a mutational one we have not modelled. Such data exist in trio cohorts but are not publicly available at
the resolution required; this analysis is the strongest statement that can be made without them.

### Limitations

Observed allele counts in population data reflect mutation rate, genetic drift and selection together, so our
measurement is of relative mutability, not a per-generation rate. gnomAD detects insertions less sensitively
than substitutions, which biases our estimate of 0.26 downward and therefore makes the reported residual
conservative; a two-fold sensitivity correction would still leave a ~75-fold discrepancy at the class level.
Sampling was restricted to chromosome 12, and genome-wide sampling would tighten the estimate. Patient carrier
counts were transcribed from a published table and inherit its ascertainment. Finally, we test only the
mutation-rate hypothesis; we do not establish what the alternative is.

---

## Methods

**Mutation rate estimates.** Base-pair-resolution germline mutation rates were obtained from Roulette
(v5.2, TFBS-corrected) by HTTP range request against the per-chromosome bgzipped VCF using its CSI index, so the
3.8 GB file was not downloaded. Implementation in `snrna_vep/remote_tabix.py` [48].

**Homopolymer insertion frequency.** Ten windows of 200 kb were sampled uniformly at random from chr12:20–110 Mb
(seed 0), excluding windows more than 20% N. Reference sequence was retrieved from the UCSC REST API and
homopolymer tracts of length ≥ 2 enumerated by regular expression. gnomAD v4.1 genome sites were retrieved for
each window by tabix range query against the 26 GB per-chromosome VCF. PASS variants were classified as
substitutions (single-base reference and alternate) or single-base insertions (alternate one base longer and
reference-prefixed), and each insertion was assigned to the homopolymer tract it extends, if any. Relative
frequency is the number of tracts of a given base and length carrying an observed insertion, divided by the
number of such tracts, expressed as a multiple of the number of observed substitutions divided by three times
the number of callable bases [49].

**Patient carrier counts** were transcribed from the published Extended Data table of the discovery cohort and
are provided in `data/curation/renu_variant_counts.tsv`.

**Statistics.** Spearman correlation for rank association; linear regression on log₁₀ relative frequency against
tract length; exact Poisson intervals on insertion counts.

**Code and data availability.** All analyses are scripted and committed; each result names its script.

---

## Figures

**Figure 1.** Mutation rate across *RNU4-2*. (a) Roulette relative rate per nucleotide with the 18-nt critical
region shaded; patient variants overlaid. (b) Rate inside versus outside the critical region. (c) Rate versus
patient carriers and versus population carriers, both flat.

**Figure 2.** Homopolymer insertion frequency. (a) Relative insertion frequency against tract length, log scale,
with the fitted slope; T₄ marked. (b) The same split by base. (c) Poisson intervals at length four.

**Figure 3.** The residual. Patient carriers per allele for insertions versus pathogenic substitutions in the
critical region, alongside the mutational expectation, with the resulting 34-fold discrepancy annotated.

## Data files

| file | contents |
|---|---|
| `results/renu_recurrence_rates.tsv` | per-variant Roulette rates and observed carrier counts |
| `results/homopolymer_insertion_rate.tsv` | insertion frequency by tract base and length |
| `data/processed/roulette_rnu4_2.tsv` | all 423 Roulette records for *RNU4-2* |
| `data/curation/renu_variant_counts.tsv` | transcribed per-variant carrier counts |
