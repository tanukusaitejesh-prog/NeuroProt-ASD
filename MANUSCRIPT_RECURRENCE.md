# The recurrent ReNU syndrome allele is mutationally disfavoured and exclusively maternal in origin

*Draft 2, 9 October 2026. Every number traces to a committed script, named in brackets. Not pre-registered;
this is a separate study from the MuST-VEP variant-effect work and is reported as such.*

**Suggested venue.** Short report / research letter: *Genetics in Medicine*, *European Journal of Human
Genetics*, or *Human Molecular Genetics*. The result directly addresses a question posed in a 2026 *Nature*
paper, so a Matters Arising there is also defensible.

---

## Abstract

De novo variants in *RNU4-2* cause ReNU syndrome, among the most frequent single-gene causes of
neurodevelopmental disorder yet described. A single-base insertion, n.64_65insT, accounts for 70–77% of all
reported cases. The reason is unresolved: saturation genome editing has shown the variant is not unusually
damaging, and the explanations named in the primary literature are an elevated local mutation rate, positive
selection in the female germline, and ascertainment bias.

We quantified the mutational supply of this allele from public data. Base-pair-resolution germline mutation
rates (Roulette) show the 18-nucleotide critical region is *less* mutable than the remainder of *RNU4-2* (mean
relative rate 0.063 versus 0.088), and do not predict which critical-region variants are observed in patients
or in population cohorts. Because n.64_65insT is a thymine inserted into a four-thymine homopolymer, we measured
the single-base insertion frequency at homopolymer tracts empirically from gnomAD v4.1, across 10 genomic
windows containing 10,992 T₄ tracts. Insertion frequency rises log-linearly with tract length (1.90-fold per
added base, R² = 0.90, ρ = 0.99), reproducing the expected signature of replication slippage and validating the
measurement. At tract length four, however, single-base insertions occur at only 0.40 times the frequency of
substitutions — measured in both short-read gnomAD (0.23) and the long-read HPRC pangenome (0.40), of which we
adopt the conservative long-read figure. Such sites are mutationally *disfavoured*.

The excess proves to be specific to a single allele rather than to insertions as a class: excluding
n.64_65insT, the per-allele insertion-to-substitution ratio across the critical region is exactly 1.00. That
permits a position-matched test, which holds assay coverage, calling behaviour and ascertainment constant by
construction. Four variants are observed at chr12:120,291,839, and n.64_65insT accounts for **89 of the 93
patient carriers** there. Tested against the substitution at the same nucleotide — the comparison to which the
measured mutational supply applies — the recurrent allele takes 89 of 91 carriers (97.8%, 95% CI 92.3–99.7%)
where 28.7% is expected: exact binomial p = 1.3 × 10⁻⁴⁵, a **110-fold** position-matched residual. The two other insertions at
that same nucleotide, which are the same class of variant call in the same homopolymer, carry two patients
between them, excluding insertion-calling bias as the explanation. We further exclude interlocus gene
conversion, the standard account of a recurrent allele with insufficient mutational supply: the tandem paralogue
*RNU4-1* lies 1.2 kb away and is 97.2% identical, but carries the *identical* T₄ tract and differs from
*RNU4-2* at no position inside the critical region, so conversion can generate none of the eleven observed
patient alleles.

What remains is an allele-specific germline excess that is, in every informative case reported, **maternal**.
The only documented mechanism capable of producing a germline excess of this magnitude — selfish spermatogonial
selection — is paternal and depends on clonal expansion across hundreds of spermatogonial divisions accumulated
over decades. Oogonial proliferation is confined to fetal development and to roughly thirty divisions, after
which the oocyte pool arrests, so an equivalent process would have to act far faster, within a bounded prenatal
window, and reproducibly across unrelated women. No documented example of germline selection operates this way.
We set out the three measurements that would identify the one that does.

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

Here we quantify the mutational supply of the allele, because that is the term that can be settled with public
data, and show that once it is measured the remaining discrepancy has a property that constrains its
explanation severely.

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
| 2 | 0.09–0.11 | 149,527 |
| 3 | 0.11–0.13 | 64,305 |
| **4** | **0.21–0.26** (short read) / **0.26–0.54** (long read) | **21,831** |
| 5 | 0.48–0.60 | 8,722 |
| 6 | 3.44–3.67 | 2,285 |
| 7 | 7.22–7.39 | 1,005 |
| 8 | 11.3–13.3 | 360 |

### A T₄ tract is mutationally disfavoured, not favoured

At tract length four — the context of n.64_65insT — short-read data gave 0.257 times the substitution frequency
for T₄ (163 insertions across 10,992 tracts; 95% Poisson CI 0.22–0.30) and 0.208 for A₄ (130 across 10,839;
CI 0.17–0.24).

Because short-read pipelines are known to under-call insertions in homopolymers, we repeated the identical
measurement in the same windows against the HPRC pangenome, whose variants derive from long-read haplotype
assemblies [50]. Long-read data gave 0.543 for T₄ and 0.262 for A₄, a mean of 0.403 — 1.73-fold above the
short-read estimate. The ratio of long-read to short-read rate grew with tract length (1.3–2.1-fold at lengths
2–4, rising to 6–7-fold at lengths 10–12), the expected signature of short-read homopolymer failure, and is
modest in the length range that matters here. We therefore adopt 0.40 as the expected relative insertion
frequency at a T₄ tract. Insertions at such sites remain roughly 2.5-fold *less* likely than substitutions.

### Interlocus gene conversion from the tandem paralogue is excluded

When a recurrent allele cannot be explained by mutation rate, the standard alternative is non-allelic gene
conversion: a near-identical donor elsewhere in the genome repeatedly overwrites the acceptor, so the "mutation"
is a copying event and its frequency is set by recombination rather than by polymerase error. This mechanism
accounts for recurrent pathogenic alleles in *CYP21A2* (from *CYP21A1P*), *PMS2* (from *PMS2CL*), *SMN1* (from
*SMN2*) and *GBA* (from *GBAP1*).

*RNU4-2* has precisely the configuration that makes this plausible. Its paralogue *RNU4-1* lies 1,194 bp away
on the same strand of chromosome 12, is the same length (141 nt), and is 97.2% identical (137/141).

The hypothesis is nonetheless decisively excluded, because conversion can only transfer sequence the donor
actually carries [56]. *RNU4-1* differs from *RNU4-2* at exactly four positions — n.37, n.88, n.99 and n.113 —
and **none lies inside the critical region**. In particular *RNU4-1* carries the *identical* four-thymine tract,
so conversion reproduces T₄ rather than extending it to T₅ and cannot generate n.64_65insT. Applying the same
test to every observed patient allele, **0 of 11 alleles and 0 of 114 patient carriers are conversion-explicable**.

### The excess is specific to one allele, and survives position-matched controls

Across the critical region, single-base insertions have 19.8 patient carriers per allele (5 alleles, 99 carriers)
against 2.5 for pathogenic substitutions (6 alleles, 15 carriers) — a 7.9-fold excess per allele. That
class-level comparison is, however, entirely attributable to one allele: **excluding n.64_65insT, the per-allele
insertion-to-substitution ratio is exactly 1.00** (2.5 versus 2.5) [58]. The excess is not a property of
insertions. It is a property of a single allele.

This permits a far stronger design than a class comparison. Because the effect is allele-specific, it can be
controlled at the level of the **nucleotide**: variants at the same position are covered by the same assays,
called by the same pipelines, reported in the same screens and ascertained together, so any explanation acting
on a region, a variant class or a sequencing difficulty is held constant by construction.

chr12:120,291,839 provides exactly this control set. Four variants are observed there:

| variant at chr12:120,291,839 | class | patient carriers | population carriers |
|---|---|---|---|
| **n.64_65insT** | insertion extending the T₄ tract | **89** | 1 |
| n.64_65insG | insertion at the same nucleotide | 2 | 0 |
| n.64_65insC | insertion at the same nucleotide | 0 | 2 |
| n.65A>G | substitution at the same nucleotide | 2 | 0 |

n.64_65insT accounts for **89 of the 93 patient carriers at this position**. The mutational supply measured
above — a T₄ insertion occurring at 0.40 times the frequency of a substitution — speaks directly to the contrast
between the recurrent insertion and the substitution at the same nucleotide, and there the recurrent allele
takes **89 of 91 carriers** (97.8%, 95% CI 92.3–99.7%) where 28.7% is expected. The departure is extreme:
exact binomial **p = 1.3 × 10⁻⁴⁵**; by Poisson, 26.1 carriers are expected and 89 observed
(P(X ≥ 89) = 5.2 × 10⁻²²). The **position-matched residual is 110-fold**.

The two other insertions at this nucleotide are the decisive control. n.64_65insG and n.64_65insC are the same
class of variant call, in the same homopolymer, in the same patients' data, and they carry 2 and 0 patients
between them. Insertion-calling behaviour in this tract therefore cannot generate the observed excess; nor can
assay coverage, since all four alleles occupy one base.

### The excess is exclusively maternal

The parental origin of this excess is reported in the primary literature and is, so far as we are aware,
without precedent in its consistency. In the discovery cohort, parental origin was resolved for 54 individuals —
46 carrying n.64_65insT, three other insertions, five substitutions — and **all 54 variants lay on the maternal
allele**. In the independent French genome-sequencing cohort, origin was determined in 50 trios and one
mother–patient duo: 47 maternal and four paternal, and **none of the four paternal variants was n.64_65insT**
(they were n.62T>C, n.68A>C, n.76del and n.92C>G). The excess quantified above is therefore not merely
associated with maternal transmission; within the resolution of the published data it occurs *only* on the
maternal allele.

---

## Discussion

The dominance of n.64_65insT in ReNU syndrome is not a consequence of an elevated mutation rate. The critical
region is less mutable than its surroundings; mutation rate does not predict which variants are seen; and the
specific sequence context of the recurrent allele is one in which insertions are disfavoured relative to
substitutions by roughly 2.5-fold. Nor is it a consequence of interlocus gene conversion, the usual explanation
for a recurrent allele with insufficient mutational supply: the only plausible donor carries the same tract. Two
of the candidate explanations can therefore be set aside.

**Ascertainment** is the more mundane of those remaining and is the reason the analysis is position-matched.
n.64_65insT was the allele through which the disorder was discovered and subsequent studies searched for it
specifically, which inflates its apparent share; a comparison across the gene, or across variant classes, cannot
separate that from a real excess. Restricting the comparison to a single nucleotide does. The three other
variants at chr12:120,291,839 are reported in the same screens, called from the same reads by the same
pipelines, and covered by any assay that covers the recurrent allele, yet carry four patients between them
against 89. An ascertainment account must therefore explain why the bias acts on one allele while sparing two
other insertions at the same base in the same homopolymer — including n.64_65insG, which differs only in the
identity of the inserted nucleotide. Two further observations point the same way: the allele's share is similar
in an independently ascertained cohort (72.6% versus 77.4%), and a single-base insertion in a homopolymer is the
harder call for short-read pipelines, not the easier one, so the technical bias runs against the observed
direction.

**That leaves germline selection, and this is where the maternal exclusivity becomes the central problem.**

The primary literature proposes positive selection in the female germline but does not quantify it. Our estimate
puts a number on what it would have to achieve: a 110-fold enrichment of this allele over the other variants at
the same nucleotide, beyond the mutational supply we measured. That magnitude is not
unprecedented in itself. Selfish spermatogonial selection — the process underlying the classical paternal
age-effect disorders, driven by activating variants in *FGFR2* (Apert), *FGFR3* (achondroplasia), *RET* (MEN2)
and other RAS–MAPK genes — produces enrichments of one to three orders of magnitude in offspring.

But that mechanism is a consequence of male germline architecture, and the female germline does not share it.
Spermatogonial stem cells divide continuously throughout adult life — several hundred divisions by middle age —
so a variant conferring even a slight proliferative advantage is amplified clonally across decades; this is why
mutant cells form discrete patches in aged testis and why the effect scales with paternal age. Oogenesis offers
a far smaller target. Oogonia do proliferate mitotically, so a proliferative compartment is not absent
altogether, but that proliferation is confined to roughly the first half of fetal development and comprises on
the order of thirty divisions, after which the oocyte pool enters meiotic arrest and the genome is not
replicated again. A selective process would therefore have to achieve, within a bounded prenatal window, an
enrichment that the male germline accumulates over hundreds of divisions and several decades — and would have
to do so reproducibly across unrelated women. No documented example of germline selection operates this way.

The maternal age effect that does exist does not fill the gap. It is substantially weaker than the paternal one
and is attributed to DNA damage accumulating in arrested oocytes and to double-strand-break-associated mutation
clusters — processes that raise mutation rate broadly rather than amplifying one allele over its immediate
neighbours. The position-matched comparison above bounds that directly: three other variants at the same
nucleotide, subject to any such process equally, carry four patients between them.

One alternative deserves explicit treatment. The authors of the French cohort suggest the maternal bias may
reflect negative selection against severely splicing-disruptive variants in the male germline. That would
explain why ReNU variants are maternal; it does not explain why *one* maternal allele is 44-fold more frequent
than the other maternal alleles at the same nucleotide, which would be purged equally. Directional selection
against paternal transmission changes the parent of origin, not the allelic spectrum within a parent.

We therefore state the problem rather than resolve it: ReNU syndrome exhibits a maternal-only germline excess of
roughly 110-fold that is specific to a single allele, at a site whose mutational supply we have measured and
found to be below that of a substitution, and that survives controls at the resolution of the individual
nucleotide. No described mechanism of germline mutation or selection produces that pattern.

Three measurements would discriminate among the possibilities. **(1) Maternal age.** If the excess arises from
damage during meiotic arrest, the n.64_65insT fraction among ReNU cases should rise with maternal age; if it
arises from a replication-independent process fixed before birth, it should not. Trio cohorts already hold the
data. **(2) Direct assay of the female germline.** Deep targeted sequencing of *RNU4-2* in ovarian tissue or
oocytes, as was done for *FGFR2* and *FGFR3* in testis, would show whether the allele is present above
expectation before fertilisation — distinguishing a germline process from a post-zygotic or transmission one.
**(3) Transmission distortion.** If the allele is favoured at fertilisation or in early embryogenesis rather
than generated more often, that is detectable in preimplantation or early-loss material, and would place the
mechanism outside mutation entirely.

### Limitations

Observed allele counts in population data reflect mutation rate, genetic drift and selection together, so our
measurement is of relative mutability, not a per-generation rate. Short-read sequencing detects insertions in
homopolymers less sensitively than substitutions; we measured that bias directly rather than assuming it,
finding a 1.73-fold under-detection at length four, and adopted the long-read-corrected rate. The HPRC estimate
rests on 232 assemblies, so its counts are sparse, and graph-based calls in repetitive sequence carry their own
uncertainty; the true rate most likely lies between the two estimates, and we have taken the one that minimises
our reported effect. Sampling was restricted to chromosome 12, and genome-wide sampling would tighten the
estimate. Patient carrier counts were transcribed from a published table and inherit its ascertainment.

The parental-origin data are not ours: they are reported in the two cohort papers and we have combined them
with our own measurement rather than re-deriving them. Origin was inferred indirectly from nearby informative
variants, and was resolvable in only a subset of cases, so the claim is that no informative counter-example has
been reported — not that none exists. The gene-conversion test excludes the proximal tandem donor; more distant
U4 pseudogenes were not examined, though homology and proximity both make them far less plausible donors.
Finally, we exclude mechanisms; we do not establish what the operative one is.

---

## Figures

**Figure 1.** Mutation rate across *RNU4-2*. (a) Roulette relative rate per nucleotide with the 18-nt critical
region shaded; patient variants overlaid. (b) Rate inside versus outside the critical region. (c) Rate versus
patient carriers and versus population carriers, both flat.

**Figure 2.** Homopolymer insertion frequency. (a) Relative insertion frequency against tract length, log scale,
with the fitted slope; T₄ marked. (b) The same split by base. (c) Short-read versus long-read estimate at each
tract length, showing the measured under-detection.

**Figure 3.** The position-matched residual and its parental origin. (a) All four variants observed at
chr12:120,291,839, with the binomial test of the recurrent insertion against the substitution at the same
nucleotide; the two other insertions at that base are the calling-bias control. (b) The *RNU4-2* / *RNU4-1* alignment over the critical region, showing the
identical T₄ tract and the four paralogue differences all lying outside it. (c) Parental origin of informative
ReNU cases in both published cohorts.

## Data files

| file | contents |
|---|---|
| `results/renu_recurrence_rates.tsv` | per-variant Roulette rates and observed carrier counts |
| `results/homopolymer_insertion_rate.tsv` | insertion frequency by tract base and length (short read) |
| `results/homopolymer_rate_lr_vs_sr.tsv` | long-read versus short-read rate at each tract length |
| `results/gene_conversion_test.tsv` | conversion-explicability of every observed patient allele |
| `data/processed/roulette_rnu4_2.tsv` | all 423 Roulette records for *RNU4-2* |
| `data/curation/renu_variant_counts.tsv` | transcribed per-variant carrier counts |

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
the number of callable bases [49]. The identical procedure, with the same seed and the same windows, was applied
to the HPRC Release 2 (v2.1) Minigraph-Cactus pangenome [50].

**Gene conversion test.** Mature transcribed-strand sequences for *RNU4-2* and *RNU4-1* were taken from the
GRCh38 reference, flanks trimmed and minus-strand genes reverse-complemented, then compared position by
position. An observed patient substitution is conversion-explicable only if the donor carries the alternate
base at the aligned position; an insertion is conversion-explicable only if the donor carries a longer
homopolymer run at that site [56].

**Parental origin** was taken from the published cohort reports and is not re-derived here.

**Patient carrier counts** were transcribed from the published Extended Data table of the discovery cohort and
are provided in `data/curation/renu_variant_counts.tsv`.

**Statistics.** Spearman correlation for rank association; linear regression on log₁₀ relative frequency against
tract length; exact Poisson intervals on insertion counts.

**Code and data availability.** All analyses are scripted and committed; each result names its script.
