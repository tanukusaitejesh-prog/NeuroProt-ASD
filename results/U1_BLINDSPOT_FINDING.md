# Two U1 snRNA genes were screened for NDD and published as negative; the negative is an artefact

## The claim
Chen et al. (Nature 2024, PMC11338827) Extended Data Table 4 screened 28 brain-expressed snRNA genes for de novo
and biallelic variant enrichment in 8,841 undiagnosed NDD probands in Genomics England. **RNU1-3 and RNU1-4
returned zero variants of every class** (de novo: 0 non-NDD / 0 NDD; homozygous and compound heterozygous:
0 / 0), so no Fisher's test was computable. In the same table RNU1-2 - the same 164 bp gene - yielded 18 and 18
de novo variants, and RNU4-2 yielded 39 (P = 2.48e-11), the discovery that defined ReNU syndrome.

Those zeros are not evidence of invariance.

## Evidence
gnomAD v4.1 (76,215 genomes, short read): RNU1-3 and RNU1-4 have **AN = 0** - no individual has a callable
genotype - while RNU1-1 and RNU1-2 sit at full cohort depth (AN ~152,000, ~1,155 PASS variants each).
(results/snrna_callability.tsv)

HPRC Release 2 v2.1 pangenome (232 long-read haplotype-assembled samples, ~464 haplotypes), which does not depend
on short-read mappability: (results/snrna_pangenome_blindspots.tsv)

| gene | gnomAD PASS | gnomAD AN | HPRC variants | HPRC alt haplotypes | var/kb |
|---|---|---|---|---|---|
| RNU1-1 | 1,154 | 152,061 | 10 | 258 / 464 | 61 |
| RNU1-2 | 1,157 | 152,134 | 22 |  40 / 464 | 134 |
| **RNU1-3** | **0** | **0** | **22** | **146 / 464** | **134** |
| **RNU1-4** | **0** | **0** | **15** | **62 / 464** | **91** |

RNU1-3 is exactly as variable as RNU1-2 (134 var/kb in both) and 31% of assembled haplotypes carry a
non-reference allele there.

## Why it matters
- U1 is the only major spliceosomal snRNA family with no established NDD gene. U2, U4, U5, U6, U11, U12,
  U4atac and U6atac all have one.
- The primary evidence used to nominate snRNA disease genes is population depletion plus de novo enrichment.
  Neither can be computed where AN = 0.
- A published systematic screen therefore contains silent blind spots, and two loci remain genuinely unassessed
  while appearing in the literature as tested and negative.
- The 2026 Nat Genet systematic screen (s41588-026-02547-5) independently reports "recurrent low VAFs in RNU1-1,
  RNU1-2, RNU2-1 ... consistent with their high sequence identity", i.e. the same mapping problem, attributed to
  artefact and set aside.

## Honest limits
- 232 HPRC samples establishes that the loci are variable; it does not give allele frequencies and says nothing
  about pathogenicity.
- Pangenome graphs can themselves err in repetitive regions, so some alt-haplotype calls may be graph artefacts.
- RNU1-1 shows 1 de novo call in GEL vs 36 for RNU1-2 despite both being fully callable in gnomAD, which suggests
  the GEL de novo pipeline may be additionally compromised at U1 beyond what gnomAD AN reveals. Worth checking.
- Chen's dismissal of the RNU1-2 / RNVU1-7 recurrent de novo variants was sound and is NOT in question: those
  variants are common (AF > 0.5% in gnomAD v4.0). This finding concerns RNU1-3 and RNU1-4 only.
- Whether RNU1-3 / RNU1-4 are transcribed at levels that would make them disease-relevant is not established here;
  Chen included them among "brain-expressed" genes, which is supportive but not decisive.

## What would close it
Long-read NDD cohort data over chr1:16,666,785-16,666,948 (RNU1-3) and chr1:16,740,516-16,740,679 (RNU1-4).
