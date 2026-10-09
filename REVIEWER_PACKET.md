# What to give a reviewer, and what to ask them to check

For the ReNU recurrence manuscript. The point of handing someone this list is to get the paper attacked where it
is actually weak, rather than praised where it is already solid. The soft spots are named below deliberately.

---

## 1. What to send

| item | file |
|---|---|
| the manuscript | `MANUSCRIPT_RECURRENCE.pdf` (9 pp, figures embedded) |
| source, if they want to diff | `MANUSCRIPT_RECURRENCE.md` |
| the one-command check | `SNRNA_LABELS=labels_v5c.tsv python scripts/59_verify_recurrence.py` → should print **73 of 73 automated checks PASS** |
| the carrier counts everything rests on | `data/curation/renu_variant_counts.tsv` |
| the analysis scripts named in the text | `scripts/48, 49, 50, 56, 58, 59, 60` |
| supporting result tables | `results/renu_recurrence_rates.tsv`, `homopolymer_insertion_rate.tsv`, `homopolymer_rate_lr_vs_sr.tsv`, `gene_conversion_test.tsv`, `position_matched.tsv`, `position_matched_sensitivity.tsv`, `paralogue_control.tsv`, `selection_requirement.tsv` |

Tell them the repository is self-contained except for remote queries to gnomAD, HPRC and Roulette, which are
made by HTTP range request and need no downloads.

---

## 2. The four claims that carry the paper

Ask them to attack these in order. If any one fails, the paper changes.

**C1. The mutational supply at a T₄ tract is below that of a substitution (0.40×).**
Everything downstream is a ratio against this number. Check: the derivation in `scripts/49` (short read) and
`scripts/50` (long read); whether relative-frequency-per-tract versus per-possible-substitution is the right
normalisation; whether 10 windows of chr12 is enough; whether adopting the long-read figure is genuinely the
conservative choice (it is the one that *shrinks* the reported effect).

**C2. The excess is allele-specific and survives position matching.**
`scripts/58`. Check: that comparing within one nucleotide really does hold ascertainment constant; that the
binomial is the right test; that the leave-one-out (ratio exactly 1.00 without n.64_65insT) is correctly
computed and correctly interpreted.

**C3. Gene conversion from *RNU4-1* cannot produce the allele.**
`scripts/56`. Check: the coordinate arithmetic for the minus-strand transcript numbering; that the flank
trimming and reverse complement are right; that "donor must carry a longer run" is the correct criterion.

**C4. The mechanism argument.**
Check the two numbers it rests on: ~30 mitotic divisions in oogenesis, and ~610 spermatogonial divisions by
paternal age 40 with up to 1000-fold enrichment. These come from the literature, not from us. The conclusion
should be tested against their preferred values — a sensitivity table is in `scripts/60`.

---

## 3. Known soft spots — ask them to confirm or refute

These are the places I would attack, and they are already stated in the manuscript. A reviewer who finds
something here is finding something we already flagged, which is the point.

1. **The carrier counts were transcribed from a published table image** (`chen2024_extdata_table4.jpg`).
   Transcription error is possible and would propagate everywhere. Ask them to re-check a few rows against the
   source table independently.
2. **Ascertainment is not fully excludable.** n.64_65insT was the discovery allele. The position-matched design
   narrows the objection — a bias would have to act on one allele while sparing n.64_65insG at the same base —
   but it does not eliminate it.
3. **The residual depends on the comparator.** 15× using every substitution in the critical region, 110× at the
   focal nucleotide alone. The manuscript now leads with the conservative 15× and tabulates the range; ask
   whether they agree that is the honest presentation.
4. **"No documented mechanism" is an argument from absence.** It is bounded quantitatively (8–15× the documented
   male per-division advantage) but it remains an absence claim.
5. **Parental origin is not our data.** Both cohort papers report it; origin was inferred indirectly from nearby
   informative variants and resolved in only a subset. The claim is that no informative counter-example has been
   published, not that none exists.
6. **Small absolute numbers.** 89 carriers of the focal allele, 2 of the comparator substitution, 4 alleles at
   the position. Ask whether the exact tests are appropriate at these counts.
7. **The gene-conversion test covers only the proximal tandem donor.** More distant U4 pseudogenes were not
   examined, on the grounds that homology and proximity both make them implausible donors. Ask whether that is
   accepted.
8. **chr12 only.** The homopolymer sampling is one chromosome; genome-wide sampling would tighten it.
9. **HPRC is 232 assemblies**, so the long-read counts are sparse and graph genotypes in repetitive sequence
   carry their own error.

---

## 4. Things a reviewer must check against the literature, not the repository

`scripts/59` cannot verify these, and says so when it runs:

- 77.4% share in the discovery cohort; 72.6% in the French cohort
- parental origin: 54/54 maternal (discovery); 47 maternal / 4 paternal of 51 (French)
- that the four paternal variants are n.62T>C, n.68A>C, n.76del, n.92C>G, and that none is n.64_65insT
- that saturation editing found n.64_65insT not unusually damaging
- the cohort sizes in the header of `renu_variant_counts.tsv`
- the spermatogonial and oogonial division counts underpinning C4

---

## 5. Three questions worth asking directly

1. Is the conservative 15-fold residual the right headline, or should the paper lead with the position-matched
   110-fold and relegate the range?
2. Is "we exclude two mechanisms and bound the third" enough for a short report, or does a reviewer expect the
   mechanism to be identified?
3. Does the paralogue control (*RNU4-1* carrying the equivalent insertion at AC 5) read as reassuring — the
   allele is generated at this context — or as undermining, since it shows the insertion is not rare in absolute
   terms? We read it the first way and say why; an independent reader should test that.

---

## 6. What not to spend their time on

- Arithmetic. `scripts/59` checks 73 numbers against source; if it passes, the numbers match the data.
- Whether the class-level insertion-versus-substitution comparison is sound. It is not, we found that
  ourselves, and the paper says so explicitly before using the position-matched design instead.
