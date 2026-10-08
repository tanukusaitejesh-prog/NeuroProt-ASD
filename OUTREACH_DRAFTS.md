# Outreach drafts — review, edit in your own voice, then send

These are the moves that change what is possible. Each is short on purpose: busy people read three paragraphs.

---

## 1. To the RNU4-2 saturation genome editing group (highest value)

**Who:** the senior authors of De Jonghe et al., *Nature* 2026, "Saturation editing of RNU4-2 reveals distinct
dominant and recessive disorders" (Crick / Whiffin / Findlay orbit — check the corresponding author line).

**Why they should care:** they have screens running on RNU2-2 and RNU5B-1. You have a model that already
predicted their RNU4-2 data without ever seeing it, plus timestamped predictions for the two genes they are
currently screening. That is a prospective test they can run at zero cost and publish either way.

> **Subject:** Structural model predicts your RNU4-2 SGE scores zero-shot — timestamped predictions for RNU2-2 and RNU5B-1
>
> Dear Dr [name],
>
> I've been working on a structure-based variant effect model for spliceosomal snRNA genes, built from contact
> profiles across 11 cryo-EM structures spanning the splicing cycle. With RNU4-2 and the entire U4 family
> excluded from training, it recovers your SGE function scores at Spearman ρ = 0.41 (CADD 0.07, phyloP −0.00),
> and reaches ClinGen PP3_moderate when calibrated against your functional classes.
>
> I understand saturation screens of RNU2-2 and RNU5B-1 are underway. I have deposited timestamped, hashed
> predictions for every possible SNV in both genes, with the evaluation protocol fixed in advance — metric,
> comparators, success criterion, and a commitment to report the result whichever way it falls. If those screens
> are completed, the comparison would be a genuinely prospective test of structure-based prediction in
> non-coding RNA, which as far as I know does not yet exist.
>
> I would be glad to share the predictions and the code now, under whatever arrangement suits you, and equally
> glad for you to run the comparison yourselves. A preprint describing the model is in preparation and I would
> welcome your comments on it before it goes up.
>
> With thanks,
> [name, affiliation]

---

## 2. To the Nava / Depienne group (RNA-seq and mechanism)

**Why:** they hold the patient RNA-seq (EGAS50000000889) and described the 5′SS usage defects. Your structural
model identifies which variants perturb the activation machinery; theirs measures what actually happens to
splicing. The combination answers a question neither side can answer alone.

> **Subject:** Multi-state structural model of snRNA variants — possible complement to your 5′SS usage data
>
> Dear Dr [name],
>
> Your 2025 Nature Genetics paper showed that RNU4-2 variants alter 5′ splice-site usage in a way that tracks
> variant location and severity. I have built a model that represents each snRNA nucleotide by its contacts
> across 11 cryo-EM structures of the splicing cycle; it independently recovers the RBM42-interacting region as
> the strongest dominant-variant signal (OR 12.0, FDR 0.003) without being given your domain assignments.
>
> The obvious next question is whether the structural state a variant perturbs predicts *which* 5′ splice sites
> fail. I cannot answer it — the RNA-seq is under controlled access and rightly so — but if it is of interest I
> would be glad to provide per-variant state-resolved predictions for your cohort and have you test them, or to
> apply for access through the proper committee with your support.
>
> [one line on who you are]
>
> With thanks,
> [name]

---

## 3. EGA data access — EGAS50000000889

**Before applying, have ready:** a one-paragraph research plan, your institutional affiliation and a signing
officer, and a data management statement. If you have no institutional affiliation this will not succeed on its
own — in that case route 2 (collaboration) is the realistic path, and should be sent first.

**Research plan paragraph:**

> We have developed a structure-based variant effect model for spliceosomal snRNA genes, using per-nucleotide
> contact profiles across 11 cryo-EM structures spanning the splicing cycle. The model predicts saturation
> genome editing function scores for RNU4-2 zero-shot (ρ = 0.41) where conservation and CADD do not, and
> identifies dominant variants as perturbing activation and catalytic machinery contacts. We wish to test
> whether the structural state a variant perturbs predicts which 5′ splice sites are mis-selected in patients
> carrying that variant. We request access to the RNA-seq data in EGAS50000000889 to compare per-variant,
> state-resolved structural predictions against measured differential 5′ splice-site usage. No attempt will be
> made to re-identify participants; analysis is restricted to aggregate splicing measurements.

---

## 4. n-Lorem / BioMarin (ReNU therapeutics) — lower priority, but cheap

They are developing an ASO against the recurrent n.64_65insT allele, which covers ~75% of ReNU cases. Your
resource scores every possible variant in the gene, including the remaining 25% and the recessive disorder. A
short note offering the scored variant set costs nothing and occasionally opens doors.

---

## Order and timing

1. **Email 1 first.** It is the only one with a deadline — its value depends on the screens not yet being public.
2. **Email 2 next**, regardless of whether 1 is answered.
3. **EGA application** only if you have an institutional signatory.
4. Post the preprint before or alongside these, so there is something citable to point at.

## Before any of them: deposit to Zenodo

`results/prospective_deposit_2026-10-09.tsv` + `.json`, as a single versioned record. Ten minutes, and it turns
"I have predictions" into "here is a DOI and a timestamp that predates your data."
