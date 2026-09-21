---
name: analysis-results-report
description: >
  Produce the IsoGraph brain-aging project's results report — a
  biological-question-centered synthesis for collaborators, lab meetings and
  manuscript development — from the completed PI review in reports/pi/ and the
  stage result files. Project-scoped: supersedes the generic
  analysis-results-report skill for this repository, because this project's
  claim calibration (four reversed conclusions, a retracted panel, a demoted
  headline) is not derivable from the outputs alone. Emits a Markdown report,
  two appendices and a LaTeX Beamer briefing.
---

# Analysis results report — IsoGraph brain aging

The `analysis-to-pi` stage is complete: `reports/pi/` holds ten stage reports, an
evidence inventory over 83 analyses, and a discrepancy register. This skill is the
**scientific reporting stage** that turns that audit into a report a collaborator can
read in two minutes and a co-author can write Results from.

Do not redo the audit. Reinspect code or outputs only to verify a number, resolve a
discrepancy, or identify the right figure.

---

## 1. Read in this order

Order matters: the first two set the ceiling on every claim in the report.

1. **`reports/pi/00_OVERVIEW.md` §1a** — what demonstrates IsoGraph's utility, the
   adopted narrative, and the features-vs-inference table that says which results need
   the *method* and which only need the *switch features*. Most over-claims this report
   could make are already refuted here.
2. **`reports/pi/00_OVERVIEW.md` §1b** — the four conclusions that reversed on the
   2026-09-19 switching-filter re-run. Legacy versions of these are still quoted in older
   prose across both repos; they must never enter the report.
3. **The stage reports** — `reports/pi/0N_*.md`, including `05a_signal_level_genetics.md`
   and `06a_allelic_switch_arm.md`. Every one has the same shape:
   - §5 Main results — the narrative and the caveats
   - §6 Key quantitative results — the numbers table, with its regeneration date
   - §7 Robustness, negative results, and limitations — the counter-evidence
   - §12 Reproducibility and audit trail — **the file paths to verify against**
4. **`ANALYSIS_MAP.md`** — analysis → CLI → wrapper → output → display item.
5. **`reports/pi/_evidence/inventory.md`** — the status board (`run`, `partially_run`,
   `partial_coverage`, `unresolved`, `not_run`) and last content-change commit.
6. **`reports/pi/_evidence/DISCREPANCY_REGISTER.md`** — numbers a repository document
   gets wrong about its own output.
7. **`manuscript/FIGURE_ORDERING.md`** — the display-item set and the one honest claim
   each figure carries. Figure legends there are already claim-calibrated; reuse their
   wording.
8. **Manuscript context** — `manuscript/MANUSCRIPT_PLAN.md`, and in the manuscript repo
   (`../../manuscript/isograph-brain-manuscript/`) `TODO.md` and `content/`.

`references/evidence-registry.md` maps each question below onto these sources.

---

## 2. Verification rule

This project's convention, inherited from the PI review: **every number in the report is
re-read from the result file that produced it, not from prose.** Locate the file via the
stage report's §12 audit trail, then read it.

Before a number goes in:

- **Status.** Check `_evidence/inventory.md`. Never quote a `partially_run` or
  `partial_coverage` analysis without stating the coverage in the text.
- **Register.** Check the discrepancy register. It records, among others, that S-LDSC must
  be quoted on `coef_p` (coefficient conditional on baselineLD), not `enrichment_p`.
- **Generation.** Do not quote pre-2026-08-29 QTL values, and do not quote
  legacy-expression-filter numbers — production has been on the switching transcript
  filter at Leiden resolution 2.0 since 2026-09-14/16. When a stage report's §6 header
  gives a regeneration date, that is the generation the report speaks from.
- **Partition provenance.** A module-id join is only valid against the fit that wrote it;
  Leiden ids are re-assigned every fit. Never combine a module id across generations.

Heavy compute does not belong in this skill. Reading parquet result files on the login
node is fine; refitting anything is not.

---

## 3. The six biological questions

The report is organized by these, not by stage order. They are fixed so the report is
re-runnable and comparable across revisions. Sources for each are in
`references/evidence-registry.md`.

1. **Does the method recover co-switching modules where ground truth is known?**
   Synthetic benchmark. The honest bound: IsoGraph wins where switching dominates, not
   everywhere.
2. **Is coordinated isoform switching a measurable genome-wide layer in the aging human
   brain, and is it separable from abundance?** The resource — switch modules across
   aging analyses in two cohorts — plus the switch-unique partial-R² result that makes
   "separable from abundance" quantitative and threshold-free.
3. **Are the modules trustworthy and reproducible?** Per-module trust funnel, split-half
   stability and driver reproducibility, the curvature test that licenses the linear age
   model, and cross-cohort eigengene projection — **the strongest surviving
   inference-specific result**. The cross-cohort module-*pair* count is null and retired;
   report it as a negative.
4. **Are the switches real molecular events?** ONT long-read confirmation, the
   allele-aware junction recount, and the switch-consequence analysis (productive UTR/CDS
   remodeling, not decay).
5. **Is the layer an artifact of cell-type composition?** Composition adjustment — which
   survives in the limbic/striatal aging arm and **collapses in the two GTEx cortical
   regions**. Confounder vs mediator is unresolvable from these data; say so.
6. **Does disease and aging genetics act through the switch layer?** S-LDSC partitioned
   heritability, colocalization (abf, signal-level SuSiE, BrainSEQ in-sample), SMR/HEIDI,
   and the set-level splicing-specificity contrast — reported **with** the per-gene
   results that run the other way, in the same section, not in a later caveat.

RBP regulation is a **secondary annotation layer** inside question 6's interpretation. It
never carries a validation claim.

---

## 4. Claim calibration

Read `references/claim-calibration.md` before writing, and run its grep pass after. The
project's north star, which the whole report is written to:

> IsoGraph is a **complementary isoform-switch network method**, not a globally superior
> one. Its defensible contribution is a small, specific **DTU-without-DGE layer** that is
> structurally invisible to any DGE or abundance-WGCNA pipeline.

Three habits that keep the report honest:

- **Put the counter-evidence in the same section as the claim.** This project's stage
  reports do it and the report should too: WGCNA's trusted rate beside IsoGraph's, the
  per-gene QTL counts beside the set-level ratio, the cortical composition collapse beside
  the limbic survival.
- **Distinguish features from inference.** §1a's table already does this. A result that
  the matched WGCNA baselines reproduce on identical features is a *representation*
  result, not a method result.
- **Validation is not replication.** Within-cohort split halves validate; a second cohort
  replicates. The cross-cohort arm is quantifier-confounded (BrainSEQ Salmon, GTEx RSEM)
  and is not a replication claim.

---

## 5. Report structure

Portable Markdown — standard headings, short paragraphs, simple tables, relative figure
paths. It must survive conversion to Word and PDF. It should read like a scientific
results document, not a chat answer.

```
Title + metadata (project, cohorts, status, date, source = the PI review)
1. Executive scientific summary     biological problem, central question, 3–6 findings
                                    with numbers, emerging model, take-home
2. Study context and design         cohorts, regions, assays, n, comparisons — only what
                                    is needed to read the results
3. Biological questions and results  the six questions; per question:
                                       Hypothesis (or "Exploratory question")
                                       Rationale
                                       Analytical approach
                                       Findings (lead with the result)
                                       Interpretation (observed vs inferred)
                                       Evidence strength: Strong / Moderate / Preliminary
                                       Supporting figures and tables
                                       Manuscript-ready result statement
4. Integrated biological model      what converges, what is independent, what the layer
                                    contributes beyond a list of associations
5. Relationship to prior knowledge  reproduced / extended / novel / contradicted
6. Robustness and supporting evidence   organized by the conclusion supported, each marked
                                        reinforces | qualifies | weakens | fails
7. Figures and presentation story   per figure: question, main message, evidence,
                                    recommended role, changes needed
8. Manuscript development           central claim, Results subsection structure,
                                    candidate language, figure-to-claim table
9. Collaborator brief               what we asked / found / it means / need
10. Scientific implications and next steps   immediate | manuscript-strengthening | future
11. Limitations                     LAST. Only limitations that change interpretation;
                                    for each: what it is, which conclusion it affects,
                                    whether it weakens, narrows or merely qualifies.
```

Writing: direct, quantitative, active voice, biologically explicit. Lead with the result.
No "Interestingly", "It is important to note", "As can be seen". Do not narrate the
project's debugging history — the re-runs and reversals appear only where they change
what a reader should believe today.

### Appendices

Two, for the questions with more material than the main report can hold:

- **`appendix_genetics.md`** — stages 05, 05a, 06a: S-LDSC, coloc.abf, signal-level SuSiE
  coloc and its MAX_SNPS sensitivity, BrainSEQ in-sample coloc, SMR/HEIDI, the swQTL arms,
  the allele-aware junction recount, and the full per-gene modality contrast.
- **`appendix_mechanism.md`** — stages 06, 07, 08: long-read confirmation, junction
  validation, switch consequence, RBP regulons and their eCLIP tissue caveat, and the
  cross-layer agreement tables where the layers disagree.

Each appendix keeps the same evidence discipline and cross-references the main report's
question number.

---

## 6. Beamer briefing

`results_briefing.tex`, 12–14 content slides, built with `pdflatex` (no pandoc on this
machine). Reuse the preamble, palette and `\fig` / `\srcnote` macros from
`reports/pi/pi_briefing.tex`, with `\fig` repointed at the local `figures/` directory so
the deck compiles from either output location.

Follow the scientific argument, not the workflow: problem → question → design → the
questions in order → integrated model → robustness → implications → **limitations last**.
Conclusion-oriented slide titles ("Aging switch modules transfer across cohorts as a
projected eigengene", not "Replication analysis"). One message per slide; one figure;
one to three numbers; one interpretation line. Use LaTeX `%` comments for speaker notes —
supporting statistics, the anticipated objection, and the source file.

16:9, large sans-serif, high contrast, left-aligned, never meaning through color alone.

---

## 7. Output

Identical tree in both locations; the manuscript-repo copy is a verbatim sync.

```
reports/results/                →  ../../manuscript/isograph-brain-manuscript/drafts/results-report/
  results_report.md
  appendix_genetics.md
  appendix_mechanism.md
  results_briefing.tex
  results_briefing.pdf
  figures/                      staged PNGs, referenced relatively
```

Mechanical steps are a committed CLI — `python -m isograph_benchmark.reporting.results_report`
with `figures`, `build` and `sync`. The prose is hand-written; the staging, the LaTeX build
and the copy are not.

---

## 8. Repository constraints

From the project working plan (`reports/markdowns/AGENTS.md`), binding here:

- Login node for light work only. This skill reads and writes text; it runs nothing heavy.
- Never `git add -A`. Stage selectively. Commit only when explicitly asked.
- **No `Co-Authored-By` trailer** on commits in this repository.
- Do not push without explicit confirmation.
- Never hardcode a repository path. `isograph_benchmark/paths.py` holds `OUTPUT_DIRS`;
  address stages through it.
- Manuscript prose belongs in the manuscript repo, not in the analysis repo.

---

## 9. Before finalizing

- [ ] Every number traced to a result file, not to prose.
- [ ] No banned phrasing from `references/claim-calibration.md` survives the grep pass.
- [ ] Each reversed claim appears only in its corrected form.
- [ ] Partial coverage stated wherever a cited analysis is not fully `run`.
- [ ] Counter-evidence sits beside its claim, not in a later caveat.
- [ ] Sample sizes, effect directions and thresholds match the figures they cite.
- [ ] "Replication" used only for a second cohort.
- [ ] Association never written as causation.
- [ ] Limitations is the last substantive section, and the last content slide.
