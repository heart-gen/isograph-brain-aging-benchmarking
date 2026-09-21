# Claim calibration

Every rule here is a claim this project **made and then retracted or bounded**. They are
not stylistic preferences; quoting the legacy version puts a withdrawn number in front of
a reviewer. Sources: `reports/pi/00_OVERVIEW.md` §1a and §1b, and the stage report named
in each row.

---

## 1. The north star

> IsoGraph is a **complementary isoform-switch network method**, not a globally superior
> one. Its defensible contribution is a small, specific **DTU-without-DGE layer** —
> coordinated isoform-usage change without a corresponding change in total gene abundance
> — that is structurally invisible to any DGE or abundance-WGCNA pipeline.

The report is written to that bound throughout. Where the project's own analyses run
against its interest, they are reported.

---

## 2. Reversed on the 2026-09-19 re-run — never quote the legacy version

| Legacy claim (withdrawn) | What the report must say | Stage |
| --- | --- | --- |
| "The splicing-specificity contrast is IsoGraph-only; no matched baseline clears 0.05" | A **representation** effect, not an inference one: `wgcna_multiplex` on identical features shows it **more** strongly — 1.084 (p = 0.0044) vs IsoGraph 1.065 (p = 0.020). `wgcna_switch_only` is null and below 1 (0.927, p = 0.065). | 05 R2 |
| "One clean method effect: IsoGraph beats matched multiplex WGCNA on phenotype rate (0.268 vs 0.189)" | **Inverted in the pooled mean**: 0.146 vs 0.189. Cohort-dependent — IsoGraph leads in BrainSEQ (0.176 vs 0.163), trails in GTEx. | 04 R1 |
| "23 of 130 module pairs replicate for aging across cohorts (P = 0.034)" | **Null and retired**: 1/38 (perm P = 0.91); WGCNA 2/52 (P = 0.15). BrainSEQ is quantified with Salmon and GTEx with RSEM, so the arm confounds processing with biology. | 03 R2 |
| "Resolution 5.0 removes the giant-module GWAS artifact (0/10 vs 17/28)" | **Reversed under the switching filter**: res 5.0 is 20/26 (77%) giant against res 2.0's 24/37 (65%). Production moved to Leiden 2.0 on 2026-09-16. | 02 R3 |

---

## 3. Demoted headlines

### SNCA is a falsification example, not a worked example

PI decision, 2026-09-18. SNCA is the project's most attractive biology and is **not**
evidence the method was needed:

- identical coloc in pools that do not depend on IsoGraph — PP4_sQTL **0.975** in the
  IsoGraph switch pool, in the all-testable-genes background pool, and in both matched
  WGCNA pools;
- HEIDI rejects the SNCA–PD shared-signal model in **all 8** instrumented tissues;
- BrainSEQ's switch axis does not colocalize;
- it fails long-read abundance qualification, and on the re-run it is not a concordant
  gene, so the junction arm no longer tests it.

Write it as what it is: a prominent disease gene that shows transcript-usage structure
without meeting the genetic evidence a causal switch interpretation needs. Do **not** draw
it as resolved in Fig 4 panel A (that panel is PRDM2).

### Splicing specificity is set-level support, not the headline

Headline decision, 2026-09-12: the lead is the **switch layer itself** — reproducible,
largely abundance-independent, orthogonally confirmed. Splicing specificity holds at set
level only and must be reported with the per-gene results beside it, which run the other
way: GTEx signal-level **22 vs 45** (P = 0.007), BrainSEQ in-sample **5 vs 33**
(P = 4.3e-6); per gene GTEx **62 splicing-only vs 161 expression-only** (P = 2.5e-11),
BrainSEQ **16 switch-only vs 94 abundance-only** (P = 1.3e-14). Two of three pre-specified
sensitivity arms no longer clear 0.05.

Also withdrawn: the effect does **not** localize to GO-invisible modules (GO-invisible
1.080, p = 0.093; GO-visible 1.049, p = 0.126 — both null; only the phenotype-associated
set clears 0.05).

### The RBP arm is a secondary annotation layer

PI decision, 2026-09-18. No sentence may lean on it for method validation — the split-half
and cross-cohort evidence carries that argument. The eCLIP data are **HepG2/K562**, not
brain; state the tissue mismatch explicitly and limit the claim to *candidate regulatory
compatibility*, never direct evidence of binding in brain.

### The allelic arm: what it does and does not license

- **KLC1 × SCZ is not "directionally inconsistent."** Its four significant rows are one cis
  effect on KLC1-213 measured through four partners; the apparent sign disagreement was an
  artifact of projecting onto a module axis the gene does not follow (06a §5a). Name it as a
  caveat pending a transcript-level estimand, not as a schizophrenia switch locus.
- **PICALM-219 × AD is a Result** (PI, 2026-09-20), stated with its prior-sensitive coloc in
  the same sentence. It remains `disease_locus_splice_linked` — never "a recovered known
  mechanism".
- **Cis control is gene-level.** It is not concentrated in the aging modules (06a §5c); no
  sentence may imply the aging-relevant modules are the cis-regulated ones.
- **No general disease-direction effect.** At the GWAS lead, directions are at chance
  (214/371). Do not write that risk alleles push isoform switching in a direction.

### "GO-invisible" is not the central biological payoff

Report it as a property of the layer, not as the result.

---

## 4. Framing rules

### Gene-wise DTU is complementary, never inadequate

Never write "why not per-gene DTU" or imply satuRn/DEXSeq are insufficient — the project
has not shown that and does not need to. satuRn asks which individual genes show
differential transcript usage; IsoGraph asks whether those changes are organized into
coordinated programs. The held-out test is honest about its size: module context adds
modestly, in **3 of 6** region/cohort analyses after conditioning on co-expression
(partial r 0.056–0.102), and not in the others. Heterogeneity is a result, not a failure.

### Validation is not replication

Within-cohort split halves **validate**. A second cohort **replicates**. The cross-cohort
module-pair arm is null and quantifier-confounded; the surviving cross-cohort result is
the **eigengene projection** (33/38 and 44/55, P ≈ 1e-6, where matched WGCNA is null at
28/52 and 3/10).

### Features versus inference

A result the matched WGCNA baselines reproduce on identical switch features belongs to the
**representation**, not to IsoGraph's network inference. §1a's table is the authority.
Currently inference-specific: cross-cohort eigengene projection, and synthetic module
recovery where switching dominates.

### Do not claim global superiority

WGCNA's chance-trusted rate (**62/72, 86%**) is *above* IsoGraph's (**93/118, 79%**).
Report both percentages wherever the trust funnel appears.

---

## 5. Retracted display item

**Fig 4E's SCZ coloc-convergence panel is retracted** — its all-genes denominator was
ascertainment (the published P = 0.0058 does not survive a CLPP-tested denominator). The
honest replacement exists: colocalizing genes do not concentrate in modules in any trait
(permutation P = 0.19–1.00, null in all 10 cells). Never re-quote the retracted version.

---

## 6. Generation guard

- Do not quote **pre-2026-08-29** QTL values (the `qtl_anchoring` refresh changed the
  Fig 3 headline).
- Do not quote **legacy expression-filter** numbers. Production has used the switching
  transcript filter since 2026-09-14 and Leiden resolution **2.0** since 2026-09-16.
- Quote S-LDSC on **`coef_p`**, the coefficient test conditional on baselineLD — not
  `enrichment_p` (discrepancy register D-2026-09-10-a).
- Treat any `*_SUMMARY.md` as hand-written prose until proven otherwise; four were found
  claiming to be machine-generated and had never moved with their analyses.
- A module-id join is valid only against the fit that wrote it. `module_genetic_anchoring`
  once pooled a stale partition and 26 of 80 rows came from a legacy fit.

---

## 7. Grep pass

Run against the finished report and both appendices. Any hit needs a reason.

```bash
grep -n -i -E \
  "outperform|superior to WGCNA|better than WGCNA|IsoGraph-only|\
23 of 130|23/130|0\.268|17/28|0/10 sig|\
worked example|why not per-gene|satuRn is|inadequate|\
genetically anchored|replicate[sd]? across cohorts|cross-cohort replication of modules|\
proves|demonstrates that .* causes|causal(ly)? drives|\
directionally inconsistent|disagrees with itself|known mechanism recovered" \
  results_report.md appendix_genetics.md appendix_mechanism.md
```

Then confirm by hand:

- every reversed claim appears only in its corrected form;
- "replication" is used only for a second cohort;
- the trust funnel carries both percentages;
- SNCA appears only as a falsification example;
- the RBP arm carries the HepG2/K562 caveat;
- partial coverage is stated wherever a cited analysis is not fully `run`.
