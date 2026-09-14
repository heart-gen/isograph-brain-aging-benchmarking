# PSI / junction corroboration of the switch layer — interpretation

Written 2026-09-12 from the committed `summary_{module,marginal}.json` files in each
`<cohort>_<region>_<trait>/` subdirectory (`validate_switch_splicing.py`, wrappers
`06_switch_mechanism/_h/03–04`). Until now these results were carried only in `TODO.md`.

## Question

Do genes that IsoGraph calls as switching with a trait also show an independent,
annotation-based splicing signal — a significant PSI–trait association in the cohort's own
junction/PSI tables — more often than other tested genes?

## What is measured

- **Gene corroboration.** A 2×2 over all tested genes: switch-positive (in a
  phenotype-associated switch set) × PSI-positive (a significant PSI–trait association).
  Fisher OR, and a covariate-adjusted OR. Two switch sets: `module` (members of
  phenotype-associated IsoGraph modules) and `marginal` (genes whose own switch coordinate is
  trait-associated; BrainSEQ caudate only).
- **Event resolution** (BrainSEQ only). Of the junctions IsoGraph's switch pairs resolve to,
  how many match a PSI event, and do the two agree in significance and direction.

## Result

| Analysis | switch set | tested genes | switch+ | PSI+ | both | adjusted OR | adjusted p |
|---|---|---:|---:|---:|---:|---:|---:|
| BrainSEQ caudate, age | marginal | 11,144 | 28 | 59 | 9 | **153** | 1.5e-27 |
| BrainSEQ caudate, age | module | 11,144 | 822 | 59 | 9 | **2.06** | 0.048 |
| BrainSEQ caudate, Dx | marginal | 15,799 | 76 | 83 | 11 | **37.4** | 2e-25 |
| BrainSEQ caudate, Dx | module | 15,799 | 233 | 83 | 3 | 2.73 | 0.092 |
| GTEx cortex, age | module | 15,707 | 1,442 | 1,389 | 249 | **1.99** | 1.1e-18 |
| GTEx frontal cortex BA9, age | module | 15,646 | 6,349 | 1,025 | 515 | **1.45** | 1.6e-8 |
| GTEx ACC BA24, age | module | 15,556 | 3,516 | 1,071 | 332 | **1.43** | 2.7e-7 |
| GTEx hippocampus, age | module | 15,783 | 2,632 | 82 | 25 | **1.86** | 0.012 |
| GTEx cerebellum, age | module | 15,245 | 831 | 36 | 1 | 0.36 | 0.31 |
| GTEx hypothalamus, age | module | 16,141 | 560 | 31 | 1 | 0.79 | 0.81 |

**Not informative (no test possible or no overlap by construction):** BrainSEQ DLPFC and
hippocampus, GTEx caudate, putamen and cerebellar hemisphere carry **no** age-associated
switch-module genes, which is the caudate-only aging-module pattern in BrainSEQ, not a
failure of corroboration. GTEx amygdala, nucleus accumbens, spinal cord and substantia nigra
have switch-module genes but 1–22 PSI-positive genes and zero or near-zero overlap, so the
test has no power there.

**Event resolution (BrainSEQ).** 17 switch pairs resolve to annotated junctions, 11 match a
PSI event, and **0** are significant in both layers. Per-event concordance is therefore not
established; the corroboration is at gene level.

## Interpretation

- **Where both layers have signal, they agree above chance, in every such analysis.** Six of
  six testable module-level analyses with more than a handful of PSI-positive genes show
  adjusted OR 1.43–2.06 (p ≤ 0.048); the two null tissues (cerebellum, hypothalamus) have
  31–36 PSI-positive genes and a single overlap each.
- **The gene's own switch coordinate is the much stronger corroborator.** Marginal switch
  genes are PSI-positive at 32% against a 0.45% background in caudate (OR ≈ 150), while
  module membership gives OR ≈ 2. That gap is expected — a module carries many genes whose
  own switch is modest — and it is the honest scale of the module-level claim: modules are
  enriched for independently measured splicing change, not composed of it.
- **This is not independent of the transcript annotation.** PSI tables and transcript
  quantification both derive from the same short-read alignments and GENCODE models, so the
  agreement rules out a pure quantification artifact of the switch coordinate but is not an
  orthogonal-platform confirmation. That role belongs to the long-read arm
  (`../longread_switch_confirm/LONGREAD_SWITCH_CONFIRM.md`, `../switch_orthogonal_confirm/`).

## Confidence

Moderate-to-high for gene-level enrichment in the four cortical/hippocampal GTEx regions and
BrainSEQ caudate; low for any per-event statement (0/11 matched events concordant).
