# Validation plan — what would move the isoform-switch layer up a journal tier

**Status:** draft for co-author discussion, 2026-09-09.
**Scope:** the two experiments that would convert a well-controlled set-level correlation into a
mechanism. Everything quoted here is recomputed from committed repository outputs; provenance
for each number is in the appendix.

---

## 0. Where the evidence actually stands

### 0.1 What carries the paper

| Result | Effect | Control that makes it credible |
| --- | --- | --- |
| Splicing-specificity contrast, phenotype-associated modules | **1.111** (95% CI 1.048–1.177), p = 3.6e-4, I² = 0.23 | Matched WGCNA baselines on *identical* features: null (switch-only 1.015 p = 0.74; multiplex 1.044 p = 0.070) |
| Same contrast, IsoGraph-only method effect (13 shared tissues) | **1.112**, p = 6.5e-4 | Isolates VAE + Leiden inference from the feature matrix |
| Long-read confirmation, genetically anchored pairs | switch-like **0.453 vs matched null 0.252**, p = 5e-4 (n = 53 detected pairs) | Abundance-decile-matched pairs from the *same genes*, *same samples*, *same code* |
| …restricted to usably expressed anchored isoforms | **0.600 vs 0.304**, p = 5e-4 (n = 30) | as above |
| Module trust funnel | **250/266** chance-trusted (94%), 6 regions | Split-half permutation null |
| Cross-cohort aging replication | **23/130**, matching-permutation p = 0.034 | Independent cohorts, different quantifiers |
| Switch consequence | UTR-remodeled **1.279** (9/9 regions), CDS-changed **1.044** (9/9); NMD, biotype and coding-status all *depleted* | Within-gene permutation null |

The consequence row is the one that dictates assay design, and it is usually skimmed: **these are
productive UTR/CDS remodeling events, not decay**. An assay that only measures protein-coding
change is aimed at the smaller of the two enriched classes.

### 0.2 The three counter-lines

These are not reviewer speculation — they are the project's own results, and they converge on the
same target from three directions.

| # | Counter-line | Number | What it threatens |
| --- | --- | --- | --- |
| **1** | Per-gene sQTL-vs-eQTL coloc contrast is null | Locus-matched switch 20/46 vs non-switch 88/188, **Fisher p = 0.88**; splicing share ~31% in all four gene pools | The set-level effect has no per-gene counterpart, yet the proposed mechanism is per-gene |
| **2** | The GO-invisible arm is null | **1.068**, p = 0.077 (was 1.172, p = 2.3e-5 before the 2026-08-29 input refresh) | The distinctive framing collapses to the blander "phenotype-associated modules" |
| **3** | S-LDSC points at expression, not splicing | sQTL coefficient clears nominal p < 0.05 in **1 of 6** trait-contexts (AD, p = 0.0355) with **no multiple-testing correction applied**; eQTL clears in 4 of 6. Joint model, SCZ: sQTL p = 0.931 vs eQTL p = **0.0064** | For a splicing paper, the partitioned-heritability evidence favours the abundance annotation |

Counter-line 3 is the least discussed internally and the most dangerous externally, because
S-LDSC is the size-robust test the field trusts and it is currently summarised in-repo as
supporting the switch layer without stating that the splicing arm does not survive correction.

**Neither experiment below rescues counter-line 3.** It is a power and annotation-resolution
problem in GTEx bulk brain sQTL, and the honest handling is disclosure plus the framing fix in
§3.3 — not a new experiment.

---

## 1. Experiment A — functional validation of a switch

**Question.** Does an anchored isoform switch *do* something to the transcript, or is it a
correlated measurement?

### 1.1 The targeting trap

The 12 genetically anchored, splicing-led genes split into two disjoint groups, and the split runs
exactly opposite to what a naive target pick would assume.

| Gene | Traits | max CLPP | max anchored IF | Long-read confirmed | Group |
| --- | --- | --- | --- | --- | --- |
| **PGS1** | ALS | **0.093** | 0.818 | yes | **dual evidence** |
| **PPP6R2** | ALS, SCZ | **0.065** | 0.904 | yes | **dual evidence** |
| CDIP1 | SCZ | 0.028 | 0.306 | yes | dual evidence |
| DLG1 | SCZ | 0.026 | 0.400 | yes | dual evidence |
| TBC1D15 | PD | 0.017 | 0.894 | yes | dual evidence |
| PRRC2B | SCZ | 0.017 | 0.281 | yes | dual evidence |
| TPCN1 | AD | 0.011 | 0.894 | yes | dual evidence |
| GGNBP2 | ALS | 0.011 | 0.930 | yes | dual evidence |
| ARVCF | SCZ | 0.010 | 0.405 | yes | dual evidence |
| **CTSH** | AD | **0.386** | **0.004** | **no** | genetics only |
| **SNCA** | LBD, PD | 0.038 | **0.003** | yes (not at usable abundance) | genetics only |
| RTEL1 | SCZ | 0.010 | 0.545 | no | neither |

**The gene with by far the strongest colocalization — CTSH, CLPP 0.386, an order of magnitude
above every other gene — has an anchored isoform at 0.4% usage and fails orthogonal
confirmation.** SNCA is the same shape: the narratively attractive locus sits at 0.3% long-read
usage. Meanwhile the genes whose switches are unambiguously real (IF 0.28–0.93, all long-read
confirmed) carry CLPP of 0.010–0.093.

Picking targets on genetic evidence alone lands on the two genes that cannot be assayed. Picking
on expression alone lands on genes a reviewer will call genetically weak. **Validate the
intersection and say why.**

### 1.2 Recommended targets

- **Primary: PGS1** (CLPP 0.093, IF 0.818, ALS) and **PPP6R2** (CLPP 0.065, IF 0.904, ALS + SCZ).
  Best joint genetic and expression evidence; PPP6R2 is additionally cross-trait.
- **Secondary: DLG1** (CLPP 0.026, IF 0.400, SCZ) — a synaptic scaffold, so a functional readout
  has interpretable neuronal meaning.
- **Narrative only: SNCA.** Already confirmed at the junction level in short read (minor-form
  usage **0.189** DLPFC n = 222, **0.234** caudate n = 238, against a pre-registered 0.05
  threshold). Do **not** build the functional assay on it — at 0.3% long-read usage the
  effect size available to any biochemical readout is marginal.

### 1.3 Assay, chosen from the consequence data

Because UTR remodeling is the dominant enriched consequence (1.279, 9/9 regions) and NMD is
*depleted*, the readout should be UTR-dependent function, not decay and not protein truncation:

1. **Isoform-resolved quantification** in an independent brain panel — targeted amplicon or
   Nanopore cDNA over the two competing isoforms. Establishes the switch exists outside BrainSEQ
   and GTEx. *This is the minimum viable experiment and could stand alone.*
2. **UTR reporter.** Clone the alternative 3′/5′ UTR of each target into a dual-luciferase or
   destabilised-GFP backbone; measure steady-state output and, with actinomycin-D chase,
   transcript half-life. Read out **relative** difference between the two UTRs, in a neuronal
   line (SH-SY5Y or iPSC-derived neurons).
3. **Allelic direction (PGS1, PPP6R2 only).** In heterozygous donors, test whether the risk
   allele shifts isoform usage in the direction the coloc predicts. This is the step that ties
   the switch to the *genetics* rather than to aging alone.

### 1.4 Pre-registered decision rule

Write this before running, in the manner already used for SNCA/CTSH:

> A target is **functionally validated** if (a) both isoforms are detected in the independent
> panel at ≥ 5% usage, **and** (b) the two UTRs differ in reporter output or half-life at
> p < 0.05 with the direction stated in advance. A target is **not** validated if either fails;
> report it and keep it off the main figure.

Two of three targets validating supports a main-figure mechanism panel. Zero of three means the
paper stays a set-level correlation and should be framed that way — which is publishable, just
not at the tier this experiment is meant to reach.

---

## 2. Experiment B — snRNA-seq and the composition question

**Question.** Do the switches persist *within* matched cell types, or are they a proportion shift?

### 2.1 What is at stake

| Contrast | comp_unique base → adjusted | Retained |
| --- | --- | --- |
| BrainSEQ SCZD caudate | 34 → **2** | **6%** |
| BrainSEQ aging caudate | 43 → 17 | 40% |
| BrainSEQ aging DLPFC | 8 → 15 | 188% (adjustment *increased* it) |
| BrainSEQ aging hippocampus | 0 → 0 | undefined — contributes nothing, must not be counted as passing |
| GTEx frontal cortex BA9 | 531 → **0** | **0%** |
| GTEx cortex | 438 → **0** | **0%** |
| GTEx anterior cingulate BA24 | 545 → 88 | 16% |
| GTEx hippocampus | 61 → 21 | 34% |
| GTEx caudate | 32 → 10 | 31% |

The disease arm is already conceded as composition-entangled. The sharper problem is that **both
GTEx cortical regions collapse to exactly zero** while limbic and striatal regions retain 16–34%.
Bulk deconvolution cannot say whether that is real regional biology or over-adjustment, because
proportion and per-cell composition are confounded in the same measurement by construction.

### 2.2 The technical constraint that decides the design

**Standard 3′-biased snRNA-seq cannot measure isoform usage.** A 10x 3′ library reads the last
few hundred bases of a transcript; most of the switches here are UTR or internal-CDS remodeling
events that a 3′ tag cannot distinguish. Running conventional snRNA-seq and hoping to recover
isoform ratios would produce a null that means nothing.

Viable options, in order of preference:

1. **10x 3′ + long-read on the same libraries** (PacBio Kinnex / MAS-ISO-seq, or ONT). Cell-type
   labels from the short reads, isoform quantification from the long reads on the identical
   nuclei. This is the design that answers the question.
2. **Targeted capture** of the specific junctions for the anchored gene set. Far cheaper, but
   only tests the genes already nominated — it cannot re-derive the layer, so it settles the
   named loci and not the composition question in general.
3. **Single-nucleus long-read alone.** Cleanest, most expensive, lowest throughput per donor.

This project has an existing 3′ bias sensitivity finding; the same physics applies here and
should be cited in the methods rather than rediscovered.

### 2.3 Power — and why the obvious design fails

Detectable ΔPSI between two groups (two-sided α = 0.05, 80% power), by donors per group and
between-donor PSI standard deviation:

| n / group | SD = 0.05 | SD = 0.08 | SD = 0.10 | SD = 0.15 |
| --- | --- | --- | --- | --- |
| 10 | 0.063 | 0.100 | 0.125 | 0.188 |
| 20 | 0.044 | 0.071 | 0.089 | 0.133 |
| 30 | 0.036 | 0.058 | 0.072 | 0.109 |
| 40 | 0.031 | 0.050 | 0.063 | 0.094 |
| 50 | 0.028 | 0.045 | 0.056 | 0.084 |

Age as a **continuous** predictor is the trap:

| Donors | Detectable \|r\| | R² |
| --- | --- | --- |
| 20 | 0.591 | 0.349 |
| 40 | 0.431 | 0.185 |
| 60 | 0.355 | 0.126 |
| 100 | 0.277 | 0.077 |
| **238 (the bulk BrainSEQ caudate n)** | **0.181** | 0.033 |

**A 40-donor snRNA-seq study detects \|r\| ≥ 0.43; the bulk study that generated the hypothesis
detects \|r\| ≥ 0.18.** Any realistically sized single-nucleus aging study is *worse powered than
the experiment it is meant to adjudicate*, so a null result would be uninterpretable. Do not run
this as an age-correlation study.

**Design that works:** frame it as a *composition* test, not an aging test. Sample donors at the
extremes of the cell-type proportion distribution (or case/control for the SCZD arm), and ask
whether the switch persists within matched cell types across that contrast. Extreme sampling buys
back effect size that continuous age cannot. n = 20–30 per group detects ΔPSI ≈ 0.04–0.09 at
plausible SDs, which brackets the bulk effects.

### 2.4 Pre-registered decision rule

> The disease/aging switch layer is **cell-intrinsic** if the anchored switches show ΔPSI in the
> stated direction within at least one matched cell type at FDR < 0.05, with the cell-type
> proportion shift explicitly modelled. It is **compositional** if the switches are present
> across the proportion contrast in pseudobulk but vanish within every matched cell type at
> adequate power. Report the power achieved per cell type; a cell type with fewer than 20 donors
> contributing nuclei is reported as untested, not as negative.

---

## 3. Recommendation

### 3.1 If you run one thing

**Experiment A, step 1 alone** — isoform-resolved quantification of PGS1, PPP6R2 and DLG1 in an
independent brain panel. It is the cheapest item on this page, it is the one a reviewer will
demand first, and it either produces a main-figure mechanism panel or tells you to stop early.

### 3.2 If you run two

Add **Experiment B option 2** (targeted junction capture on snRNA-seq nuclei) for the same gene
set. Together these answer "is the switch real and does it live in a cell type" for the named
loci, without committing to a full long-read single-nucleus atlas.

### 3.3 Framing fixes that cost nothing and should happen regardless

1. `MANUSCRIPT_PLAN.md` §12 still lists **"Splicing-QTL specifically anchor the GO-invisible
   switch layer"** as *Established* and as the Fig 3 headline. It is null (1.068, p = 0.077).
   `FIGURE_ORDERING.md` was corrected to the phenotype-associated claim; the claims table was
   not. **These two documents currently disagree about the paper's headline.**
2. State the S-LDSC splicing arm as it is: one of six nominal, uncorrected. Presenting it as
   support invites counter-line 3 to be discovered by a reviewer instead of disclosed by you.
3. `reports/pi/00_OVERVIEW.md` §8 is stale — 89% should be 94%, and SNCA is no longer "uncertain".

---

## Appendix — provenance

| Number | Source |
| --- | --- |
| 1.111 / 1.112 contrasts, baseline nulls | `05_genetic_anchoring/_m/qtl_anchoring_meta/QTL_ANCHORING_META.md` (re-run 2026-09-09 after the matched-baseline re-fit) |
| 0.453 / 0.252, 0.600 / 0.304 | `06_switch_mechanism/_m/switch_orthogonal_confirm/anchored_summary.json` |
| 12-gene anchored table, CLPP, IF | `06_switch_mechanism/_m/switch_orthogonal_confirm/anchored_gene_confirmation.parquet` |
| UTR 1.279 / CDS 1.044 / NMD depleted | `06_switch_mechanism/_m/switch_consequence_meta.parquet`, `analysis_class = aging, stratum = all` |
| 250/266, 23/130 p = 0.034 | `03_module_trust/_m/stability/module_trust/`, re-run 2026-09-09 on corrected split halves |
| SNCA 0.189 / 0.234 | `06_switch_mechanism/_m/junction_coloc_confirm/` |
| Composition table | `04_module_characterization/_m/COMPOSITION_ADJUSTMENT_SUMMARY.md` |
| Per-gene coloc Fisher p = 0.88 | `05_genetic_anchoring/_m/coloc_modality_contrast/` |
| S-LDSC tau p-values | `05_genetic_anchoring/_m/ldsc/LDSC_SUMMARY.md` |
| Power tables | Computed 2026-09-09; two-sample t and Fisher-z, α = 0.05 two-sided, power = 0.80 |
