# Candidate RBP regulons among IsoGraph co-switch modules (mature scope)

Motifs are scanned on the **mature transcript** (exonic + UTR) sequence.

For each module and motif family, whether binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool. Two arms are reported for every module x motif family cell:

1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, BH-corrected across all cells. It conditions on nothing else.
2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as covariates, so a module cannot score simply because its genes are long, GC-rich or UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the estimable cells only.

The adjusted arm is the stricter reading and the two can disagree in both directions; where they do, the disagreement is reported below rather than resolved in favour of the larger number.

- switch genes tested: **11859** over 14 regions
- module x motif family cells: **29104**
- hypergeometric q<0.05: **555** (GO-invisible 108)
- estimable GLM cells (BH family): **28285** of 29104
- covariate-adjusted q<0.05: **261**
- of the 555 hypergeometric hits: 555 estimable, **36** also adjusted-q<0.05, 243 with an adjusted CI entirely above 1, **22** entirely *below* it (enriched on raw counts, depleted once opportunity is adjusted for)

### Estimability

A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including such cells would put an arbitrary point estimate into the BH family and dilute every real test. They are excluded from the correction and reported here instead of silently carrying a p-value.

| GLM cell status | n | meaning |
|---|---|---|
| `fit` | 28285 | estimable; carries a p-value and enters the BH family |
| `separated_zero_cell` | 756 | an empty cell in the (in-module x switched) 2x2 -- complete separation, no finite MLE |
| `outcome_or_predictor_constant` | 44 | no variation to model in the region |
| `quasi_separated` | 19 | fitted probabilities pinned at 0/1, or a degenerate coefficient/standard error |

### Top candidate regulons (by hypergeometric q)

| region | module | motif family | module size | switched | enrichment | q | adj. OR (95% CI) | adj. q | GO-inv |
|---|---|---|---|---|---|---|---|---|---|
| frontal_cortex_ba9 | M008 | F0112 | 273.0 | 105.0 | 3.12 | 3.81e-25 | 1.55 (1.16-2.07) | 1.29e-01 | no |
| frontal_cortex_ba9 | M008 | F0056 | 273.0 | 161.0 | 1.93 | 3.24e-19 | 1.24 (0.94-1.64) | 5.98e-01 | no |
| frontal_cortex_ba9 | M007 | F0116 | 406.0 | 168.0 | 1.99 | 1.58e-18 | 1.29 (1.03-1.63) | 3.52e-01 | no |
| frontal_cortex_ba9 | M007 | F0065 | 406.0 | 289.0 | 1.42 | 3.10e-15 | 1.32 (1.04-1.66) | 3.08e-01 | no |
| frontal_cortex_ba9 | M007 | F0003 | 406.0 | 207.0 | 1.67 | 3.29e-15 | 1.30 (1.04-1.63) | 3.09e-01 | no |
| frontal_cortex_ba9 | M008 | F0114 | 273.0 | 54.0 | 3.62 | 9.21e-14 | 1.66 (1.16-2.36) | 1.67e-01 | no |
| frontal_cortex_ba9 | M007 | F0085 | 406.0 | 83.0 | 2.60 | 2.23e-13 | 1.58 (1.19-2.10) | 9.85e-02 | no |
| frontal_cortex_ba9 | M008 | F0054 | 273.0 | 96.0 | 2.29 | 4.02e-13 | 0.99 (0.74-1.32) | 9.87e-01 | no |
| frontal_cortex_ba9 | M001 | F0065 | 677.0 | 439.0 | 1.29 | 7.51e-13 | 1.30 (1.09-1.56) | 1.41e-01 | no |
| frontal_cortex_ba9 | M008 | F0113 | 273.0 | 101.0 | 2.13 | 4.21e-12 | 1.40 (1.06-1.85) | 2.90e-01 | no |
| frontal_cortex_ba9 | M007 | F0032 | 406.0 | 214.0 | 1.54 | 4.21e-12 | 1.07 (0.85-1.34) | 8.89e-01 | no |
| frontal_cortex_ba9 | M007 | F0001 | 406.0 | 253.0 | 1.43 | 5.88e-12 | 1.19 (0.95-1.49) | 5.99e-01 | no |
| frontal_cortex_ba9 | M008 | F0108 | 273.0 | 100.0 | 2.13 | 6.15e-12 | 1.41 (1.07-1.87) | 2.73e-01 | no |
| frontal_cortex_ba9 | M008 | F0064 | 273.0 | 127.0 | 1.85 | 7.68e-12 | 1.59 (1.21-2.08) | 6.72e-02 | no |
| frontal_cortex_ba9 | M008 | F0053 | 273.0 | 190.0 | 1.48 | 1.77e-11 | 1.39 (1.04-1.86) | 3.44e-01 | no |
| frontal_cortex_ba9 | M001 | F0014 | 677.0 | 341.0 | 1.35 | 2.66e-10 | 1.20 (1.01-1.42) | 3.88e-01 | no |
| frontal_cortex_ba9 | M007 | F0111 | 406.0 | 282.0 | 1.33 | 3.47e-10 | 1.42 (1.13-1.79) | 1.18e-01 | no |
| frontal_cortex_ba9 | M008 | F0127 | 273.0 | 79.0 | 2.21 | 1.82e-09 | 1.45 (1.08-1.95) | 2.64e-01 | no |
| frontal_cortex_ba9 | M001 | F0130 | 677.0 | 311.0 | 1.37 | 2.09e-09 | 1.16 (0.98-1.38) | 5.40e-01 | no |
| frontal_cortex_ba9 | M008 | F0069 | 273.0 | 162.0 | 1.53 | 2.34e-09 | 1.22 (0.93-1.60) | 6.32e-01 | no |
| frontal_cortex_ba9 | M001 | F0062 | 677.0 | 280.0 | 1.41 | 2.34e-09 | 1.21 (1.02-1.44) | 3.57e-01 | no |
| frontal_cortex_ba9 | M007 | F0004 | 406.0 | 116.0 | 1.85 | 2.61e-09 | 1.52 (1.19-1.94) | 6.90e-02 | no |
| frontal_cortex_ba9 | M007 | F0091 | 406.0 | 101.0 | 1.96 | 3.66e-09 | 1.62 (1.26-2.10) | 3.41e-02 | no |
| frontal_cortex_ba9 | M001 | F0039 | 677.0 | 458.0 | 1.22 | 4.01e-09 | 1.13 (0.94-1.35) | 6.68e-01 | no |
| frontal_cortex_ba9 | M008 | F0051 | 273.0 | 58.0 | 2.58 | 6.65e-09 | 1.15 (0.83-1.61) | 8.24e-01 | no |
| frontal_cortex_ba9 | M007 | F0039 | 406.0 | 289.0 | 1.29 | 1.03e-08 | 1.49 (1.18-1.88) | 6.93e-02 | no |
| frontal_cortex_ba9 | M007 | F0086 | 406.0 | 286.0 | 1.29 | 1.14e-08 | 1.61 (1.27-2.04) | 1.87e-02 | no |
| frontal_cortex_ba9 | M007 | F0087 | 406.0 | 66.0 | 2.37 | 1.14e-08 | 1.36 (1.00-1.84) | 4.40e-01 | no |
| frontal_cortex_ba9 | M007 | F0006 | 406.0 | 57.0 | 2.57 | 1.34e-08 | 1.60 (1.15-2.22) | 1.64e-01 | no |
| frontal_cortex_ba9 | M007 | F0021 | 406.0 | 184.0 | 1.50 | 2.29e-08 | 1.08 (0.87-1.35) | 8.64e-01 | no |

### Contradicted by opportunity adjustment

These 22 cells are hypergeometric-significant yet have an adjusted 95% CI entirely below 1 -- the raw enrichment is explained, and then some, by the motif opportunity their genes carry. They must not be described as candidate regulons.

| region | module | motif family | enrichment | q | adj. OR (95% CI) |
|---|---|---|---|---|---|
| frontal_cortex_ba9 | M005 | F0052 | 1.34 | 7.32e-06 | 0.81 (0.66-0.99) |
| frontal_cortex_ba9 | M005 | F0053 | 1.26 | 7.95e-06 | 0.73 (0.59-0.90) |
| frontal_cortex_ba9 | M005 | F0056 | 1.37 | 2.52e-05 | 0.78 (0.63-0.97) |
| frontal_cortex_ba9 | M005 | F0070 | 1.34 | 3.43e-05 | 0.77 (0.62-0.95) |
| frontal_cortex_ba9 | M005 | F0045 | 1.30 | 1.55e-04 | 0.81 (0.66-1.00) |
| frontal_cortex_ba9 | M005 | F0103 | 1.26 | 4.19e-04 | 0.79 (0.64-0.96) |
| frontal_cortex_ba9 | M005 | F0079 | 1.23 | 8.58e-04 | 0.75 (0.61-0.92) |
| frontal_cortex_ba9 | M005 | F0050 | 1.21 | 2.11e-03 | 0.74 (0.60-0.91) |
| frontal_cortex_ba9 | M005 | F0067 | 1.17 | 3.21e-03 | 0.75 (0.61-0.92) |
| frontal_cortex_ba9 | M005 | F0058 | 1.13 | 1.43e-02 | 0.74 (0.59-0.91) |
| frontal_cortex_ba9 | M005 | F0134 | 1.43 | 1.75e-02 | 0.74 (0.57-0.97) |
| frontal_cortex_ba9 | M005 | F0126 | 1.13 | 1.76e-02 | 0.79 (0.64-0.98) |
| anterior_cingulate_cortex_ba24 | M011 | F0048 | 1.34 | 1.77e-02 | 0.67 (0.47-0.96) |
| anterior_cingulate_cortex_ba24 | M011 | F0079 | 1.32 | 2.07e-02 | 0.59 (0.41-0.84) |
| anterior_cingulate_cortex_ba24 | M011 | F0068 | 1.32 | 2.07e-02 | 0.64 (0.45-0.91) |
| anterior_cingulate_cortex_ba24 | M011 | F0080 | 1.58 | 2.09e-02 | 0.60 (0.41-0.90) |
| anterior_cingulate_cortex_ba24 | M011 | F0069 | 1.34 | 3.01e-02 | 0.58 (0.40-0.84) |
| frontal_cortex_ba9 | M005 | F0064 | 1.26 | 3.21e-02 | 0.65 (0.52-0.81) |
| anterior_cingulate_cortex_ba24 | M011 | F0108 | 1.64 | 3.33e-02 | 0.58 (0.38-0.88) |
| anterior_cingulate_cortex_ba24 | M011 | F0103 | 1.33 | 3.61e-02 | 0.61 (0.43-0.87) |

