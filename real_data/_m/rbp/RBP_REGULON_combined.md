# Candidate RBP regulons among IsoGraph co-switch modules (combined scope)

Motif presence unions the **mature-transcript** and **intronic splice-site flank** scans.

For each module and RBP, whether binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool. Two arms are reported for every module x RBP cell:

1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, BH-corrected across all cells. It conditions on nothing else.
2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as covariates, so a module cannot score simply because its genes are long, GC-rich or UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the estimable cells only.

The adjusted arm is the stricter reading and the two can disagree in both directions; where they do, the disagreement is reported below rather than resolved in favour of the larger number.

- switch genes tested: **11859** over 14 regions
- module x RBP cells: **34240**
- hypergeometric q<0.05: **725** (GO-invisible 163)
- estimable GLM cells (BH family): **29207** of 34240
- covariate-adjusted q<0.05: **223**
- of the 725 hypergeometric hits: 724 estimable, **38** also adjusted-q<0.05, 304 with an adjusted CI entirely above 1, **12** entirely *below* it (enriched on raw counts, depleted once opportunity is adjusted for)

### Estimability

A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including such cells would put an arbitrary point estimate into the BH family and dilute every real test. They are excluded from the correction and reported here instead of silently carrying a p-value.

| GLM cell status | n | meaning |
|---|---|---|
| `fit` | 29207 | estimable; carries a p-value and enters the BH family |
| `separated_zero_cell` | 3789 | an empty cell in the (in-module x switched) 2x2 -- complete separation, no finite MLE |
| `outcome_or_predictor_constant` | 1146 | no variation to model in the region |
| `quasi_separated` | 98 | fitted probabilities pinned at 0/1, or a degenerate coefficient/standard error |

### Top candidate regulons (by hypergeometric q)

| region | module | RBP | module size | switched | enrichment | q | adj. OR (95% CI) | adj. q | GO-inv |
|---|---|---|---|---|---|---|---|---|---|
| frontal_cortex_ba9 | M008 | PPIE | 273.0 | 169.0 | 2.35 | 3.55e-32 | 1.56 (1.17-2.06) | 1.19e-01 | no |
| frontal_cortex_ba9 | M008 | KHDRBS1 | 273.0 | 130.0 | 2.78 | 1.16e-28 | 1.25 (0.94-1.66) | 5.84e-01 | no |
| frontal_cortex_ba9 | M008 | U2AF2 | 273.0 | 142.0 | 2.33 | 4.78e-24 | 1.27 (0.97-1.67) | 5.16e-01 | no |
| frontal_cortex_ba9 | M007 | PPRC1 | 406.0 | 286.0 | 1.59 | 5.29e-24 | 1.97 (1.56-2.49) | 2.63e-05 | no |
| frontal_cortex_ba9 | M008 | ELAVL3 | 273.0 | 154.0 | 2.15 | 1.77e-23 | 1.24 (0.94-1.63) | 5.92e-01 | no |
| frontal_cortex_ba9 | M008 | CPEB1 | 273.0 | 153.0 | 2.16 | 2.35e-23 | 1.21 (0.92-1.59) | 6.38e-01 | no |
| frontal_cortex_ba9 | M008 | PABPC1 | 273.0 | 156.0 | 1.99 | 6.33e-20 | 1.37 (1.04-1.81) | 3.27e-01 | no |
| frontal_cortex_ba9 | M008 | ELAVL4 | 273.0 | 64.0 | 4.00 | 1.43e-19 | 1.57 (1.12-2.19) | 2.25e-01 | no |
| frontal_cortex_ba9 | M008 | HNRNPC | 273.0 | 82.0 | 2.90 | 7.63e-17 | 1.27 (0.94-1.71) | 5.86e-01 | no |
| frontal_cortex_ba9 | M001 | RBM14 | 677.0 | 470.0 | 1.31 | 7.80e-17 | 1.49 (1.25-1.79) | 6.11e-03 | no |
| frontal_cortex_ba9 | M001 | PPRC1 | 677.0 | 411.0 | 1.37 | 5.21e-16 | 1.29 (1.08-1.54) | 1.60e-01 | no |
| frontal_cortex_ba9 | M008 | CPEB4 | 273.0 | 148.0 | 1.85 | 3.01e-15 | 1.23 (0.94-1.61) | 5.84e-01 | no |
| frontal_cortex_ba9 | M007 | DDX58 | 406.0 | 188.0 | 1.70 | 2.73e-14 | 1.28 (1.03-1.60) | 3.44e-01 | no |
| frontal_cortex_ba9 | M008 | RNASEL | 273.0 | 77.0 | 2.73 | 3.64e-14 | 0.99 (0.72-1.35) | 9.86e-01 | no |
| frontal_cortex_ba9 | M008 | TIAL1 | 273.0 | 58.0 | 3.41 | 3.64e-14 | 1.21 (0.86-1.71) | 7.35e-01 | no |
| frontal_cortex_ba9 | M008 | HNRNPA0 | 273.0 | 146.0 | 1.81 | 5.08e-14 | 1.03 (0.79-1.35) | 9.62e-01 | no |
| frontal_cortex_ba9 | M007 | ADAR | 406.0 | 227.0 | 1.55 | 7.60e-14 | 1.21 (0.97-1.51) | 5.27e-01 | no |
| frontal_cortex_ba9 | M007 | RBM14 | 406.0 | 294.0 | 1.37 | 1.63e-13 | 1.80 (1.42-2.27) | 1.11e-03 | no |
| frontal_cortex_ba9 | M008 | HNRNPD | 273.0 | 107.0 | 2.12 | 3.28e-13 | 1.08 (0.82-1.43) | 8.87e-01 | no |
| frontal_cortex_ba9 | M007 | GRSF1 | 406.0 | 257.0 | 1.44 | 4.48e-13 | 1.28 (1.03-1.60) | 3.51e-01 | no |
| frontal_cortex_ba9 | M008 | ELAVL1 | 273.0 | 37.0 | 4.63 | 1.98e-12 | 1.68 (1.11-2.55) | 2.68e-01 | no |
| frontal_cortex_ba9 | M007 | DHX9 | 406.0 | 205.0 | 1.57 | 2.95e-12 | 1.10 (0.88-1.38) | 8.05e-01 | no |
| frontal_cortex_ba9 | M007 | RBMS1 | 406.0 | 207.0 | 1.53 | 3.27e-11 | 1.60 (1.29-1.99) | 7.87e-03 | no |
| anterior_cingulate_cortex_ba24 | M023 | PABPC1 | 77.0 | 50.0 | 2.55 | 2.27e-10 | 2.51 (1.51-4.18) | 5.01e-02 | no |
| frontal_cortex_ba9 | M001 | RBM41 | 677.0 | 392.0 | 1.30 | 2.27e-10 | 1.13 (0.95-1.34) | 6.45e-01 | no |
| anterior_cingulate_cortex_ba24 | M023 | ELAVL3 | 77.0 | 48.0 | 2.65 | 2.43e-10 | 2.12 (1.27-3.53) | 1.53e-01 | no |
| frontal_cortex_ba9 | M008 | DDX19B | 273.0 | 167.0 | 1.54 | 3.44e-10 | 1.16 (0.88-1.51) | 7.48e-01 | no |
| frontal_cortex_ba9 | M008 | NUDT21 | 273.0 | 70.0 | 2.44 | 3.60e-10 | 1.27 (0.93-1.73) | 5.99e-01 | no |
| anterior_cingulate_cortex_ba24 | M023 | CPEB1 | 77.0 | 47.0 | 2.67 | 3.71e-10 | 1.99 (1.20-3.30) | 2.07e-01 | no |
| frontal_cortex_ba9 | M001 | G3BP1 | 677.0 | 466.0 | 1.23 | 5.19e-10 | 1.30 (1.08-1.55) | 1.67e-01 | no |

### Contradicted by opportunity adjustment

These 12 cells are hypergeometric-significant yet have an adjusted 95% CI entirely below 1 -- the raw enrichment is explained, and then some, by the motif opportunity their genes carry. They must not be described as candidate regulons.

| region | module | RBP | enrichment | q | adj. OR (95% CI) |
|---|---|---|---|---|---|
| frontal_cortex_ba9 | M005 | PABPC1 | 1.40 | 5.42e-06 | 0.78 (0.63-0.97) |
| frontal_cortex_ba9 | M005 | KHDRBS3 | 1.25 | 7.38e-05 | 0.81 (0.66-0.99) |
| frontal_cortex_ba9 | M005 | PUM1 | 1.28 | 5.36e-04 | 0.81 (0.66-1.00) |
| frontal_cortex_ba9 | M005 | HNRNPA0 | 1.30 | 1.45e-03 | 0.80 (0.64-0.99) |
| frontal_cortex_ba9 | M005 | KHDRBS2 | 1.18 | 3.17e-03 | 0.79 (0.64-0.97) |
| frontal_cortex_ba9 | M005 | HNRNPU | 1.23 | 3.41e-03 | 0.77 (0.63-0.94) |
| cerebellar_hemisphere | M013 | KHDRBS1 | 2.22 | 3.56e-03 | 0.51 (0.28-0.95) |
| anterior_cingulate_cortex_ba24 | M011 | HNRNPDL | 1.48 | 1.19e-02 | 0.65 (0.45-0.94) |
| frontal_cortex_ba9 | M005 | CPEB4 | 1.25 | 1.64e-02 | 0.69 (0.56-0.86) |
| frontal_cortex_ba9 | M005 | NOVA2 | 1.32 | 1.81e-02 | 0.76 (0.60-0.95) |
| frontal_cortex_ba9 | M005 | RBFOX1 | 1.18 | 3.37e-02 | 0.80 (0.66-0.98) |
| caudate_basal_ganglia | M014 | TARDBP | 1.40 | 4.54e-02 | 0.59 (0.39-0.90) |

