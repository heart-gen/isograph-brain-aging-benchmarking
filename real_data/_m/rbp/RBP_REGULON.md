# Candidate RBP regulons among IsoGraph co-switch modules (mature scope)

Motifs are scanned on the **mature transcript** (exonic + UTR) sequence.

For each module and RBP, whether binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool. Two arms are reported for every module x RBP cell:

1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, BH-corrected across all cells. It conditions on nothing else.
2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as covariates, so a module cannot score simply because its genes are long, GC-rich or UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the estimable cells only.

The adjusted arm is the stricter reading and the two can disagree in both directions; where they do, the disagreement is reported below rather than resolved in favour of the larger number.

- switch genes tested: **11859** over 14 regions
- module x RBP cells: **34240**
- hypergeometric q<0.05: **714** (GO-invisible 149)
- estimable GLM cells (BH family): **31452** of 34240
- covariate-adjusted q<0.05: **310**
- of the 714 hypergeometric hits: 713 estimable, **43** also adjusted-q<0.05, 344 with an adjusted CI entirely above 1, **19** entirely *below* it (enriched on raw counts, depleted once opportunity is adjusted for)

### Estimability

A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including such cells would put an arbitrary point estimate into the BH family and dilute every real test. They are excluded from the correction and reported here instead of silently carrying a p-value.

| GLM cell status | n | meaning |
|---|---|---|
| `fit` | 31452 | estimable; carries a p-value and enters the BH family |
| `separated_zero_cell` | 2529 | an empty cell in the (in-module x switched) 2x2 -- complete separation, no finite MLE |
| `outcome_or_predictor_constant` | 179 | no variation to model in the region |
| `quasi_separated` | 80 | fitted probabilities pinned at 0/1, or a degenerate coefficient/standard error |

### Top candidate regulons (by hypergeometric q)

| region | module | RBP | module size | switched | enrichment | q | adj. OR (95% CI) | adj. q | GO-inv |
|---|---|---|---|---|---|---|---|---|---|
| frontal_cortex_ba9 | M008 | ELAVL1 | 273.0 | 97.0 | 3.37 | 2.93e-25 | 1.82 (1.36-2.44) | 1.52e-02 | no |
| frontal_cortex_ba9 | M008 | ELAVL4 | 273.0 | 120.0 | 2.69 | 4.73e-24 | 1.44 (1.09-1.91) | 2.25e-01 | no |
| frontal_cortex_ba9 | M008 | KHDRBS1 | 273.0 | 179.0 | 1.82 | 4.98e-20 | 1.37 (1.03-1.83) | 3.48e-01 | no |
| frontal_cortex_ba9 | M008 | TIAL1 | 273.0 | 109.0 | 2.59 | 6.75e-20 | 1.37 (1.03-1.82) | 3.48e-01 | no |
| frontal_cortex_ba9 | M008 | RNASEL | 273.0 | 115.0 | 2.29 | 6.99e-17 | 1.05 (0.79-1.39) | 9.42e-01 | no |
| frontal_cortex_ba9 | M008 | PPIE | 273.0 | 200.0 | 1.51 | 1.07e-13 | 1.45 (1.07-1.96) | 2.73e-01 | no |
| frontal_cortex_ba9 | M008 | ZFP36 | 273.0 | 114.0 | 2.10 | 1.11e-13 | 1.37 (1.04-1.81) | 3.29e-01 | no |
| frontal_cortex_ba9 | M007 | ADAR | 406.0 | 258.0 | 1.45 | 3.90e-13 | 1.25 (1.00-1.57) | 4.28e-01 | no |
| frontal_cortex_ba9 | M008 | OAS1 | 273.0 | 72.0 | 2.78 | 4.20e-13 | 1.57 (1.15-2.14) | 1.60e-01 | no |
| frontal_cortex_ba9 | M001 | G3BP1 | 677.0 | 429.0 | 1.30 | 7.03e-13 | 1.23 (1.04-1.47) | 2.88e-01 | no |
| frontal_cortex_ba9 | M007 | DDX58 | 406.0 | 223.0 | 1.54 | 8.95e-13 | 1.24 (1.00-1.55) | 4.36e-01 | no |
| frontal_cortex_ba9 | M008 | HNRNPC | 273.0 | 136.0 | 1.83 | 1.13e-12 | 1.18 (0.90-1.55) | 6.89e-01 | no |
| frontal_cortex_ba9 | M008 | TIA1 | 273.0 | 78.0 | 2.56 | 1.32e-12 | 1.24 (0.92-1.68) | 6.37e-01 | no |
| frontal_cortex_ba9 | M007 | DHX9 | 406.0 | 237.0 | 1.48 | 2.28e-12 | 1.13 (0.90-1.41) | 7.49e-01 | no |
| frontal_cortex_ba9 | M008 | IGF2BP3 | 273.0 | 97.0 | 2.20 | 2.52e-12 | 1.44 (1.09-1.91) | 2.34e-01 | no |
| frontal_cortex_ba9 | M008 | HNRNPD | 273.0 | 171.0 | 1.60 | 2.70e-12 | 1.25 (0.95-1.66) | 5.56e-01 | no |
| frontal_cortex_ba9 | M007 | MBNL1 | 406.0 | 110.0 | 2.05 | 2.06e-11 | 1.34 (1.04-1.73) | 3.11e-01 | no |
| frontal_cortex_ba9 | M008 | ELAVL2 | 273.0 | 55.0 | 3.09 | 2.07e-11 | 1.17 (0.83-1.66) | 7.98e-01 | no |
| frontal_cortex_ba9 | M001 | RBM14 | 677.0 | 390.0 | 1.32 | 2.48e-11 | 1.31 (1.10-1.55) | 1.01e-01 | no |
| frontal_cortex_ba9 | M007 | ZNF346 | 406.0 | 161.0 | 1.70 | 4.47e-11 | 1.41 (1.12-1.76) | 1.34e-01 | no |
| frontal_cortex_ba9 | M007 | PPRC1 | 406.0 | 213.0 | 1.51 | 4.88e-11 | 1.44 (1.16-1.78) | 6.74e-02 | no |
| frontal_cortex_ba9 | M005 | KHDRBS1 | 471.0 | 247.0 | 1.45 | 5.30e-11 | 0.92 (0.74-1.14) | 8.30e-01 | no |
| frontal_cortex_ba9 | M001 | RBM41 | 677.0 | 360.0 | 1.33 | 3.33e-10 | 1.08 (0.91-1.29) | 7.91e-01 | no |
| frontal_cortex_ba9 | M008 | PABPC1 | 273.0 | 174.0 | 1.49 | 1.01e-09 | 1.22 (0.92-1.61) | 6.38e-01 | no |
| frontal_cortex_ba9 | M008 | RBMX | 273.0 | 38.0 | 3.51 | 6.89e-09 | 1.68 (1.12-2.51) | 2.35e-01 | no |
| anterior_cingulate_cortex_ba24 | M023 | TIAL1 | 77.0 | 34.0 | 3.41 | 9.60e-09 | 2.36 (1.43-3.89) | 5.97e-02 | no |
| frontal_cortex_ba9 | M005 | U2AF2 | 471.0 | 302.0 | 1.30 | 1.03e-08 | 0.98 (0.79-1.21) | 9.70e-01 | no |
| frontal_cortex_ba9 | M001 | SNRPB2 | 677.0 | 367.0 | 1.29 | 1.10e-08 | 1.13 (0.96-1.34) | 6.10e-01 | no |
| frontal_cortex_ba9 | M008 | RC3H1 | 273.0 | 71.0 | 2.23 | 1.82e-08 | 1.30 (0.96-1.77) | 5.28e-01 | no |
| frontal_cortex_ba9 | M005 | PPIE | 471.0 | 298.0 | 1.30 | 2.03e-08 | 0.83 (0.67-1.02) | 5.02e-01 | no |

### Contradicted by opportunity adjustment

These 19 cells are hypergeometric-significant yet have an adjusted 95% CI entirely below 1 -- the raw enrichment is explained, and then some, by the motif opportunity their genes carry. They must not be described as candidate regulons.

| region | module | RBP | enrichment | q | adj. OR (95% CI) |
|---|---|---|---|---|---|
| frontal_cortex_ba9 | M005 | PABPC1 | 1.28 | 1.81e-05 | 0.77 (0.63-0.95) |
| frontal_cortex_ba9 | M005 | HNRNPD | 1.30 | 1.98e-05 | 0.72 (0.59-0.89) |
| frontal_cortex_ba9 | M005 | HNRNPC | 1.40 | 2.40e-05 | 0.78 (0.63-0.97) |
| frontal_cortex_ba9 | M005 | DAZAP1 | 1.38 | 1.87e-04 | 0.79 (0.64-0.99) |
| frontal_cortex_ba9 | M005 | CPEB4 | 1.19 | 5.43e-04 | 0.80 (0.65-0.98) |
| frontal_cortex_ba9 | M005 | ZFP36 | 1.42 | 5.81e-04 | 0.73 (0.58-0.92) |
| frontal_cortex_ba9 | M005 | SF1 | 1.22 | 2.82e-03 | 0.79 (0.65-0.97) |
| frontal_cortex_ba9 | M005 | NUDT21 | 1.31 | 4.23e-03 | 0.68 (0.55-0.85) |
| anterior_cingulate_cortex_ba24 | M011 | HNRNPDL | 1.33 | 4.26e-03 | 0.65 (0.45-0.94) |
| anterior_cingulate_cortex_ba24 | M011 | KHDRBS1 | 1.42 | 8.43e-03 | 0.56 (0.38-0.82) |
| anterior_cingulate_cortex_ba24 | M011 | SF1 | 1.36 | 8.51e-03 | 0.69 (0.48-0.98) |
| anterior_cingulate_cortex_ba24 | M011 | SRSF10 | 1.34 | 1.17e-02 | 0.68 (0.48-0.98) |
| anterior_cingulate_cortex_ba24 | M011 | PABPC1 | 1.34 | 1.30e-02 | 0.66 (0.46-0.95) |
| frontal_cortex_ba9 | M005 | NOVA2 | 1.19 | 1.41e-02 | 0.72 (0.59-0.88) |
| anterior_cingulate_cortex_ba24 | M011 | IGF2BP3 | 1.71 | 1.70e-02 | 0.61 (0.40-0.94) |
| anterior_cingulate_cortex_ba24 | M011 | HNRNPK | 1.35 | 1.93e-02 | 0.60 (0.42-0.86) |
| frontal_cortex_ba9 | M005 | HNRNPA2B1 | 1.34 | 2.00e-02 | 0.70 (0.55-0.89) |
| anterior_cingulate_cortex_ba24 | M011 | PCBP1 | 1.52 | 2.34e-02 | 0.64 (0.44-0.94) |
| anterior_cingulate_cortex_ba24 | M011 | CPEB4 | 1.24 | 3.14e-02 | 0.62 (0.43-0.90) |

