# Candidate RBP regulons among IsoGraph co-switch modules (mature scope)

Motifs are scanned on the **mature transcript** (exonic + UTR) sequence.

For each module and RBP, whether binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool. Two arms are reported for every module x RBP cell:

1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, BH-corrected across all cells. It conditions on nothing else.
2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as covariates, so a module cannot score simply because its genes are long, GC-rich or UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the estimable cells only.

The adjusted arm is the stricter reading and the two can disagree in both directions; where they do, the disagreement is reported below rather than resolved in favour of the larger number.

- switch genes tested: **11637** over 16 regions
- module x RBP cells: **22560**
- hypergeometric q<0.05: **486** (GO-invisible 71)
- estimable GLM cells (BH family): **21095** of 22560
- covariate-adjusted q<0.05: **543**
- of the 486 hypergeometric hits: 485 estimable, **131** also adjusted-q<0.05, 258 with an adjusted CI entirely above 1, **2** entirely *below* it (enriched on raw counts, depleted once opportunity is adjusted for)

### Estimability

A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including such cells would put an arbitrary point estimate into the BH family and dilute every real test. They are excluded from the correction and reported here instead of silently carrying a p-value.

| GLM cell status | n | meaning |
|---|---|---|
| `fit` | 21095 | estimable; carries a p-value and enters the BH family |
| `separated_zero_cell` | 1348 | an empty cell in the (in-module x switched) 2x2 -- complete separation, no finite MLE |
| `quasi_separated` | 62 | fitted probabilities pinned at 0/1, or a degenerate coefficient/standard error |
| `outcome_or_predictor_constant` | 55 | no variation to model in the region |

### Top candidate regulons (by hypergeometric q)

| region | module | RBP | module size | switched | enrichment | q | adj. OR (95% CI) | adj. q | GO-inv |
|---|---|---|---|---|---|---|---|---|---|
| caudate | M004 | DHX9 | 444.0 | 248.0 | 1.37 | 4.47e-08 | 1.53 (1.23-1.91) | 1.55e-02 | no |
| caudate | M004 | ADAR | 444.0 | 271.0 | 1.33 | 4.47e-08 | 1.60 (1.28-2.00) | 6.04e-03 | no |
| hippocampus | M000 | ELAVL4 | 1481.0 | 348.0 | 1.31 | 1.82e-07 | 1.20 (0.98-1.46) | 4.21e-01 | no |
| cerebellar_hemisphere | M000 | RBM41 | 1667.0 | 819.0 | 1.13 | 1.52e-06 | 1.11 (0.96-1.29) | 5.49e-01 | no |
| hippocampus | M000 | ELAVL1 | 1481.0 | 251.0 | 1.36 | 2.78e-06 | 1.08 (0.86-1.35) | 8.39e-01 | no |
| cerebellar_hemisphere | M000 | RBM6 | 1667.0 | 738.0 | 1.14 | 5.01e-06 | 1.09 (0.95-1.27) | 6.34e-01 | no |
| frontal_cortex_ba9 | M001 | IGF2BP1 | 1182.0 | 469.0 | 1.21 | 7.55e-06 | 1.29 (1.12-1.49) | 2.68e-02 | yes |
| nucleus_accumbens_basal_ganglia | M010 | TIAL1 | 69.0 | 29.0 | 3.12 | 7.58e-06 | 0.78 (0.44-1.38) | 7.66e-01 | no |
| caudate | M004 | DDX58 | 444.0 | 219.0 | 1.34 | 7.58e-06 | 1.38 (1.11-1.72) | 1.02e-01 | no |
| caudate | M004 | EIF4A3 | 444.0 | 357.0 | 1.16 | 1.06e-05 | 1.90 (1.47-2.46) | 7.41e-04 | no |
| frontal_cortex_ba9 | M001 | G3BP2 | 1182.0 | 633.0 | 1.15 | 1.14e-05 | 1.23 (1.07-1.41) | 8.60e-02 | yes |
| caudate | M004 | ZNF346 | 444.0 | 172.0 | 1.42 | 1.33e-05 | 1.61 (1.28-2.02) | 6.94e-03 | no |
| frontal_cortex_ba9 | M001 | RBM6 | 1182.0 | 519.0 | 1.19 | 1.90e-05 | 1.20 (1.04-1.38) | 1.64e-01 | yes |
| cerebellar_hemisphere | M000 | FXR1 | 1667.0 | 602.0 | 1.15 | 2.40e-05 | 1.20 (1.03-1.39) | 2.33e-01 | no |
| caudate | M004 | HNRNPCL1 | 444.0 | 374.0 | 1.14 | 3.59e-05 | 1.95 (1.48-2.57) | 1.16e-03 | no |
| putamen_basal_ganglia | M001 | ACO1 | 947.0 | 591.0 | 1.11 | 3.90e-05 | 1.76 (1.44-2.13) | 4.82e-05 | no |
| nucleus_accumbens_basal_ganglia | M010 | ELAVL4 | 69.0 | 28.0 | 2.89 | 5.09e-05 | 0.72 (0.41-1.27) | 6.62e-01 | no |
| nucleus_accumbens_basal_ganglia | M001 | SAMD4A | 1253.0 | 756.0 | 1.11 | 6.15e-05 | 1.38 (1.20-1.59) | 1.92e-03 | no |
| putamen_basal_ganglia | M001 | G3BP2 | 947.0 | 496.0 | 1.13 | 6.15e-05 | 1.63 (1.34-1.97) | 7.41e-04 | no |
| nucleus_accumbens_basal_ganglia | M001 | YTHDC1 | 1253.0 | 860.0 | 1.09 | 6.51e-05 | 1.42 (1.23-1.65) | 1.10e-03 | no |
| cerebellar_hemisphere | M000 | SNRNP70 | 1667.0 | 596.0 | 1.14 | 6.51e-05 | 1.14 (0.98-1.32) | 4.54e-01 | no |
| nucleus_accumbens_basal_ganglia | M010 | NUDT21 | 69.0 | 36.0 | 2.30 | 6.51e-05 | 1.07 (0.63-1.80) | 9.49e-01 | no |
| nucleus_accumbens_basal_ganglia | M010 | OAS1 | 69.0 | 21.0 | 3.64 | 7.33e-05 | 1.04 (0.58-1.89) | 9.69e-01 | no |
| amygdala | M001 | G3BP1 | 997.0 | 564.0 | 1.14 | 9.63e-05 | 0.99 (0.83-1.17) | 9.70e-01 | no |
| caudate | M004 | FMR1 | 444.0 | 137.0 | 1.46 | 1.07e-04 | 1.55 (1.22-1.97) | 2.50e-02 | no |
| cerebellar_hemisphere | M000 | RBM14 | 1667.0 | 834.0 | 1.10 | 1.07e-04 | 1.15 (1.00-1.33) | 3.42e-01 | no |
| amygdala | M001 | G3BP2 | 997.0 | 508.0 | 1.15 | 1.12e-04 | 1.09 (0.92-1.28) | 7.14e-01 | no |
| caudate | M004 | RBMS3 | 444.0 | 282.0 | 1.22 | 1.16e-04 | 1.45 (1.16-1.81) | 4.49e-02 | no |
| hippocampus | M000 | ZFP36 | 1481.0 | 404.0 | 1.21 | 1.17e-04 | 1.12 (0.94-1.34) | 5.98e-01 | no |
| caudate | M004 | NOVA2 | 444.0 | 261.0 | 1.24 | 1.21e-04 | 1.92 (1.54-2.39) | 3.01e-05 | no |

### Contradicted by opportunity adjustment

These 2 cells are hypergeometric-significant yet have an adjusted 95% CI entirely below 1 -- the raw enrichment is explained, and then some, by the motif opportunity their genes carry. They must not be described as candidate regulons.

| region | module | RBP | enrichment | q | adj. OR (95% CI) |
|---|---|---|---|---|---|
| nucleus_accumbens_basal_ganglia | M010 | KHDRBS1 | 1.72 | 6.53e-03 | 0.48 (0.28-0.84) |
| nucleus_accumbens_basal_ganglia | M010 | HNRNPC | 1.80 | 1.95e-02 | 0.49 (0.29-0.85) |

