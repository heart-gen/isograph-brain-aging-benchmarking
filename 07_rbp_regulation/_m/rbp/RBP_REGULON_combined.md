# Candidate RBP regulons among IsoGraph co-switch modules (combined scope)

Motif presence unions the **mature-transcript** and **intronic splice-site flank** scans.

For each module and RBP, whether binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool. Two arms are reported for every module x RBP cell:

1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, BH-corrected across all cells. It conditions on nothing else.
2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as covariates, so a module cannot score simply because its genes are long, GC-rich or UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the estimable cells only.

The adjusted arm is the stricter reading and the two can disagree in both directions; where they do, the disagreement is reported below rather than resolved in favour of the larger number.

- switch genes tested: **11637** over 16 regions
- module x RBP cells: **22560**
- hypergeometric q<0.05: **618** (GO-invisible 98)
- estimable GLM cells (BH family): **19700** of 22560
- covariate-adjusted q<0.05: **496**
- of the 618 hypergeometric hits: 618 estimable, **154** also adjusted-q<0.05, 324 with an adjusted CI entirely above 1, **12** entirely *below* it (enriched on raw counts, depleted once opportunity is adjusted for)

### Estimability

A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including such cells would put an arbitrary point estimate into the BH family and dilute every real test. They are excluded from the correction and reported here instead of silently carrying a p-value.

| GLM cell status | n | meaning |
|---|---|---|
| `fit` | 19700 | estimable; carries a p-value and enters the BH family |
| `separated_zero_cell` | 2043 | an empty cell in the (in-module x switched) 2x2 -- complete separation, no finite MLE |
| `outcome_or_predictor_constant` | 681 | no variation to model in the region |
| `quasi_separated` | 136 | fitted probabilities pinned at 0/1, or a degenerate coefficient/standard error |

### Top candidate regulons (by hypergeometric q)

| region | module | RBP | module size | switched | enrichment | q | adj. OR (95% CI) | adj. q | GO-inv |
|---|---|---|---|---|---|---|---|---|---|
| caudate | M004 | DHX9 | 444.0 | 208.0 | 1.43 | 1.39e-07 | 1.50 (1.20-1.88) | 2.38e-02 | no |
| caudate | M004 | ADAR | 444.0 | 227.0 | 1.39 | 1.39e-07 | 1.51 (1.21-1.89) | 2.00e-02 | no |
| cerebellar_hemisphere | M000 | RBM6 | 1667.0 | 883.0 | 1.13 | 1.39e-07 | 1.22 (1.06-1.40) | 1.23e-01 | no |
| cerebellar_hemisphere | M000 | RBM8A | 1667.0 | 752.0 | 1.15 | 2.23e-07 | 1.24 (1.07-1.43) | 9.67e-02 | no |
| frontal_cortex_ba9 | M001 | SAMD4A | 1182.0 | 789.0 | 1.14 | 2.23e-07 | 1.41 (1.22-1.62) | 1.47e-03 | yes |
| frontal_cortex_ba9 | M001 | G3BP2 | 1182.0 | 699.0 | 1.15 | 9.55e-07 | 1.30 (1.13-1.49) | 1.91e-02 | yes |
| nucleus_accumbens_basal_ganglia | M001 | ENOX1 | 1253.0 | 809.0 | 1.12 | 1.17e-06 | 1.46 (1.26-1.68) | 3.74e-04 | no |
| frontal_cortex_ba9 | M001 | IGF2BP1 | 1182.0 | 567.0 | 1.19 | 1.66e-06 | 1.25 (1.09-1.43) | 5.78e-02 | yes |
| cerebellar_hemisphere | M000 | SNRNP70 | 1667.0 | 736.0 | 1.14 | 3.11e-06 | 1.19 (1.03-1.37) | 2.15e-01 | no |
| substantia_nigra | M000 | G3BP2 | 1357.0 | 771.0 | 1.11 | 3.11e-06 | 1.25 (1.06-1.46) | 1.32e-01 | no |
| nucleus_accumbens_basal_ganglia | M001 | SAMD4A | 1253.0 | 840.0 | 1.11 | 3.13e-06 | 1.50 (1.29-1.73) | 1.79e-04 | no |
| nucleus_accumbens_basal_ganglia | M010 | PPIE | 69.0 | 37.0 | 2.54 | 3.50e-06 | 0.86 (0.50-1.50) | 8.82e-01 | no |
| amygdala | M001 | SAMD4A | 997.0 | 654.0 | 1.13 | 3.83e-06 | 1.33 (1.12-1.57) | 4.69e-02 | no |
| caudate | M004 | DDX58 | 444.0 | 175.0 | 1.43 | 3.83e-06 | 1.37 (1.09-1.73) | 1.25e-01 | no |
| hippocampus | M000 | ELAVL4 | 1481.0 | 141.0 | 1.51 | 3.90e-06 | 1.38 (1.03-1.84) | 2.79e-01 | no |
| amygdala | M001 | RBM6 | 997.0 | 512.0 | 1.17 | 4.24e-06 | 1.15 (0.97-1.35) | 4.82e-01 | no |
| cerebellar_hemisphere | M000 | ZC3H10 | 1667.0 | 867.0 | 1.11 | 4.24e-06 | 1.14 (0.99-1.31) | 3.93e-01 | no |
| caudate | M004 | G3BP1 | 444.0 | 316.0 | 1.21 | 4.27e-06 | 1.63 (1.29-2.05) | 5.63e-03 | no |
| caudate | M004 | ZNF346 | 444.0 | 102.0 | 1.68 | 4.50e-06 | 1.76 (1.34-2.30) | 6.81e-03 | no |
| amygdala | M001 | PPRC1 | 997.0 | 537.0 | 1.16 | 4.53e-06 | 1.09 (0.92-1.29) | 7.17e-01 | no |
| amygdala | M001 | ZC3H10 | 997.0 | 531.0 | 1.16 | 4.99e-06 | 1.22 (1.04-1.44) | 2.13e-01 | no |
| nucleus_accumbens_basal_ganglia | M001 | CNOT4 | 1253.0 | 748.0 | 1.12 | 6.57e-06 | 1.41 (1.23-1.63) | 1.07e-03 | no |
| amygdala | M001 | RBM14 | 997.0 | 625.0 | 1.13 | 6.57e-06 | 1.14 (0.96-1.34) | 5.41e-01 | no |
| frontal_cortex_ba9 | M001 | CNOT4 | 1182.0 | 718.0 | 1.13 | 6.57e-06 | 1.24 (1.08-1.42) | 7.77e-02 | yes |
| hippocampus | M000 | PPIE | 1481.0 | 461.0 | 1.21 | 1.28e-05 | 1.05 (0.88-1.26) | 8.72e-01 | no |
| hippocampus | M000 | CPEB1 | 1481.0 | 495.0 | 1.19 | 1.28e-05 | 1.02 (0.86-1.22) | 9.49e-01 | no |
| substantia_nigra | M000 | PPRC1 | 1357.0 | 745.0 | 1.11 | 1.28e-05 | 1.13 (0.96-1.33) | 5.63e-01 | no |
| cerebellar_hemisphere | M000 | RBM42 | 1667.0 | 552.0 | 1.16 | 1.28e-05 | 1.11 (0.95-1.30) | 5.95e-01 | no |
| caudate | M004 | RBMS3 | 444.0 | 297.0 | 1.22 | 1.28e-05 | 1.53 (1.22-1.91) | 1.66e-02 | no |
| cerebellar_hemisphere | M000 | PPRC1 | 1667.0 | 896.0 | 1.10 | 1.31e-05 | 1.09 (0.94-1.26) | 6.80e-01 | no |

### Contradicted by opportunity adjustment

These 12 cells are hypergeometric-significant yet have an adjusted 95% CI entirely below 1 -- the raw enrichment is explained, and then some, by the motif opportunity their genes carry. They must not be described as candidate regulons.

| region | module | RBP | enrichment | q | adj. OR (95% CI) |
|---|---|---|---|---|---|
| nucleus_accumbens_basal_ganglia | M010 | ELAVL3 | 2.25 | 3.26e-05 | 0.54 (0.30-0.95) |
| nucleus_accumbens_basal_ganglia | M010 | U2AF2 | 2.39 | 5.44e-05 | 0.54 (0.31-0.96) |
| nucleus_accumbens_basal_ganglia | M010 | HNRNPA0 | 1.97 | 4.77e-04 | 0.54 (0.31-0.94) |
| nucleus_accumbens_basal_ganglia | M010 | KHDRBS1 | 2.52 | 6.42e-04 | 0.37 (0.20-0.67) |
| nucleus_accumbens_basal_ganglia | M010 | HNRNPD | 2.39 | 8.98e-04 | 0.54 (0.31-0.96) |
| nucleus_accumbens_basal_ganglia | M010 | HNRNPC | 2.97 | 1.07e-03 | 0.50 (0.27-0.95) |
| amygdala | M009 | KHDRBS1 | 1.90 | 1.33e-02 | 0.54 (0.31-0.95) |
| frontal_cortex_ba9 | M004 | KHDRBS1 | 1.34 | 1.73e-02 | 0.68 (0.51-0.91) |
| nucleus_accumbens_basal_ganglia | M010 | CPEB1 | 1.82 | 2.02e-02 | 0.30 (0.17-0.53) |
| frontal_cortex_ba9 | M004 | CPEB1 | 1.25 | 2.91e-02 | 0.76 (0.59-0.97) |
| cerebellar_hemisphere | M003 | KHDRBS1 | 1.31 | 3.20e-02 | 0.69 (0.52-0.92) |
| frontal_cortex_ba9 | M004 | U2AF2 | 1.27 | 3.56e-02 | 0.75 (0.57-0.97) |

