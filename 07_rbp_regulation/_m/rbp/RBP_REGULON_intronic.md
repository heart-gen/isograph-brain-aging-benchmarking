# Candidate RBP regulons among IsoGraph co-switch modules (intronic scope)

Motifs are scanned on **intronic splice-site flanks** (pre-mRNA sense; up to 100 nt into each intron), the binding niche for splicing-regulatory RBPs invisible to the mature-transcript scan.

For each module and RBP, whether binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool. Two arms are reported for every module x RBP cell:

1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, BH-corrected across all cells. It conditions on nothing else.
2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as covariates, so a module cannot score simply because its genes are long, GC-rich or UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the estimable cells only.

The adjusted arm is the stricter reading and the two can disagree in both directions; where they do, the disagreement is reported below rather than resolved in favour of the larger number.

- switch genes tested: **11637** over 16 regions
- module x RBP cells: **22560**
- hypergeometric q<0.05: **1031** (GO-invisible 235)
- estimable GLM cells (BH family): **22470** of 22560
- covariate-adjusted q<0.05: **412**
- of the 1031 hypergeometric hits: 1031 estimable, **203** also adjusted-q<0.05, 542 with an adjusted CI entirely above 1, **2** entirely *below* it (enriched on raw counts, depleted once opportunity is adjusted for)

### Estimability

A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including such cells would put an arbitrary point estimate into the BH family and dilute every real test. They are excluded from the correction and reported here instead of silently carrying a p-value.

| GLM cell status | n | meaning |
|---|---|---|
| `fit` | 22470 | estimable; carries a p-value and enters the BH family |
| `separated_zero_cell` | 90 | an empty cell in the (in-module x switched) 2x2 -- complete separation, no finite MLE |

### Top candidate regulons (by hypergeometric q)

| region | module | RBP | module size | switched | enrichment | q | adj. OR (95% CI) | adj. q | GO-inv |
|---|---|---|---|---|---|---|---|---|---|
| substantia_nigra | M000 | SAMD4A | 1357.0 | 582.0 | 1.21 | 3.60e-12 | 1.40 (1.18-1.66) | 1.19e-02 | no |
| amygdala | M001 | SAMD4A | 997.0 | 434.0 | 1.30 | 3.60e-12 | 1.34 (1.13-1.59) | 5.02e-02 | no |
| substantia_nigra | M000 | HNRNPA3 | 1357.0 | 589.0 | 1.20 | 3.78e-11 | 1.29 (1.09-1.53) | 9.73e-02 | no |
| frontal_cortex_ba9 | M001 | HNRNPA3 | 1182.0 | 523.0 | 1.27 | 5.61e-11 | 1.37 (1.19-1.58) | 3.60e-03 | yes |
| frontal_cortex_ba9 | M001 | CNOT4 | 1182.0 | 423.0 | 1.31 | 3.49e-10 | 1.40 (1.21-1.63) | 2.65e-03 | yes |
| frontal_cortex_ba9 | M001 | SAMD4A | 1182.0 | 498.0 | 1.26 | 1.89e-09 | 1.31 (1.13-1.51) | 2.34e-02 | yes |
| frontal_cortex_ba9 | M001 | ENOX1 | 1182.0 | 417.0 | 1.29 | 3.47e-09 | 1.42 (1.22-1.64) | 1.90e-03 | yes |
| frontal_cortex_ba9 | M001 | SRSF4 | 1182.0 | 648.0 | 1.19 | 3.76e-09 | 1.37 (1.19-1.58) | 2.74e-03 | yes |
| frontal_cortex_ba9 | M001 | RBMS3 | 1182.0 | 358.0 | 1.33 | 6.51e-09 | 1.30 (1.11-1.53) | 5.90e-02 | yes |
| frontal_cortex_ba9 | M001 | LIN28A | 1182.0 | 670.0 | 1.18 | 1.05e-08 | 1.28 (1.11-1.48) | 3.60e-02 | yes |
| amygdala | M001 | HNRNPA3 | 997.0 | 424.0 | 1.25 | 1.33e-08 | 1.23 (1.04-1.46) | 2.45e-01 | no |
| cerebellum | M000 | G3BP1 | 1285.0 | 594.0 | 1.15 | 2.73e-08 | 1.34 (1.10-1.61) | 9.07e-02 | no |
| amygdala | M001 | HNRNPCL1 | 997.0 | 551.0 | 1.18 | 2.96e-08 | 1.18 (0.99-1.39) | 4.28e-01 | no |
| substantia_nigra | M000 | HNRNPAB | 1357.0 | 745.0 | 1.13 | 4.17e-08 | 1.33 (1.13-1.56) | 3.80e-02 | no |
| cerebellar_hemisphere | M000 | RBMS3 | 1667.0 | 493.0 | 1.21 | 5.20e-08 | 1.22 (1.03-1.45) | 2.49e-01 | no |
| substantia_nigra | M000 | HNRNPM | 1357.0 | 730.0 | 1.13 | 5.20e-08 | 1.21 (1.03-1.42) | 2.76e-01 | no |
| substantia_nigra | M000 | RBM14 | 1357.0 | 583.0 | 1.17 | 5.49e-08 | 1.26 (1.06-1.49) | 1.55e-01 | no |
| amygdala | M001 | RBM14 | 997.0 | 432.0 | 1.23 | 5.49e-08 | 1.19 (1.00-1.41) | 3.85e-01 | no |
| nucleus_accumbens_basal_ganglia | M001 | HNRNPAB | 1253.0 | 709.0 | 1.15 | 5.72e-08 | 1.48 (1.29-1.71) | 1.00e-04 | no |
| substantia_nigra | M000 | RBM41 | 1357.0 | 353.0 | 1.24 | 5.87e-08 | 1.28 (1.04-1.56) | 2.42e-01 | no |
| cerebellar_hemisphere | M000 | A1CF | 1667.0 | 569.0 | 1.18 | 5.87e-08 | 1.29 (1.10-1.51) | 6.43e-02 | no |
| cerebellar_hemisphere | M000 | G3BP1 | 1667.0 | 755.0 | 1.14 | 5.87e-08 | 1.20 (1.03-1.39) | 2.30e-01 | no |
| cerebellum | M000 | CELF4 | 1285.0 | 554.0 | 1.15 | 5.87e-08 | 1.36 (1.13-1.64) | 6.29e-02 | no |
| cerebellar_hemisphere | M000 | YTHDC1 | 1667.0 | 714.0 | 1.15 | 5.87e-08 | 1.26 (1.09-1.46) | 7.61e-02 | no |
| cerebellar_hemisphere | M000 | CNOT4 | 1667.0 | 553.0 | 1.19 | 5.93e-08 | 1.28 (1.09-1.50) | 8.09e-02 | no |
| nucleus_accumbens_basal_ganglia | M001 | ZRANB2 | 1253.0 | 804.0 | 1.13 | 6.17e-08 | 1.56 (1.35-1.80) | 8.10e-06 | no |
| frontal_cortex_ba9 | M001 | RBFOX2 | 1182.0 | 697.0 | 1.16 | 6.26e-08 | 1.36 (1.18-1.56) | 3.97e-03 | yes |
| frontal_cortex_ba9 | M001 | YTHDC1 | 1182.0 | 512.0 | 1.22 | 8.27e-08 | 1.29 (1.12-1.49) | 3.02e-02 | yes |
| anterior_cingulate_cortex_ba24 | M003 | RC3H1 | 277.0 | 97.0 | 1.80 | 8.27e-08 | 1.48 (1.10-1.98) | 1.70e-01 | yes |
| nucleus_accumbens_basal_ganglia | M001 | ENOX1 | 1253.0 | 434.0 | 1.24 | 8.27e-08 | 1.47 (1.26-1.71) | 7.14e-04 | no |

### Contradicted by opportunity adjustment

These 2 cells are hypergeometric-significant yet have an adjusted 95% CI entirely below 1 -- the raw enrichment is explained, and then some, by the motif opportunity their genes carry. They must not be described as candidate regulons.

| region | module | RBP | enrichment | q | adj. OR (95% CI) |
|---|---|---|---|---|---|
| nucleus_accumbens_basal_ganglia | M010 | ELAVL1 | 1.96 | 1.75e-02 | 0.55 (0.31-0.97) |
| nucleus_accumbens_basal_ganglia | M010 | TIAL1 | 1.68 | 2.01e-02 | 0.59 (0.34-1.00) |

