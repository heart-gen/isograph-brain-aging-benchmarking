# Candidate RBP regulons among IsoGraph co-switch modules (mature scope)

Motifs are scanned on the **mature transcript** (exonic + UTR) sequence.

For each module and motif family, whether binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool. Two arms are reported for every module x motif family cell:

1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, BH-corrected across all cells. It conditions on nothing else.
2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as covariates, so a module cannot score simply because its genes are long, GC-rich or UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the estimable cells only.

The adjusted arm is the stricter reading and the two can disagree in both directions; where they do, the disagreement is reported below rather than resolved in favour of the larger number.

- switch genes tested: **11790** over 16 regions
- module x motif family cells: **9384**
- hypergeometric q<0.05: **336** (GO-invisible 41)
- estimable GLM cells (BH family): **8833** of 9384
- covariate-adjusted q<0.05: **467**
- of the 336 hypergeometric hits: 335 estimable, **99** also adjusted-q<0.05, 183 with an adjusted CI entirely above 1, **1** entirely *below* it (enriched on raw counts, depleted once opportunity is adjusted for)

### Estimability

A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including such cells would put an arbitrary point estimate into the BH family and dilute every real test. They are excluded from the correction and reported here instead of silently carrying a p-value.

| GLM cell status | n | meaning |
|---|---|---|
| `fit` | 8833 | estimable; carries a p-value and enters the BH family |
| `outcome_or_predictor_constant` | 414 | no variation to model in the region |
| `separated_zero_cell` | 131 | an empty cell in the (in-module x switched) 2x2 -- complete separation, no finite MLE |
| `quasi_separated` | 6 | fitted probabilities pinned at 0/1, or a degenerate coefficient/standard error |

### Top candidate regulons (by hypergeometric q)

| region | module | motif family | module size | switched | enrichment | q | adj. OR (95% CI) | adj. q | GO-inv |
|---|---|---|---|---|---|---|---|---|---|
| cerebellar_hemisphere | M000 | F0039 | 2700.0 | 1729.0 | 1.10 | 1.46e-16 | 1.18 (1.04-1.34) | 1.16e-01 | no |
| cerebellar_hemisphere | M000 | F0031 | 2700.0 | 1408.0 | 1.13 | 9.14e-16 | 1.21 (1.07-1.37) | 5.00e-02 | no |
| cerebellar_hemisphere | M000 | F0111 | 2700.0 | 1590.0 | 1.10 | 5.82e-13 | 1.13 (0.99-1.28) | 3.23e-01 | no |
| cerebellar_hemisphere | M000 | F0065 | 2700.0 | 1526.0 | 1.10 | 1.84e-12 | 1.10 (0.97-1.25) | 4.96e-01 | no |
| cerebellar_hemisphere | M000 | F0130 | 2700.0 | 1099.0 | 1.14 | 1.84e-12 | 1.13 (0.99-1.29) | 3.29e-01 | no |
| cerebellar_hemisphere | M000 | F0096 | 2700.0 | 1095.0 | 1.14 | 2.43e-12 | 1.12 (0.98-1.27) | 3.97e-01 | no |
| cerebellar_hemisphere | M000 | F0073 | 2700.0 | 1529.0 | 1.09 | 1.53e-10 | 1.12 (0.99-1.27) | 3.56e-01 | no |
| cerebellar_hemisphere | M000 | F0014 | 2700.0 | 1296.0 | 1.11 | 1.53e-10 | 1.13 (0.99-1.27) | 3.29e-01 | no |
| cerebellar_hemisphere | M000 | F0011 | 2700.0 | 1149.0 | 1.12 | 2.09e-10 | 1.03 (0.91-1.17) | 8.94e-01 | no |
| cerebellar_hemisphere | M000 | F0042 | 2700.0 | 1519.0 | 1.09 | 4.27e-10 | 1.17 (1.03-1.32) | 1.28e-01 | no |
| cerebellar_hemisphere | M000 | F0047 | 2700.0 | 1724.0 | 1.07 | 5.49e-09 | 1.19 (1.05-1.35) | 7.93e-02 | no |
| cerebellar_hemisphere | M000 | F0060 | 2700.0 | 1843.0 | 1.06 | 1.59e-08 | 1.23 (1.08-1.40) | 3.96e-02 | no |
| cerebellar_hemisphere | M000 | F0040 | 2700.0 | 606.0 | 1.18 | 1.59e-08 | 1.37 (1.17-1.61) | 5.32e-03 | no |
| cerebellar_hemisphere | M000 | F0029 | 2700.0 | 1680.0 | 1.07 | 1.66e-08 | 1.17 (1.03-1.33) | 1.29e-01 | no |
| cerebellar_hemisphere | M000 | F0062 | 2700.0 | 861.0 | 1.13 | 1.09e-07 | 1.11 (0.97-1.28) | 4.56e-01 | no |
| cerebellar_hemisphere | M000 | F0035 | 2700.0 | 1619.0 | 1.07 | 1.40e-07 | 1.26 (1.12-1.43) | 9.02e-03 | no |
| cerebellar_hemisphere | M000 | F0104 | 2700.0 | 1142.0 | 1.10 | 1.40e-07 | 1.01 (0.88-1.14) | 9.84e-01 | no |
| cerebellar_hemisphere | M000 | F0135 | 2700.0 | 1655.0 | 1.07 | 1.40e-07 | 1.19 (1.05-1.35) | 8.00e-02 | no |
| cerebellar_hemisphere | M000 | F0136 | 2700.0 | 1793.0 | 1.06 | 2.03e-07 | 1.15 (1.01-1.31) | 2.23e-01 | no |
| cerebellar_hemisphere | M000 | F0121 | 2700.0 | 817.0 | 1.13 | 2.03e-07 | 1.09 (0.95-1.25) | 6.01e-01 | no |
| cerebellar_hemisphere | M000 | F0016 | 2700.0 | 1094.0 | 1.10 | 2.33e-07 | 1.09 (0.96-1.24) | 5.43e-01 | no |
| cerebellar_hemisphere | M003 | F0044 | 933.0 | 162.0 | 1.50 | 1.17e-06 | 1.30 (1.06-1.61) | 1.35e-01 | no |
| cerebellar_hemisphere | M000 | F0030 | 2700.0 | 1962.0 | 1.05 | 1.17e-06 | 1.21 (1.06-1.39) | 6.60e-02 | no |
| cerebellar_hemisphere | M000 | F0033 | 2700.0 | 1903.0 | 1.05 | 2.87e-06 | 1.17 (1.03-1.34) | 1.57e-01 | no |
| putamen_basal_ganglia | M001 | F0133 | 1357.0 | 865.0 | 1.10 | 2.87e-06 | 1.60 (1.37-1.85) | 1.71e-06 | no |
| cerebellar_hemisphere | M000 | F0023 | 2700.0 | 1276.0 | 1.08 | 3.74e-06 | 1.05 (0.92-1.19) | 8.01e-01 | no |
| cerebellar_hemisphere | M000 | F0124 | 2700.0 | 1727.0 | 1.06 | 4.13e-06 | 1.12 (0.99-1.27) | 3.42e-01 | no |
| cerebellar_hemisphere | M000 | F0002 | 2700.0 | 1707.0 | 1.06 | 4.31e-06 | 1.22 (1.07-1.39) | 4.99e-02 | no |
| cerebellum | M000 | F0031 | 1590.0 | 849.0 | 1.10 | 4.43e-06 | 1.18 (1.02-1.37) | 2.23e-01 | no |
| cerebellar_hemisphere | M000 | F0074 | 2700.0 | 1660.0 | 1.06 | 6.60e-06 | 1.01 (0.89-1.15) | 9.56e-01 | no |

### Contradicted by opportunity adjustment

These 1 cells are hypergeometric-significant yet have an adjusted 95% CI entirely below 1 -- the raw enrichment is explained, and then some, by the motif opportunity their genes carry. They must not be described as candidate regulons.

| region | module | motif family | enrichment | q | adj. OR (95% CI) |
|---|---|---|---|---|---|
| cerebellar_hemisphere | M003 | F0069 | 1.11 | 4.40e-02 | 0.80 (0.68-0.94) |

