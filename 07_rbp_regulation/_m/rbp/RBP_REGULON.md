# Candidate RBP regulons among IsoGraph co-switch modules (mature scope)

Motifs are scanned on the **mature transcript** (exonic + UTR) sequence.

For each module and RBP, whether binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool. Two arms are reported for every module x RBP cell:

1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, BH-corrected across all cells. It conditions on nothing else.
2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as covariates, so a module cannot score simply because its genes are long, GC-rich or UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the estimable cells only.

The adjusted arm is the stricter reading and the two can disagree in both directions; where they do, the disagreement is reported below rather than resolved in favour of the larger number.

- switch genes tested: **11790** over 16 regions
- module x RBP cells: **11040**
- hypergeometric q<0.05: **384** (GO-invisible 61)
- estimable GLM cells (BH family): **9955** of 11040
- covariate-adjusted q<0.05: **416**
- of the 384 hypergeometric hits: 383 estimable, **89** also adjusted-q<0.05, 216 with an adjusted CI entirely above 1, **1** entirely *below* it (enriched on raw counts, depleted once opportunity is adjusted for)

### Estimability

A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including such cells would put an arbitrary point estimate into the BH family and dilute every real test. They are excluded from the correction and reported here instead of silently carrying a p-value.

| GLM cell status | n | meaning |
|---|---|---|
| `fit` | 9955 | estimable; carries a p-value and enters the BH family |
| `outcome_or_predictor_constant` | 530 | no variation to model in the region |
| `separated_zero_cell` | 510 | an empty cell in the (in-module x switched) 2x2 -- complete separation, no finite MLE |
| `quasi_separated` | 45 | fitted probabilities pinned at 0/1, or a degenerate coefficient/standard error |

### Top candidate regulons (by hypergeometric q)

| region | module | RBP | module size | switched | enrichment | q | adj. OR (95% CI) | adj. q | GO-inv |
|---|---|---|---|---|---|---|---|---|---|
| cerebellar_hemisphere | M000 | RBM14 | 2700.0 | 1394.0 | 1.13 | 7.76e-16 | 1.30 (1.15-1.47) | 3.61e-03 | no |
| cerebellar_hemisphere | M000 | RBM41 | 2700.0 | 1217.0 | 1.14 | 8.41e-15 | 1.10 (0.97-1.25) | 5.13e-01 | no |
| cerebellar_hemisphere | M000 | G3BP1 | 2700.0 | 1460.0 | 1.11 | 3.97e-14 | 1.06 (0.93-1.21) | 7.21e-01 | no |
| cerebellar_hemisphere | M000 | RBM6 | 2700.0 | 1231.0 | 1.13 | 1.65e-12 | 1.18 (1.04-1.34) | 1.19e-01 | no |
| cerebellar_hemisphere | M000 | SNRPB2 | 2700.0 | 1322.0 | 1.12 | 2.20e-12 | 1.11 (0.98-1.27) | 4.26e-01 | no |
| cerebellar_hemisphere | M000 | RBMS3 | 2700.0 | 1311.0 | 1.11 | 2.25e-11 | 1.04 (0.91-1.18) | 8.40e-01 | no |
| cerebellar_hemisphere | M000 | PPRC1 | 2700.0 | 1080.0 | 1.13 | 4.04e-11 | 1.13 (0.99-1.29) | 3.73e-01 | no |
| cerebellar_hemisphere | M000 | G3BP2 | 2700.0 | 1259.0 | 1.10 | 4.22e-09 | 1.10 (0.97-1.24) | 5.08e-01 | no |
| cerebellar_hemisphere | M000 | ENOX1 | 2700.0 | 1372.0 | 1.10 | 4.45e-09 | 1.11 (0.98-1.26) | 4.38e-01 | no |
| cerebellar_hemisphere | M000 | ZC3H10 | 2700.0 | 1124.0 | 1.11 | 4.79e-09 | 1.14 (1.00-1.30) | 2.84e-01 | no |
| cerebellar_hemisphere | M000 | IGF2BP1 | 2700.0 | 974.0 | 1.13 | 4.79e-09 | 1.15 (1.01-1.31) | 2.68e-01 | no |
| cerebellar_hemisphere | M000 | CNOT4 | 2700.0 | 1328.0 | 1.10 | 4.79e-09 | 1.02 (0.90-1.15) | 9.42e-01 | no |
| cerebellar_hemisphere | M000 | HNRNPA3 | 2700.0 | 1658.0 | 1.08 | 7.30e-09 | 1.09 (0.96-1.24) | 5.55e-01 | no |
| cerebellar_hemisphere | M000 | CELF4 | 2700.0 | 1593.0 | 1.08 | 1.60e-08 | 1.22 (1.08-1.39) | 3.90e-02 | no |
| cerebellar_hemisphere | M000 | CELF5 | 2700.0 | 1577.0 | 1.08 | 3.22e-08 | 1.19 (1.05-1.35) | 9.00e-02 | no |
| cerebellar_hemisphere | M000 | FXR1 | 2700.0 | 944.0 | 1.12 | 6.36e-08 | 1.19 (1.04-1.36) | 1.35e-01 | no |
| cerebellar_hemisphere | M000 | CELF6 | 2700.0 | 1577.0 | 1.08 | 1.16e-07 | 1.18 (1.05-1.34) | 1.03e-01 | no |
| cerebellar_hemisphere | M003 | HNRNPDL | 933.0 | 561.0 | 1.18 | 1.37e-07 | 1.14 (0.97-1.33) | 4.55e-01 | no |
| putamen_basal_ganglia | M001 | G3BP2 | 1357.0 | 655.0 | 1.15 | 1.38e-07 | 1.51 (1.30-1.75) | 7.90e-05 | no |
| cerebellar_hemisphere | M000 | MATR3 | 2700.0 | 1169.0 | 1.10 | 1.62e-07 | 1.07 (0.94-1.21) | 6.92e-01 | no |
| cerebellar_hemisphere | M000 | SNRNP70 | 2700.0 | 957.0 | 1.12 | 1.75e-07 | 1.16 (1.02-1.33) | 2.25e-01 | no |
| cerebellar_hemisphere | M000 | PABPC3 | 2700.0 | 1260.0 | 1.09 | 2.04e-07 | 1.13 (1.00-1.28) | 3.09e-01 | no |
| cerebellar_hemisphere | M000 | RBM46 | 2700.0 | 1112.0 | 1.10 | 2.67e-07 | 1.07 (0.94-1.22) | 6.62e-01 | no |
| cerebellar_hemisphere | M000 | PABPC5 | 2700.0 | 1682.0 | 1.07 | 4.14e-07 | 1.16 (1.03-1.31) | 1.76e-01 | no |
| cerebellar_hemisphere | M000 | YTHDC1 | 2700.0 | 1795.0 | 1.06 | 8.55e-07 | 1.21 (1.06-1.37) | 6.95e-02 | no |
| cerebellar_hemisphere | M000 | EIF4B | 2700.0 | 1058.0 | 1.10 | 1.29e-06 | 1.14 (1.00-1.30) | 2.89e-01 | no |
| cerebellar_hemisphere | M000 | RBMS1 | 2700.0 | 882.0 | 1.12 | 1.39e-06 | 0.91 (0.79-1.05) | 5.56e-01 | no |
| cerebellar_hemisphere | M000 | A1CF | 2700.0 | 1492.0 | 1.07 | 2.26e-06 | 1.02 (0.91-1.16) | 9.06e-01 | no |
| cerebellar_hemisphere | M003 | HNRNPA0 | 933.0 | 547.0 | 1.17 | 2.26e-06 | 1.09 (0.93-1.27) | 6.67e-01 | no |
| putamen_basal_ganglia | M001 | ACO1 | 1357.0 | 865.0 | 1.10 | 2.71e-06 | 1.60 (1.37-1.85) | 3.84e-06 | no |

### Contradicted by opportunity adjustment

These 1 cells are hypergeometric-significant yet have an adjusted 95% CI entirely below 1 -- the raw enrichment is explained, and then some, by the motif opportunity their genes carry. They must not be described as candidate regulons.

| region | module | RBP | enrichment | q | adj. OR (95% CI) |
|---|---|---|---|---|---|
| cerebellar_hemisphere | M003 | TIAL1 | 1.25 | 6.03e-03 | 0.81 (0.66-1.00) |

