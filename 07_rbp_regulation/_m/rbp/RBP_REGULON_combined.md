# Candidate RBP regulons among IsoGraph co-switch modules (combined scope)

Motif presence unions the **mature-transcript** and **intronic splice-site flank** scans.

For each module and RBP, whether binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool. Two arms are reported for every module x RBP cell:

1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, BH-corrected across all cells. It conditions on nothing else.
2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as covariates, so a module cannot score simply because its genes are long, GC-rich or UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the estimable cells only.

The adjusted arm is the stricter reading and the two can disagree in both directions; where they do, the disagreement is reported below rather than resolved in favour of the larger number.

- switch genes tested: **11790** over 16 regions
- module x RBP cells: **11040**
- hypergeometric q<0.05: **408** (GO-invisible 60)
- estimable GLM cells (BH family): **9387** of 11040
- covariate-adjusted q<0.05: **381**
- of the 408 hypergeometric hits: 407 estimable, **105** also adjusted-q<0.05, 237 with an adjusted CI entirely above 1, **3** entirely *below* it (enriched on raw counts, depleted once opportunity is adjusted for)

### Estimability

A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including such cells would put an arbitrary point estimate into the BH family and dilute every real test. They are excluded from the correction and reported here instead of silently carrying a p-value.

| GLM cell status | n | meaning |
|---|---|---|
| `fit` | 9387 | estimable; carries a p-value and enters the BH family |
| `outcome_or_predictor_constant` | 791 | no variation to model in the region |
| `separated_zero_cell` | 774 | an empty cell in the (in-module x switched) 2x2 -- complete separation, no finite MLE |
| `quasi_separated` | 88 | fitted probabilities pinned at 0/1, or a degenerate coefficient/standard error |

### Top candidate regulons (by hypergeometric q)

| region | module | RBP | module size | switched | enrichment | q | adj. OR (95% CI) | adj. q | GO-inv |
|---|---|---|---|---|---|---|---|---|---|
| cerebellar_hemisphere | M000 | RBM14 | 2700.0 | 1706.0 | 1.12 | 4.48e-22 | 1.36 (1.20-1.54) | 6.84e-04 | no |
| cerebellar_hemisphere | M000 | PPRC1 | 2700.0 | 1410.0 | 1.14 | 4.85e-19 | 1.16 (1.02-1.32) | 1.95e-01 | no |
| cerebellar_hemisphere | M000 | G3BP1 | 2700.0 | 1680.0 | 1.10 | 2.82e-16 | 1.25 (1.10-1.41) | 2.41e-02 | no |
| cerebellar_hemisphere | M000 | ENOX1 | 2700.0 | 1592.0 | 1.11 | 2.61e-15 | 1.28 (1.13-1.45) | 7.58e-03 | no |
| cerebellar_hemisphere | M000 | SNRPB2 | 2700.0 | 1469.0 | 1.12 | 3.86e-15 | 1.18 (1.04-1.34) | 1.13e-01 | no |
| cerebellar_hemisphere | M000 | RBM8A | 2700.0 | 1142.0 | 1.15 | 7.41e-15 | 1.21 (1.07-1.38) | 6.17e-02 | no |
| cerebellar_hemisphere | M000 | ZC3H10 | 2700.0 | 1330.0 | 1.13 | 1.45e-14 | 1.21 (1.07-1.37) | 5.45e-02 | no |
| cerebellar_hemisphere | M000 | MSI1 | 2700.0 | 885.0 | 1.17 | 1.90e-13 | 1.13 (0.98-1.30) | 4.24e-01 | no |
| cerebellar_hemisphere | M000 | RBM28 | 2700.0 | 1022.0 | 1.15 | 7.90e-13 | 1.12 (0.98-1.28) | 4.42e-01 | no |
| cerebellar_hemisphere | M000 | RBMS3 | 2700.0 | 1476.0 | 1.10 | 1.06e-12 | 1.11 (0.98-1.26) | 4.47e-01 | no |
| cerebellar_hemisphere | M000 | RBM46 | 2700.0 | 1362.0 | 1.11 | 1.25e-12 | 1.19 (1.05-1.34) | 1.03e-01 | no |
| cerebellar_hemisphere | M000 | RBMS1 | 2700.0 | 1065.0 | 1.14 | 2.62e-12 | 1.00 (0.87-1.14) | 9.91e-01 | no |
| cerebellar_hemisphere | M000 | RBM41 | 2700.0 | 1349.0 | 1.11 | 3.21e-12 | 1.06 (0.93-1.20) | 7.47e-01 | no |
| cerebellar_hemisphere | M000 | RBM6 | 2700.0 | 1442.0 | 1.10 | 3.21e-12 | 1.18 (1.04-1.34) | 1.16e-01 | no |
| cerebellar_hemisphere | M000 | G3BP2 | 2700.0 | 1393.0 | 1.11 | 5.63e-12 | 1.19 (1.05-1.34) | 1.01e-01 | no |
| cerebellar_hemisphere | M000 | IGF2BP1 | 2700.0 | 1196.0 | 1.12 | 1.42e-11 | 1.11 (0.98-1.27) | 4.18e-01 | no |
| cerebellar_hemisphere | M000 | CELF6 | 2700.0 | 1725.0 | 1.08 | 1.42e-11 | 1.25 (1.11-1.42) | 1.82e-02 | no |
| cerebellar_hemisphere | M000 | RBM42 | 2700.0 | 783.0 | 1.16 | 1.43e-10 | 1.24 (1.08-1.43) | 6.36e-02 | no |
| cerebellar_hemisphere | M000 | SAMD4A | 2700.0 | 1721.0 | 1.08 | 1.07e-09 | 1.18 (1.04-1.34) | 1.16e-01 | no |
| cerebellar_hemisphere | M000 | SNRNP70 | 2700.0 | 1099.0 | 1.12 | 1.08e-09 | 1.20 (1.05-1.36) | 9.19e-02 | no |
| cerebellar_hemisphere | M000 | CNOT4 | 2700.0 | 1558.0 | 1.08 | 3.55e-09 | 1.09 (0.97-1.24) | 5.23e-01 | no |
| cerebellar_hemisphere | M000 | EIF4B | 2700.0 | 1235.0 | 1.10 | 4.05e-09 | 1.16 (1.03-1.32) | 1.74e-01 | no |
| putamen_basal_ganglia | M001 | G3BP2 | 1357.0 | 735.0 | 1.14 | 4.88e-09 | 1.57 (1.35-1.82) | 8.12e-06 | no |
| cerebellar_hemisphere | M000 | ESRP2 | 2700.0 | 1469.0 | 1.09 | 5.38e-09 | 1.07 (0.94-1.21) | 6.97e-01 | no |
| cerebellar_hemisphere | M000 | FXR1 | 2700.0 | 1102.0 | 1.11 | 9.69e-09 | 1.16 (1.02-1.32) | 2.00e-01 | no |
| cerebellar_hemisphere | M000 | IGF2BP2 | 2700.0 | 1604.0 | 1.08 | 1.13e-08 | 1.20 (1.06-1.36) | 6.91e-02 | no |
| substantia_nigra | M000 | RBM14 | 1657.0 | 1061.0 | 1.10 | 1.64e-08 | 1.20 (1.03-1.40) | 1.70e-01 | no |
| cerebellar_hemisphere | M000 | PABPC3 | 2700.0 | 1516.0 | 1.08 | 1.66e-08 | 1.23 (1.09-1.39) | 3.19e-02 | no |
| substantia_nigra | M000 | PPRC1 | 1657.0 | 877.0 | 1.11 | 7.62e-08 | 1.10 (0.95-1.29) | 5.92e-01 | no |
| cerebellar_hemisphere | M000 | AKAP1 | 2700.0 | 1451.0 | 1.08 | 1.15e-07 | 1.07 (0.95-1.21) | 6.57e-01 | no |

### Contradicted by opportunity adjustment

These 3 cells are hypergeometric-significant yet have an adjusted 95% CI entirely below 1 -- the raw enrichment is explained, and then some, by the motif opportunity their genes carry. They must not be described as candidate regulons.

| region | module | RBP | enrichment | q | adj. OR (95% CI) |
|---|---|---|---|---|---|
| cerebellar_hemisphere | M003 | HNRNPD | 1.24 | 4.12e-03 | 0.82 (0.68-1.00) |
| cerebellar_hemisphere | M003 | RNASEL | 1.28 | 2.17e-02 | 0.78 (0.61-1.00) |
| substantia_nigra | M004 | U2AF2 | 1.20 | 4.70e-02 | 0.78 (0.62-1.00) |

