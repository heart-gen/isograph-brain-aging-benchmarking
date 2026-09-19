# Candidate RBP regulons among IsoGraph co-switch modules (intronic scope)

Motifs are scanned on **intronic splice-site flanks** (pre-mRNA sense; up to 100 nt into each intron), the binding niche for splicing-regulatory RBPs invisible to the mature-transcript scan.

For each module and RBP, whether binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool. Two arms are reported for every module x RBP cell:

1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, BH-corrected across all cells. It conditions on nothing else.
2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as covariates, so a module cannot score simply because its genes are long, GC-rich or UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the estimable cells only.

The adjusted arm is the stricter reading and the two can disagree in both directions; where they do, the disagreement is reported below rather than resolved in favour of the larger number.

- switch genes tested: **11790** over 16 regions
- module x RBP cells: **11040**
- hypergeometric q<0.05: **488** (GO-invisible 77)
- estimable GLM cells (BH family): **10514** of 11040
- covariate-adjusted q<0.05: **188**
- of the 488 hypergeometric hits: 488 estimable, **82** also adjusted-q<0.05, 247 with an adjusted CI entirely above 1, **0** entirely *below* it (enriched on raw counts, depleted once opportunity is adjusted for)

### Estimability

A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including such cells would put an arbitrary point estimate into the BH family and dilute every real test. They are excluded from the correction and reported here instead of silently carrying a p-value.

| GLM cell status | n | meaning |
|---|---|---|
| `fit` | 10514 | estimable; carries a p-value and enters the BH family |
| `outcome_or_predictor_constant` | 480 | no variation to model in the region |
| `separated_zero_cell` | 46 | an empty cell in the (in-module x switched) 2x2 -- complete separation, no finite MLE |

### Top candidate regulons (by hypergeometric q)

| region | module | RBP | module size | switched | enrichment | q | adj. OR (95% CI) | adj. q | GO-inv |
|---|---|---|---|---|---|---|---|---|---|
| cerebellar_hemisphere | M000 | G3BP1 | 2700.0 | 1158.0 | 1.19 | 1.13e-24 | 1.35 (1.18-1.54) | 9.45e-03 | no |
| cerebellar_hemisphere | M000 | PPRC1 | 2700.0 | 897.0 | 1.23 | 3.35e-24 | 1.27 (1.09-1.47) | 6.88e-02 | no |
| cerebellar_hemisphere | M000 | RBM14 | 2700.0 | 1127.0 | 1.17 | 5.14e-20 | 1.29 (1.13-1.48) | 2.27e-02 | no |
| cerebellar_hemisphere | M000 | SRSF11 | 2700.0 | 1176.0 | 1.17 | 5.23e-20 | 1.34 (1.18-1.53) | 9.45e-03 | no |
| cerebellar_hemisphere | M000 | CELF4 | 2700.0 | 1135.0 | 1.17 | 3.02e-19 | 1.34 (1.18-1.53) | 9.45e-03 | no |
| cerebellar_hemisphere | M000 | RBM8A | 2700.0 | 672.0 | 1.24 | 2.87e-18 | 1.30 (1.11-1.53) | 6.06e-02 | no |
| cerebellar_hemisphere | M000 | CELF5 | 2700.0 | 1106.0 | 1.16 | 1.25e-17 | 1.31 (1.15-1.50) | 1.37e-02 | no |
| cerebellar_hemisphere | M000 | RALY | 2700.0 | 1345.0 | 1.14 | 2.51e-17 | 1.34 (1.18-1.52) | 9.45e-03 | no |
| cerebellar_hemisphere | M000 | SRSF4 | 2700.0 | 1415.0 | 1.13 | 2.57e-17 | 1.42 (1.25-1.61) | 4.85e-04 | no |
| cerebellar_hemisphere | M000 | RBMS1 | 2700.0 | 455.0 | 1.31 | 3.52e-17 | 1.29 (1.06-1.57) | 1.82e-01 | no |
| cerebellar_hemisphere | M000 | RBMS3 | 2700.0 | 738.0 | 1.22 | 7.51e-17 | 1.24 (1.06-1.44) | 1.42e-01 | no |
| cerebellar_hemisphere | M000 | CPEB2 | 2700.0 | 1425.0 | 1.13 | 7.51e-17 | 1.30 (1.15-1.48) | 1.11e-02 | no |
| cerebellar_hemisphere | M000 | HNRNPM | 2700.0 | 1366.0 | 1.13 | 8.19e-16 | 1.23 (1.09-1.40) | 5.59e-02 | no |
| substantia_nigra | M000 | SAMD4A | 1657.0 | 660.0 | 1.22 | 8.31e-16 | 1.40 (1.19-1.65) | 1.11e-02 | no |
| cerebellar_hemisphere | M000 | HNRNPCL1 | 2700.0 | 1472.0 | 1.11 | 4.79e-15 | 1.20 (1.06-1.36) | 1.17e-01 | no |
| cerebellar_hemisphere | M000 | SNRPB2 | 2700.0 | 572.0 | 1.24 | 5.97e-15 | 1.32 (1.11-1.57) | 6.19e-02 | no |
| cerebellar_hemisphere | M000 | LIN28A | 2700.0 | 1469.0 | 1.11 | 1.41e-14 | 1.21 (1.06-1.37) | 9.77e-02 | no |
| cerebellar_hemisphere | M000 | PABPC3 | 2700.0 | 926.0 | 1.17 | 1.51e-14 | 1.30 (1.13-1.49) | 2.68e-02 | no |
| cerebellar_hemisphere | M000 | SAMD4A | 2700.0 | 985.0 | 1.16 | 1.51e-14 | 1.15 (1.00-1.32) | 3.25e-01 | no |
| cerebellar_hemisphere | M000 | ESRP2 | 2700.0 | 954.0 | 1.16 | 2.36e-14 | 1.17 (1.02-1.34) | 2.62e-01 | no |
| cerebellar_hemisphere | M000 | RBM28 | 2700.0 | 507.0 | 1.25 | 1.56e-13 | 1.48 (1.23-1.77) | 1.10e-02 | no |
| cerebellar_hemisphere | M000 | MSI1 | 2700.0 | 468.0 | 1.26 | 2.16e-13 | 1.34 (1.11-1.62) | 7.70e-02 | no |
| cerebellar_hemisphere | M000 | YTHDC1 | 2700.0 | 1071.0 | 1.14 | 4.23e-13 | 1.27 (1.11-1.45) | 3.39e-02 | no |
| cerebellar_hemisphere | M000 | ENOX1 | 2700.0 | 845.0 | 1.17 | 5.75e-13 | 1.28 (1.11-1.48) | 4.44e-02 | no |
| cerebellar_hemisphere | M000 | HNRNPA3 | 2700.0 | 1030.0 | 1.14 | 1.17e-12 | 1.15 (1.01-1.31) | 3.05e-01 | no |
| cerebellar_hemisphere | M000 | A1CF | 2700.0 | 824.0 | 1.17 | 1.19e-12 | 1.22 (1.05-1.40) | 1.43e-01 | no |
| cerebellar_hemisphere | M000 | CSTF2 | 2700.0 | 1505.0 | 1.10 | 1.32e-12 | 1.27 (1.13-1.44) | 1.97e-02 | no |
| cerebellar_hemisphere | M000 | RBM41 | 2700.0 | 635.0 | 1.20 | 3.71e-12 | 1.12 (0.95-1.31) | 5.81e-01 | no |
| cerebellar_hemisphere | M000 | PABPC5 | 2700.0 | 1079.0 | 1.13 | 4.46e-12 | 1.16 (1.02-1.32) | 2.65e-01 | no |
| cerebellar_hemisphere | M000 | SART3 | 2700.0 | 1239.0 | 1.12 | 7.21e-12 | 1.20 (1.06-1.36) | 1.16e-01 | no |
