# Candidate RBP regulons among IsoGraph co-switch modules (intronic scope)

Motifs are scanned on **intronic splice-site flanks** (pre-mRNA sense; up to 100 nt into each intron), the binding niche for splicing-regulatory RBPs invisible to the mature-transcript scan.

For each module and RBP, whether binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool. Two arms are reported for every module x RBP cell:

1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, BH-corrected across all cells. It conditions on nothing else.
2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as covariates, so a module cannot score simply because its genes are long, GC-rich or UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the estimable cells only.

The adjusted arm is the stricter reading and the two can disagree in both directions; where they do, the disagreement is reported below rather than resolved in favour of the larger number.

- switch genes tested: **11859** over 14 regions
- module x RBP cells: **34240**
- hypergeometric q<0.05: **646** (GO-invisible 236)
- estimable GLM cells (BH family): **33934** of 34240
- covariate-adjusted q<0.05: **49**
- of the 646 hypergeometric hits: 646 estimable, **10** also adjusted-q<0.05, 273 with an adjusted CI entirely above 1, **0** entirely *below* it (enriched on raw counts, depleted once opportunity is adjusted for)

### Estimability

A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including such cells would put an arbitrary point estimate into the BH family and dilute every real test. They are excluded from the correction and reported here instead of silently carrying a p-value.

| GLM cell status | n | meaning |
|---|---|---|
| `fit` | 33934 | estimable; carries a p-value and enters the BH family |
| `separated_zero_cell` | 298 | an empty cell in the (in-module x switched) 2x2 -- complete separation, no finite MLE |
| `quasi_separated` | 8 | fitted probabilities pinned at 0/1, or a degenerate coefficient/standard error |

### Top candidate regulons (by hypergeometric q)

| region | module | RBP | module size | switched | enrichment | q | adj. OR (95% CI) | adj. q | GO-inv |
|---|---|---|---|---|---|---|---|---|---|
| frontal_cortex_ba9 | M007 | PPRC1 | 406.0 | 208.0 | 1.95 | 3.49e-24 | 1.83 (1.47-2.28) | 4.83e-04 | no |
| frontal_cortex_ba9 | M008 | ELAVL4 | 273.0 | 146.0 | 2.12 | 3.68e-20 | 1.37 (1.05-1.79) | 4.08e-01 | no |
| frontal_cortex_ba9 | M007 | RBM14 | 406.0 | 220.0 | 1.66 | 2.14e-16 | 1.93 (1.56-2.40) | 2.36e-05 | no |
| frontal_cortex_ba9 | M001 | RBMS3 | 677.0 | 242.0 | 1.63 | 1.40e-14 | 1.34 (1.12-1.61) | 1.77e-01 | no |
| frontal_cortex_ba9 | M001 | PPRC1 | 677.0 | 276.0 | 1.55 | 1.40e-14 | 1.23 (1.03-1.48) | 4.22e-01 | no |
| frontal_cortex_ba9 | M001 | HNRNPCL1 | 677.0 | 408.0 | 1.33 | 4.07e-13 | 1.33 (1.12-1.58) | 1.63e-01 | no |
| frontal_cortex_ba9 | M008 | TIAL1 | 273.0 | 138.0 | 1.82 | 1.70e-12 | 1.11 (0.85-1.45) | 8.61e-01 | no |
| frontal_cortex_ba9 | M008 | ELAVL1 | 273.0 | 107.0 | 2.07 | 5.16e-12 | 1.08 (0.82-1.42) | 9.05e-01 | no |
| frontal_cortex_ba9 | M001 | HNRNPA3 | 677.0 | 317.0 | 1.42 | 9.59e-12 | 1.28 (1.07-1.52) | 2.73e-01 | no |
| frontal_cortex_ba9 | M007 | RBM8A | 406.0 | 144.0 | 1.81 | 2.43e-11 | 1.52 (1.21-1.91) | 8.50e-02 | no |
| frontal_cortex_ba9 | M003 | PTBP1 | 514.0 | 90.0 | 2.29 | 2.63e-11 | 1.59 (1.20-2.11) | 1.68e-01 | yes |
| frontal_cortex_ba9 | M001 | RBM14 | 677.0 | 312.0 | 1.41 | 2.63e-11 | 1.27 (1.07-1.51) | 2.80e-01 | no |
| frontal_cortex_ba9 | M003 | HNRNPA1 | 514.0 | 112.0 | 2.03 | 4.62e-11 | 1.50 (1.17-1.94) | 1.82e-01 | yes |
| amygdala | M001 | PPRC1 | 412.0 | 185.0 | 1.48 | 7.86e-11 | 1.60 (1.22-2.09) | 1.15e-01 | yes |
| frontal_cortex_ba9 | M001 | G3BP1 | 677.0 | 329.0 | 1.37 | 1.17e-10 | 1.21 (1.02-1.44) | 4.46e-01 | no |
| frontal_cortex_ba9 | M008 | HNRNPC | 273.0 | 142.0 | 1.68 | 2.17e-10 | 1.10 (0.84-1.42) | 8.82e-01 | no |
| frontal_cortex_ba9 | M001 | RBM8A | 677.0 | 208.0 | 1.57 | 3.46e-10 | 1.26 (1.04-1.52) | 3.94e-01 | no |
| frontal_cortex_ba9 | M008 | TIA1 | 273.0 | 118.0 | 1.82 | 5.62e-10 | 1.16 (0.89-1.51) | 7.92e-01 | no |
| frontal_cortex_ba9 | M001 | RALY | 677.0 | 356.0 | 1.33 | 8.98e-10 | 1.24 (1.05-1.48) | 3.43e-01 | no |
| frontal_cortex_ba9 | M001 | SAMD4A | 677.0 | 276.0 | 1.42 | 1.31e-09 | 1.12 (0.94-1.34) | 7.42e-01 | no |
| frontal_cortex_ba9 | M008 | CELF1 | 273.0 | 103.0 | 1.93 | 1.36e-09 | 1.35 (1.02-1.77) | 4.71e-01 | no |
| frontal_cortex_ba9 | M007 | RALY | 406.0 | 229.0 | 1.42 | 2.44e-09 | 1.79 (1.44-2.23) | 7.54e-04 | no |
| frontal_cortex_ba9 | M011 | PPRC1 | 175.0 | 89.0 | 1.94 | 2.45e-09 | 1.35 (0.98-1.86) | 5.61e-01 | yes |
| frontal_cortex_ba9 | M008 | CELF2 | 273.0 | 106.0 | 1.86 | 3.93e-09 | 1.26 (0.96-1.66) | 6.10e-01 | no |
| frontal_cortex_ba9 | M003 | SSB | 514.0 | 134.0 | 1.76 | 3.93e-09 | 1.25 (0.99-1.60) | 5.56e-01 | yes |
| frontal_cortex_ba9 | M001 | CPEB2 | 677.0 | 389.0 | 1.28 | 4.41e-09 | 1.25 (1.06-1.49) | 3.19e-01 | no |
| frontal_cortex_ba9 | M001 | LIN28A | 677.0 | 404.0 | 1.26 | 9.51e-09 | 1.23 (1.04-1.47) | 3.83e-01 | no |
| frontal_cortex_ba9 | M001 | HNRNPA1L2 | 677.0 | 151.0 | 1.67 | 9.51e-09 | 1.40 (1.13-1.73) | 1.85e-01 | no |
| frontal_cortex_ba9 | M001 | RBM41 | 677.0 | 186.0 | 1.56 | 9.51e-09 | 1.15 (0.94-1.40) | 7.05e-01 | no |
| frontal_cortex_ba9 | M003 | CELF1 | 514.0 | 162.0 | 1.61 | 1.46e-08 | 1.17 (0.93-1.46) | 7.06e-01 | yes |
