# Module-level cis control of isoform choice

Does the allele-specific signal concentrate in the co-switching modules that carry
the aging association? Aggregation of the fitted within-donor allelic test; nothing
is refitted and no module eigengene is mapped.

Gene counts as cis-controlled if any fitted pair clears q < 0.05.
**Primary arm:** every module with a fitted gene is kept.
Module age strength is ranked within (region, trait); null is 10,000 permutations of the gene -> module labels (seed 13), which hold module
sizes and the number of cis-controlled genes fixed.

## Per region

| region | modules | ranked | genes | cis genes | rate range | rho | null mean | perm p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| caudate | 7 | 7 | 481 | 185 | 0.30–0.50 | -0.691 | +0.007 | 0.070 |
| dlpfc | 3 | 2 | — | — | — | — | — | *too few modules in a multi-module trait stratum to correlate* |
| hippocampus | 15 | 14 | 279 | 89 | 0.00–1.00 | +0.326 | +0.004 | 0.286 |

## Pooled

Mean rho over caudate, hippocampus: **-0.183** against a null mean of +0.005 (95% -0.460 to +0.473), permutation p = **0.470**.

## Modules

| region | module | trait | age effect | genes fitted | cis | rate | aging rank |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| caudate | M000 | Age_linear | +0.26 | 84 | 29 | 0.35 | 1.00 |
| caudate | M002 | Age_linear | +0.15 | 92 | 34 | 0.37 | 0.50 |
| caudate | M003 | Age_spline | +2.94 | 100 | 30 | 0.30 | 1.00 |
| caudate | M004 | Age_linear | +0.15 | 75 | 30 | 0.40 | 0.00 |
| caudate | M006 | Age_spline | +0.84 | 80 | 40 | 0.50 | 0.00 |
| caudate | M008 | Age_spline | +1.67 | 38 | 17 | 0.45 | 0.67 |
| caudate | M011 | Age_spline | +1.42 | 12 | 5 | 0.42 | 0.33 |
| dlpfc | M000 | Age_spline | -2.05 | 58 | 20 | 0.34 | nan |
| dlpfc | M007 | Age_spline | -1.51 | 32 | 11 | 0.34 | nan |
| dlpfc | M008 | Age_linear | +0.23 | 25 | 6 | 0.24 | nan |
| hippocampus | M000 | Age_linear | +0.27 | 28 | 12 | 0.43 | 0.15 |
| hippocampus | M001 | Age_linear | +0.44 | 44 | 16 | 0.36 | 0.92 |
| hippocampus | M002 | Age_linear | +0.35 | 28 | 6 | 0.21 | 0.46 |
| hippocampus | M003 | Age_linear | +0.41 | 13 | 2 | 0.15 | 0.77 |
| hippocampus | M004 | Age_linear | +0.42 | 77 | 22 | 0.29 | 0.85 |
| hippocampus | M006 | Age_linear | +0.17 | 7 | 1 | 0.14 | 0.00 |
| hippocampus | M007 | Age_linear | +0.35 | 31 | 6 | 0.19 | 0.54 |
| hippocampus | M008 | Age_linear | +0.39 | 22 | 9 | 0.41 | 0.69 |
| hippocampus | M009 | Age_linear | +0.30 | 12 | 5 | 0.42 | 0.38 |
| hippocampus | M010 | Age_linear | +0.28 | 3 | 3 | 1.00 | 0.23 |
| hippocampus | M011 | Age_linear | +0.39 | 8 | 6 | 0.75 | 0.62 |
| hippocampus | M012 | Age_linear | +0.28 | 4 | 0 | 0.00 | 0.31 |
| hippocampus | M013 | Age_linear | +0.25 | 1 | 0 | 0.00 | 0.08 |
| hippocampus | M015 | Age_linear | +0.45 | 1 | 1 | 1.00 | 1.00 |

**Reading.** A positive rho means the modules with the stronger age association also
carry proportionally more cis-controlled genes. With 3-15 modules per region this arm
is descriptive: the permutation null is the only honest reference, and an interval
spanning zero means the data do not separate the two possibilities.
