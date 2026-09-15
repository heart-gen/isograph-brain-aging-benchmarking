# Cross-cohort eigengene projection

Each trusted module's eigengene weights are frozen in the source cohort and applied unchanged to the target cohort's features. **Preservation** is the mean signed kME against a null that applies the same weights to size- and type-matched random target features (BH across modules). **Aging** is module-specific: in each cohort the frozen eigengene's age correlation is expressed as a z against the same-weight random projections, and the two cohorts' z are compared. Pooled over the three matched region pairs.

**Why the raw age correlation is not the statistic.** Both cohorts carry cohort-wide structure (RNA quality, ischemic time, composition) that is itself age-correlated, so any weighted projection inherits one age sign regardless of the module. On the raw correlations the three BrainSEQ→GTEx pairs gave sign agreement of 44/44, 1/50 and 32/36 — a property of each target cohort, not of the modules. The raw counts are kept in the last columns to show the confound.

| method | direction | modules | preserved q<0.05 | median signed kME (null) | median |r| vs native PC1 | age-z sign match | among source-age-sig | both sig | Spearman age z | raw sign match | raw both sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| isograph | brainseq_to_gtex | 94 | 77 (0.82) | 0.381 (0.250) | 0.958 | 76/94 (p = 5.9e-10) | 36/41 | 30 | 0.637 (p = 5.14e-12) | 57/94 | 60 |
| isograph | gtex_to_brainseq | 63 | 61 (0.97) | 0.337 (0.181) | 0.948 | 48/63 (p = 1.88e-05) | 20/27 | 16 | 0.447 (p = 0.000244) | 15/63 | 21 |
| wgcna | brainseq_to_gtex | 51 | 47 (0.92) | 0.717 (0.598) | 0.999 | 23/51 (p = 0.799) | 19/38 | 33 | 0.120 (p = 0.4) | 32/51 | 12 |
| wgcna | gtex_to_brainseq | 10 | 10 (1.00) | 0.687 (0.255) | 0.999 | 3/10 (p = 0.945) | 3/10 | 9 | 0.806 (p = 0.00486) | 4/10 | 2 |

## Switch-axis orientation

| pair | shared genes | orientable | flipped | median |cosine| |
|---|---|---|---|---|
| caudate | 11,476 | 5,739 | 2,616 | 0.500 |
| dlpfc_ba9 | 11,340 | 5,660 | 2,522 | 0.499 |
| hippocampus | 10,944 | 5,450 | 2,505 | 0.497 |

**Scope.** Granularity still matters: WGCNA's larger modules average more features, which raises kME stability mechanically, so compare each method with its own null rather than the two methods' raw preservation rates. No covariates enter the age test, matching the published linear arm; cohort and quantifier remain confounded.
