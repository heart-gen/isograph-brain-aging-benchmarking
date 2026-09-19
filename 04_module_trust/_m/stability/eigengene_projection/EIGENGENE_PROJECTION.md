# Cross-cohort eigengene projection

Each trusted module's eigengene weights are frozen in the source cohort and applied unchanged to the target cohort's features. **Preservation** is the mean signed kME against a null that applies the same weights to size- and type-matched random target features (BH across modules). **Aging** is module-specific: in each cohort the frozen eigengene's age correlation is expressed as a z against the same-weight random projections, and the two cohorts' z are compared. Pooled over the three matched region pairs.

**Why the raw age correlation is not the statistic.** Both cohorts carry cohort-wide structure (RNA quality, ischemic time, composition) that is itself age-correlated, so any weighted projection inherits one age sign regardless of the module. On the raw correlations the three BrainSEQ→GTEx pairs gave sign agreement of 44/44, 1/50 and 32/36 — a property of each target cohort, not of the modules. The raw counts are kept in the last columns to show the confound.

| method | direction | modules | preserved q<0.05 | median signed kME (null) | median |r| vs native PC1 | age-z sign match | among source-age-sig | both sig | Spearman age z | raw sign match | raw both sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| isograph | brainseq_to_gtex | 38 | 37 (0.97) | 0.456 (0.269) | 0.981 | 33/38 (p = 2.13e-06) | 21/25 | 19 | 0.490 (p = 0.00181) | 22/38 | 29 |
| isograph | gtex_to_brainseq | 55 | 53 (0.96) | 0.332 (0.161) | 0.939 | 44/55 (p = 4.35e-06) | 18/22 | 17 | 0.555 (p = 1.09e-05) | 17/55 | 17 |
| wgcna | brainseq_to_gtex | 52 | 49 (0.94) | 0.733 (0.599) | 1.000 | 28/52 (p = 0.339) | 23/38 | 32 | 0.395 (p = 0.00379) | 28/52 | 8 |
| wgcna | gtex_to_brainseq | 10 | 10 (1.00) | 0.672 (0.253) | 0.999 | 3/10 (p = 0.945) | 3/9 | 8 | 0.685 (p = 0.0289) | 6/10 | 3 |

## Switch-axis orientation

| pair | shared genes | orientable | flipped | median |cosine| |
|---|---|---|---|---|
| caudate | 11,476 | 5,739 | 2,616 | 0.500 |
| dlpfc_ba9 | 11,340 | 5,660 | 2,522 | 0.499 |
| hippocampus | 10,944 | 5,450 | 2,505 | 0.497 |

**Scope.** Granularity still matters: WGCNA's larger modules average more features, which raises kME stability mechanically, so compare each method with its own null rather than the two methods' raw preservation rates. No covariates enter the age test, matching the published linear arm; cohort and quantifier remain confounded.
