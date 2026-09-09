# Empirical null for the cross-cohort aging-replication count

`T_obs` counts matched BrainSEQ<->GTEx module pairs whose age association is nominally significant (p<0.05) in **both** cohorts with a concordant sign. The two nulls answer different questions and both are reported: **age** (Freedman-Lane) regenerates the GTEx age association under the age-null while holding module construction, the trusted sets and the gene-Jaccard matching fixed; **matching** holds every age statistic fixed and permutes which GTEx module each BrainSEQ module is matched to. `matching` is the stricter of the two.

The statistic is computed **the same way on both cohorts**. Scoring BrainSEQ with the published covariate-free Pearson while scoring GTEx with a covariate-adjusted model would make `both_sig` a hybrid of two age models and its p-value would calibrate a statistic that is never reported.

Covariate modes: `full` adds every covariate; `complement` adds only those the fit did not already residualize, so each is adjusted exactly once; `none` adds nothing. IsoGraph residualizes its DISCOVERY covariates inside the fit, so `complement` is the non-double-adjusting choice there — but only if the persisted feature_scores carry residualized values, which `covariate_leakage()` measures. WGCNA's eigengenes are not residualized at all, so it needs `full`.

## `isograph`

| covariates | statistic | T_obs | n pairs | null | null mean ± sd | p_emp |
|---|---|---|---|---|---|---|
| complement | pearson | 25 | 130 | age | 2.18 ± 3.39 | 0.0007999 |
| complement | pearson | 23 | 130 | matching | 16.08 ± 3.47 | 0.034 |
| complement | partial_linear | 24 | 130 | age | 2.01 ± 3.28 | 0.0013 |
| complement | partial_linear | 24 | 130 | matching | 15.35 ± 3.41 | 0.0123 |
| complement | spline_f | 15 | 130 | age | 1.91 ± 3.10 | 0.0103 |
| complement | spline_f | 15 | 130 | matching | 11.31 ± 3.01 | 0.1463 |
| none | pearson | 25 | 130 | age | 2.10 ± 3.29 | 0.0007999 |
| none | pearson | 25 | 130 | matching | 16.47 ± 3.50 | 0.0142 |
| none | partial_linear | 25 | 130 | age | 2.10 ± 3.29 | 0.0007999 |
| none | partial_linear | 25 | 130 | matching | 16.47 ± 3.50 | 0.0142 |
| none | spline_f | 16 | 130 | age | 2.04 ± 3.17 | 0.008699 |
| none | spline_f | 16 | 130 | matching | 12.65 ± 3.15 | 0.1825 |
| full | pearson | 25 | 130 | age | 1.82 ± 2.48 | 9.999e-05 |
| full | pearson | 25 | 130 | matching | 16.47 ± 3.50 | 0.0142 |
| full | partial_linear | 2 | 130 | age | 0.59 ± 1.27 | 0.1344 |
| full | partial_linear | 2 | 130 | matching | 2.25 ± 1.37 | 0.6892 |
| full | spline_f | 3 | 130 | age | 0.49 ± 1.10 | 0.05349 |
| full | spline_f | 3 | 130 | matching | 1.66 ± 1.18 | 0.2185 |

Statistics: `pearson` = Pearson r, no covariates (the published statistic); `partial_linear` = linear age term, covariate-adjusted; `spline_f` = df=3 natural cubic spline block F-test, covariate-adjusted

## `wgcna`

| covariates | statistic | T_obs | n pairs | null | null mean ± sd | p_emp |
|---|---|---|---|---|---|---|
| complement | pearson | 3 | 53 | age | 0.67 ± 1.53 | 0.1159 |
| complement | pearson | 3 | 53 | matching | 0.85 ± 0.85 | 0.0403 |
| complement | partial_linear | 2 | 53 | age | 0.65 ± 1.53 | 0.1704 |
| complement | partial_linear | 2 | 53 | matching | 0.71 ± 0.78 | 0.1473 |
| complement | spline_f | 7 | 53 | age | 0.64 ± 1.46 | 0.009399 |
| complement | spline_f | 7 | 53 | matching | 2.83 ± 1.46 | 0.008399 |
| none | pearson | 3 | 53 | age | 0.68 ± 1.53 | 0.1168 |
| none | pearson | 3 | 53 | matching | 0.85 ± 0.85 | 0.0403 |
| none | partial_linear | 3 | 53 | age | 0.68 ± 1.53 | 0.1168 |
| none | partial_linear | 3 | 53 | matching | 0.85 ± 0.85 | 0.0403 |
| none | spline_f | 8 | 53 | age | 0.71 ± 1.54 | 0.007499 |
| none | spline_f | 8 | 53 | matching | 3.11 ± 1.53 | 0.0037 |
| full | pearson | 3 | 53 | age | 0.34 ± 0.97 | 0.05369 |
| full | pearson | 3 | 53 | matching | 0.85 ± 0.85 | 0.0403 |
| full | partial_linear | 3 | 53 | age | 0.44 ± 1.12 | 0.06279 |
| full | partial_linear | 3 | 53 | matching | 1.55 ± 1.13 | 0.1912 |
| full | spline_f | 2 | 53 | age | 0.32 ± 0.99 | 0.06009 |
| full | spline_f | 2 | 53 | matching | 0.57 ± 0.63 | 0.07509 |

Statistics: `pearson` = Pearson r, no covariates (the published statistic); `partial_linear` = linear age term, covariate-adjusted; `spline_f` = df=3 natural cubic spline block F-test, covariate-adjusted

## How to write this up

The primary (covariate-adjusted spline) count is **15/130**, p_emp=0.1463 against the matching null. The published covariate-free Pearson count is **23/130**, p_emp=0.034.

**The count does not survive under the primary model.** Per the pre-registered decision rule the word "replication" must not be used; the honest wording is "matched modules with concordant age effects", reported alongside the covariate-free sensitivity analysis and the note that the two disagree.
