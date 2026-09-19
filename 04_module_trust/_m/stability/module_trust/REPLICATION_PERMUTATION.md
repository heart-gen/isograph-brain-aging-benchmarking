# Empirical null for the cross-cohort aging-replication count

`T_obs` counts matched BrainSEQ<->GTEx module pairs whose age association is nominally significant (p<0.05) in **both** cohorts with a concordant sign. The two nulls answer different questions and both are reported: **age** (Freedman-Lane) regenerates the GTEx age association under the age-null while holding module construction, the trusted sets and the gene-Jaccard matching fixed; **matching** holds every age statistic fixed and permutes which GTEx module each BrainSEQ module is matched to. `matching` is the stricter of the two.

The statistic is computed **the same way on both cohorts**. Scoring BrainSEQ with the published covariate-free Pearson while scoring GTEx with a covariate-adjusted model would make `both_sig` a hybrid of two age models and its p-value would calibrate a statistic that is never reported.

Covariate modes: `full` adds every covariate; `complement` adds only those the fit did not already residualize, so each is adjusted exactly once; `none` adds nothing. IsoGraph residualizes its DISCOVERY covariates inside the fit, so `complement` is the non-double-adjusting choice there — but only if the persisted feature_scores carry residualized values, which `covariate_leakage()` measures. WGCNA's eigengenes are not residualized at all, so it needs `full`.

## `isograph`

| covariates | statistic | T_obs | n pairs | null | null mean ± sd | p_emp |
|---|---|---|---|---|---|---|
| complement | pearson | 1 | 38 | age | 0.56 ± 1.39 | 0.2495 |
| complement | pearson | 1 | 38 | matching | 2.24 ± 1.38 | 0.9116 |
| complement | partial_linear | 0 | 38 | age | 0.50 ± 1.34 | 1 |
| complement | partial_linear | 0 | 38 | matching | 2.64 ± 1.48 | 1 |
| complement | spline_f | 0 | 38 | age | 0.49 ± 1.29 | 1 |
| complement | spline_f | 0 | 38 | matching | 1.44 ± 1.11 | 1 |
| none | pearson | 1 | 38 | age | 0.54 ± 1.37 | 0.2462 |
| none | pearson | 1 | 38 | matching | 2.24 ± 1.38 | 0.9116 |
| none | partial_linear | 1 | 38 | age | 0.54 ± 1.37 | 0.2462 |
| none | partial_linear | 1 | 38 | matching | 2.24 ± 1.38 | 0.9116 |
| none | spline_f | 0 | 38 | age | 0.49 ± 1.29 | 1 |
| none | spline_f | 0 | 38 | matching | 1.94 ± 1.28 | 1 |
| full | pearson | 1 | 38 | age | 0.24 ± 0.60 | 0.1867 |
| full | pearson | 1 | 38 | matching | 2.24 ± 1.38 | 0.9116 |
| full | partial_linear | 0 | 38 | age | 0.10 ± 0.34 | 1 |
| full | partial_linear | 0 | 38 | matching | 0.21 ± 0.43 | 1 |
| full | spline_f | 0 | 38 | age | 0.08 ± 0.28 | 1 |
| full | spline_f | 0 | 38 | matching | 0.17 ± 0.40 | 1 |

Statistics: `pearson` = Pearson r, no covariates (the published statistic); `partial_linear` = linear age term, covariate-adjusted; `spline_f` = df=3 natural cubic spline block F-test, covariate-adjusted

## `wgcna`

| covariates | statistic | T_obs | n pairs | null | null mean ± sd | p_emp |
|---|---|---|---|---|---|---|
| complement | pearson | 2 | 52 | age | 0.59 ± 1.34 | 0.1587 |
| complement | pearson | 2 | 52 | matching | 0.72 ± 0.78 | 0.1506 |
| complement | partial_linear | 2 | 52 | age | 0.52 ± 1.20 | 0.1408 |
| complement | partial_linear | 2 | 52 | matching | 0.72 ± 0.78 | 0.1506 |
| complement | spline_f | 5 | 52 | age | 0.53 ± 1.24 | 0.0256 |
| complement | spline_f | 5 | 52 | matching | 3.08 ± 1.44 | 0.1611 |
| none | pearson | 2 | 52 | age | 0.58 ± 1.32 | 0.1543 |
| none | pearson | 2 | 52 | matching | 0.72 ± 0.78 | 0.1506 |
| none | partial_linear | 2 | 52 | age | 0.58 ± 1.32 | 0.1543 |
| none | partial_linear | 2 | 52 | matching | 0.72 ± 0.78 | 0.1506 |
| none | spline_f | 5 | 52 | age | 0.56 ± 1.24 | 0.0246 |
| none | spline_f | 5 | 52 | matching | 3.08 ± 1.44 | 0.1611 |
| full | pearson | 2 | 52 | age | 0.37 ± 0.79 | 0.07889 |
| full | pearson | 2 | 52 | matching | 0.72 ± 0.78 | 0.1506 |
| full | partial_linear | 3 | 52 | age | 0.30 ± 0.80 | 0.0391 |
| full | partial_linear | 3 | 52 | matching | 1.44 ± 1.07 | 0.1585 |
| full | spline_f | 2 | 52 | age | 0.21 ± 0.65 | 0.05049 |
| full | spline_f | 2 | 52 | matching | 1.24 ± 0.98 | 0.3647 |

Statistics: `pearson` = Pearson r, no covariates (the published statistic); `partial_linear` = linear age term, covariate-adjusted; `spline_f` = df=3 natural cubic spline block F-test, covariate-adjusted

## How to write this up

The primary (covariate-adjusted spline) count is **0/38**, p_emp=1 against the matching null. The published covariate-free Pearson count is **1/38**, p_emp=0.9116.

**The count does not survive under the primary model.** Per the pre-registered decision rule the word "replication" must not be used; the honest wording is "matched modules with concordant age effects", reported alongside the covariate-free sensitivity analysis and the note that the two disagree.
