# Empirical null for the cross-cohort aging-replication count

`T_obs` counts matched BrainSEQ<->GTEx module pairs whose age association is nominally significant (p<0.05) in **both** cohorts with a concordant sign. The two nulls answer different questions and both are reported: **age** (Freedman-Lane) regenerates the GTEx age association under the age-null while holding module construction, the trusted sets and the gene-Jaccard matching fixed; **matching** holds every age statistic fixed and permutes which GTEx module each BrainSEQ module is matched to. `matching` is the stricter of the two.

The statistic is computed **the same way on both cohorts**. Scoring BrainSEQ with the published covariate-free Pearson while scoring GTEx with a covariate-adjusted model would make `both_sig` a hybrid of two age models and its p-value would calibrate a statistic that is never reported.

Covariate modes: `full` adds every covariate; `complement` adds only those the fit did not already residualize, so each is adjusted exactly once; `none` adds nothing. IsoGraph residualizes its DISCOVERY covariates inside the fit, so `complement` is the non-double-adjusting choice there — but only if the persisted feature_scores carry residualized values, which `covariate_leakage()` measures. WGCNA's eigengenes are not residualized at all, so it needs `full`.

## `isograph`

| covariates | statistic | T_obs | n pairs | null | null mean ± sd | p_emp |
|---|---|---|---|---|---|---|
| complement | pearson | 4 | 94 | age | 1.23 ± 2.59 | 0.1064 |
| complement | pearson | 4 | 94 | matching | 8.58 ± 2.64 | 0.9821 |
| complement | partial_linear | 2 | 94 | age | 1.10 ± 2.48 | 0.2153 |
| complement | partial_linear | 2 | 94 | matching | 7.94 ± 2.54 | 0.9985 |
| complement | spline_f | 2 | 94 | age | 1.08 ± 2.48 | 0.2035 |
| complement | spline_f | 2 | 94 | matching | 6.14 ± 2.26 | 0.9913 |
| none | pearson | 4 | 94 | age | 1.23 ± 2.57 | 0.1042 |
| none | pearson | 4 | 94 | matching | 8.58 ± 2.64 | 0.9821 |
| none | partial_linear | 4 | 94 | age | 1.23 ± 2.57 | 0.1042 |
| none | partial_linear | 4 | 94 | matching | 8.58 ± 2.64 | 0.9821 |
| none | spline_f | 1 | 94 | age | 1.13 ± 2.50 | 0.3659 |
| none | spline_f | 1 | 94 | matching | 6.38 ± 2.30 | 0.9993 |
| full | pearson | 4 | 94 | age | 0.66 ± 1.08 | 0.0247 |
| full | pearson | 4 | 94 | matching | 8.58 ± 2.64 | 0.9821 |
| full | partial_linear | 0 | 94 | age | 0.25 ± 0.63 | 1 |
| full | partial_linear | 0 | 94 | matching | 1.46 ± 1.09 | 1 |
| full | spline_f | 0 | 94 | age | 0.25 ± 0.59 | 1 |
| full | spline_f | 0 | 94 | matching | 1.63 ± 1.16 | 1 |

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

The primary (covariate-adjusted spline) count is **2/94**, p_emp=0.9913 against the matching null. The published covariate-free Pearson count is **4/94**, p_emp=0.9821.

**The count does not survive under the primary model.** Per the pre-registered decision rule the word "replication" must not be used; the honest wording is "matched modules with concordant age effects", reported alongside the covariate-free sensitivity analysis and the note that the two disagree.
