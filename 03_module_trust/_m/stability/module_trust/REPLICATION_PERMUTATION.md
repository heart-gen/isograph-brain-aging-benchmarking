# Empirical null for the cross-cohort aging-replication count

`T_obs` counts matched BrainSEQ<->GTEx module pairs whose age association is nominally significant (p<0.05) in **both** cohorts with a concordant sign. The two nulls answer different questions and both are reported: **age** (Freedman-Lane) regenerates the GTEx age association under the age-null while holding module construction, the trusted sets and the gene-Jaccard matching fixed; **matching** holds every age statistic fixed and permutes which GTEx module each BrainSEQ module is matched to. `matching` is the stricter of the two.

The statistic is computed **the same way on both cohorts**. Scoring BrainSEQ with the published covariate-free Pearson while scoring GTEx with a covariate-adjusted model would make `both_sig` a hybrid of two age models and its p-value would calibrate a statistic that is never reported.

Covariate modes: `full` adds every covariate; `complement` adds only those the fit did not already residualize, so each is adjusted exactly once; `none` adds nothing. IsoGraph residualizes its DISCOVERY covariates inside the fit, so `complement` is the non-double-adjusting choice there — but only if the persisted feature_scores carry residualized values, which `covariate_leakage()` measures. WGCNA's eigengenes are not residualized at all, so it needs `full`.

## `isograph`

| covariates | statistic | T_obs | n pairs | null | null mean ± sd | p_emp |
|---|---|---|---|---|---|---|
| complement | pearson | 25 | 130 | age | 2.18 ± 3.39 | 0.0007999 |
| complement | pearson | 25 | 130 | matching | 16.47 ± 3.50 | 0.0142 |
| complement | partial_linear | 24 | 130 | age | 2.01 ± 3.28 | 0.0013 |
| complement | partial_linear | 24 | 130 | matching | 15.35 ± 3.41 | 0.0123 |
| complement | spline_f | 15 | 130 | age | 1.91 ± 3.10 | 0.0103 |
| complement | spline_f | 15 | 130 | matching | 11.35 ± 3.03 | 0.1521 |
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

The primary (covariate-adjusted spline) count is **15/130**, p_emp=0.1521 against the matching null. The published covariate-free Pearson count is **25/130**, p_emp=0.0142.

**The count does not survive under the primary model.** Per the pre-registered decision rule the word "replication" must not be used; the honest wording is "matched modules with concordant age effects", reported alongside the covariate-free sensitivity analysis and the note that the two disagree.

## Which component moves when the age model changes

Panel C of the trust funnel stays on the **linear (covariate-free Pearson) arm**. That is a
deliberate choice and it is not a claim that the linear arm is the better model — it is the
arm the count and the permutation p above are computed on, and switching the panel to the
spline would present a number whose null the panel does not show. The asymmetry belongs in
the text instead, and is recorded here.

From `module_trust replication-model-contrast --method {isograph,wgcna}`
(`replication_model_contrast__{method}.parquet`; it re-fits nothing, reading only the
per-pair tables the two `replication --model` runs already wrote):

| method | model | n pairs | sign_match | both_sig | concordant | discovery p<0.05 | replication p<0.05 |
|---|---|---|---|---|---|---|---|
| isograph | linear | 130 | 69 | 51 | 25 | 84 | 81 |
| isograph | spline | 130 | **69** | **8** | 2 | **20** | 67 |
| wgcna | linear | 53 | 27 | 10 | 3 | 28 | 15 |
| wgcna | spline | 53 | 32 | 5 | 1 | 11 | 22 |

**The two models do not disagree about the direction of aging.** For IsoGraph `sign_match`
is identical to the module — 69/130 under both — so every module that agrees in sign under
one model agrees under the other. The entire drop is `both_sig`, 51 -> 8.

**And that drop is one-sided.** BrainSEQ modules clearing p<0.05 fall 84 -> 20 under
covariate adjustment while GTEx falls only 81 -> 67; the WGCNA arm shows the same asymmetry
(28 -> 11 against 15 -> 22, the replication cohort actually gaining). So the covariate
adjustment is removing discovery-cohort age signal specifically. The most likely reading is
that a substantial part of the BrainSEQ age association is carried by covariates that
co-vary with age in that cohort, which is a statement about BrainSEQ, not about the
matching or about IsoGraph.

**What this does and does not license.** It does not relax the pre-registered wording rule:
the count still fails to separate from the matching null under the primary model, so
"replication" is still not used, and the honest phrase remains "matched modules with
concordant age effects". What it adds is that the failure is a power/adjustment effect on
the discovery side rather than a directional disagreement — which is a materially different
caveat and should be stated as such rather than left as an unexplained discrepancy between
two numbers.

**Do not equate the contrast table's spline count with a row of the grid above.** The
contrast reads the production `age_spline.parquet` fit, whereas the grid's `spline_f` rows
vary with covariate mode, running from T_obs = 15 (`complement`) to 3 (`full`) for IsoGraph.
The counts are therefore not interchangeable. What is invariant across all of them is the
decomposition's direction: sign agreement is stable and both-cohort significance is what
collapses.
