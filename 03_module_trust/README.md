# 03 — Module trust

**Question:** are the discovered modules reproducible, or is a fine-grained partition
just noise?

Two arms that answer one question, so they share a stage: within-cohort split-half
stability, and cross-cohort (BrainSEQ ↔ GTEx) aging replication. BrainSEQ and GTEx are
independent cohorts — different donors, libraries and quantifiers — run through the same
pipeline, which is what makes the second arm a real replication.

## Order

| Step | Wrapper | Produces |
|---|---|---|
| 01–03 | `stability_isograph`, `stability_wgcna`, `stability_aggregate` | Split-half partitions + `stability_summary` |
| 04 | `module_meta` | Per-module meta tables, driver loadings |
| 05–07 | `lr_validation`, `lr_validation_launch`, `lr_aggregate` | Learning-rate / software-robustness validation |
| 08 | `gcap_ab` | Giant-cap ablation (`_m/stability_gcap_ab/`) |
| 09 | `module_trust_replication` | Q3 cross-cohort aging replication of trusted modules |
| 10–11 | `replication`, `replication_go` | Cross-cohort module matching + GO consistency |
| 12–13 | `replication_permutation`, `replication_functional` | Permutation null + functional preservation |

`_m/` holds `stability/` and `replication/`, plus the resolution-2.0 and giant-cap
variants as named siblings (`stability_res2`, `stability_gcap_ab`, `replication_res2`).

**CLIs:** `isograph_benchmark/real_data/{stability,module_trust,replication,replication_go,replication_permutation,replication_functional}.py`.

## Caveat

Cross-cohort **module** Jaccard is granularity-confounded — a method with few giant
modules wins it mechanically — so it is retired as a quality gate in favour of
gene-level comparison. See `docs/MODULE_TRUST_PLAN.md`.

## Display items

Main **Fig 2** `figTrustFunnel`; supplementary `tableS7_module_trust_funnel`.
Design notes: `docs/MODULE_TRUST_PLAN.md`, `docs/SOFTWARE_ROBUSTNESS_PLAN.md`.

## Age model: why the linear arm is quoted

The cross-cohort aging replication count is **25/130 matched module pairs, `complement`
covariate mode, linear age model** (matching permutation p = 0.014). The df = 3 spline
arm gives 15/130 (p = 0.15) at the same covariate mode.

That gap is a power cost, not a correction, and `14.age_model_curvature.sh` is the test
that says so. It applies a Wald second-difference contrast inside the spline's own fit —
is each module's projected age trajectory distinguishable from a straight line? Across all
266 modules in all six cohort x region fits, none is (0 at nominal p < 0.05 against 13.3
expected by chance; smallest p = 0.71). The contrast is well powered because the projected
effects share nearly all their uncertainty and every contrast row sums to zero, so the
shared part cancels. The spline also finds no *different* modules: 4/266 spline-only
significant against 78/266 linear-only, and 4 is the chance rate.

So aging is linear at module resolution here, the spline spends 2 df on noise, and the
spline arm belongs in the supplement as a declared sensitivity.

The test needs the full projection covariance (`cov_age_ij`) and refuses to run on the
diagonal alone — dropping the off-diagonal terms would inflate the variance by exactly the
shared term the contrast cancels, and so would manufacture a null. The WGCNA baseline's
`age_spline.parquet` carries only per-point standard errors, so it cannot be tested here.

    sbatch 03_module_trust/_h/14.age_model_curvature.sh
    # -> _m/stability/module_trust/age_model_curvature__isograph.{parquet,json}
