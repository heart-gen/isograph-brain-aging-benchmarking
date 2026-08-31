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
