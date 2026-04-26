# IsoGraph Brain Aging Benchmarking

Benchmarking and brain aging analysis repository for IsoGraph.

Raw source files are in `inputs/raw/` for local reproducibility and are
intentionally ignored. They will be posted to Zenodo.

Terminology: synthetic nonlinear settings refer to interactions within feature space.
Real-data spline aging analyses refer to spline models of age against module eigengenes.

## Layout

- `inputs/raw/` - ignored local copies of BrainSEQ and GTEx source files.
- `inputs/processed/` - compressed parquet inputs used by analyses.
- `inputs/bundles/` - IsoGraph dataset bundles.
- `benchmark/` - synthetic benchmark runs and metrics.
- `real_data/` - BrainSEQ and GTEx aging analyses.
- `figures/` - manuscript-ready figure panels.
- `reports/` - statistical summaries and manifests.

## Synthetic Benchmark Execution

The synthetic benchmark uses a reduced paired grid. All methods share the same
`dataset_id` values so confidence intervals and method deltas can be computed
over paired synthetic datasets.
