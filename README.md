# IsoGraph Brain Aging Benchmarking

Self-contained benchmarking and brain aging analysis repository for IsoGraph.

Raw source files are copied into `inputs/raw/` for local reproducibility, but that
directory is intentionally ignored and must never be tracked. Downstream steps convert
raw text/gzip inputs into compressed parquet under `inputs/processed/` and IsoGraph
dataset bundles under `inputs/bundles/`.

## Layout

- `inputs/raw/` - ignored local copies of BrainSEQ and GTEx source files.
- `inputs/processed/` - compressed parquet inputs used by analyses.
- `inputs/bundles/` - IsoGraph dataset bundles.
- `benchmark/` - synthetic benchmark runs and metrics.
- `real_data/` - BrainSEQ and GTEx aging analyses.
- `figures/` - manuscript-ready figure panels.
- `reports/` - statistical summaries and manifests.

## Data Policy

Run this before committing:

```bash
python -m isograph_benchmark.checks.no_tracked_raw
```

It fails if any file under `inputs/raw/` is tracked by git.

## Main Commands

```bash
python -m isograph_benchmark.inputs.copy_raw
python -m isograph_benchmark.inputs.build_parquet
python -m isograph_benchmark.inputs.build_bundles
python -m isograph_benchmark.benchmark.run_synthetic
python -m isograph_benchmark.real_data.run_models
python -m isograph_benchmark.stats.summarize
```

