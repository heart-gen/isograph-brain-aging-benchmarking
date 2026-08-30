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
- `01_synthetic_benchmark/` - synthetic benchmark runs and metrics.
- `real_data/` - BrainSEQ and GTEx aging analyses.
- `figures/` - manuscript-ready figure panels.
- `reports/` - statistical summaries and manifests.

## Synthetic Benchmark Execution

The synthetic benchmark uses a reduced paired grid. All methods share the same
`dataset_id` values so confidence intervals and method deltas can be computed
over paired synthetic datasets.

## Module Interpretation

Synthetic module interpretation accuracy is evaluated in `01_synthetic_benchmark/02_interpret/`
against ground-truth synthetic modules and switching genes.

Real-data module interpretation is run with:

```bash
sbatch real_data/brainseq/_h/04.interpret_modules.sh
sbatch real_data/gtex/_h/04.interpret_modules.sh
```

By default, real-data interpretation uses the GENCODE v47 primary-assembly GTF
at `/ocean/projects/bio250020p/shared/resources/genomes/human/gencode-v47/gtf/gencode.v47.primary_assembly.annotation.gtf`
to add structural transcript annotations for selected module genes. The GTF
parse cache is written under ignored `real_data/_m/tmp/`.

BrainSEQ `caudate_sczd` is diagnosis-focused: IsoGraph is fit with `Dx` as the
trait and downstream module associations are written to `diagnosis_assoc.parquet`.
