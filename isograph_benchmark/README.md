# `isograph_benchmark` — benchmarking and brain-aging analysis package

`isograph_benchmark` is the Python package that drives this repository: it generates
the synthetic benchmark, scores module recovery and interpretation accuracy against
ground truth, runs the real-data BrainSEQ and GTEx aging analyses, and produces the
manuscript figures. It is the *driver* and *evaluation* layer; the method under
evaluation, **IsoGraph**, is a separate installable package (see
[The IsoGraph software](#the-isograph-software) below).

This README mirrors the scientific summary in `01_synthetic_benchmark/README.md` and adds the
software and pipeline details needed to run or extend the analyses.

## Package layout

| Module | Responsibility |
|---|---|
| `01_synthetic_benchmark/` | Synthetic data generation (`synthetic_data.py`), per-run execution (`run_one.py`), metric computation, and module interpretation (`interpret_modules.py`) |
| `stats/` | Bootstrap CIs and per-(scenario, method, metric) summaries (`summarize.py`); paired Wilcoxon tests with rank-biserial + Cliff's δ effect sizes and BH-FDR (`hypothesis_tests.py`) |
| `figures/` | Manuscript figure generation (`synthetic_benchmark.R`): Fig 1, Table 1, and supplementary panels figS1–figS11 |
| `real_data/` | BrainSEQ and GTEx model fits and downstream analysis (`run_models.py`, `incremental_association.py`, `module_enrichment.py`, `sweep_leiden.py`) |
| `inputs/` | Dataset-bundle construction from processed counts (`build_bundles.py`, `build_parquet.py`) |
| `gwas/`, `checks/` | GWAS-overlap utilities and pipeline integrity checks |
| `paths.py` | Repo-relative path resolution (`rel`, `ensure_dir`) used throughout |

The on-disk pipeline lives under `01_synthetic_benchmark/` in the repository root and runs in four
ordered stages:

```
00_design   → enumerate the run grid from configs/synthetic_grid.yaml
01_synthetic → generate datasets + run every method × scenario × seed cell
02_interpret → score switch-driver interpretation vs ground truth   (runs BEFORE metrics)
03_metrics  → collect, bootstrap-summarize, paired-test, and plot
```

`02_interpret` precedes `03_metrics` because the metrics stage's interpretation panel
(figS11) reads `02_interpret/_m/synthetic_interpret_summary.parquet`.

## What the benchmark measures

The benchmark evaluates IsoGraph's variational backend (`isograph_vae`) against a
gene-level WGCNA baseline (`wgcna_gene`) [@doi:10.1186/1471-2105-9-559], plus IsoGraph's
multiplex, residualizing, reliability, GPU, and linear variants where the scenario calls
for them. It spans **15 scenarios** and **12,530 paired runs** at **15–30 seeds per
cell**, on shared `dataset_id` values so every comparison is paired on identical data.
Metrics include module recovery, switch-gene detection rate, predicted-module count,
edge count, runtime, and (multiplex scenarios) per-gene channel-role recovery;
interpretation is scored by top-1 switch-driver transcript accuracy and switch-magnitude
Spearman correlation.

### Headline results

- **Switch-defined modules:** IsoGraph VAE beats WGCNA on module recovery when isoform
  switching is the dominant signal — e.g. `idealized_switching` 0.968 vs 0.572,
  `noise_stress` 0.974 vs 0.629, `non_switching_background` 0.841 vs 0.535 (all
  FDR < 0.05, large Cliff's δ).
- **Scale:** the two reach parity at BrainSEQ scale (16k genes × 300 samples: 0.625 vs
  0.625, n.s.).
- **Specificity:** under pure noise, IsoGraph returns ~0.05 recovery / ~9 modules vs
  WGCNA's spurious 0.405 / ~282 modules.
- **Limits (reported honestly):** WGCNA is stronger when modules are abundance-dominated
  (`abundance_switch_mixed` 0.675 vs VAE 0.331; multiplex closes to 0.492) and when 3′
  degradation corrupts the switch channel without an abundance fallback
  (`rna_degradation` 0.613 vs 0.466/0.472).
- **Confounds:** the residualizing variant restores or exceeds WGCNA under
  cell-composition, library-depth, and batch confounds.
- **Interpretation:** top-1 switch-driver transcript accuracy ≈ 1.00 (chance 0.25);
  switch-magnitude recovery is modest.

Full numbers, effect sizes, and the per-scenario narrative are in
[`01_synthetic_benchmark/README.md`](../benchmark/README.md).

## Running the pipeline

```bash
# Stage 00–01: build grid + run all cells (SLURM array jobs)
sbatch 01_synthetic_benchmark/00_design/_h/*.sh
sbatch 01_synthetic_benchmark/01_synthetic/_h/*.sh

# Stage 02: interpretation scoring
sbatch 01_synthetic_benchmark/02_interpret/_h/step_1.sh
sbatch 01_synthetic_benchmark/02_interpret/_h/step_2.sh

# Stage 03: collect → summarize (CIs + paired tests + effect sizes) → figures
sbatch 01_synthetic_benchmark/03_metrics/_h/step_1_collect.sh
python -m isograph_benchmark.stats.summarize
sbatch 01_synthetic_benchmark/03_metrics/_h/step_3_figures.sh
```

Environment on Bridges-2: the Python pipeline uses the `isograph` env at
`/ocean/projects/bio260021p/shared/opt/envs/isograph`; the R figure script uses the
R env at `/ocean/projects/bio250020p/shared/opt/env/R_env`. Set
`PYTHONPATH="$PWD:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src"` so the
driver can import both this package and IsoGraph.

## The IsoGraph software

The method under evaluation is **IsoGraph** (v0.1.4, Apache-2.0), developed in this
environment at `/ocean/projects/bio260021p/kbenjamin/software/IsoGraph` and published at
<https://github.com/heart-gen/IsoGraph>. IsoGraph discovers co-regulated transcript
programs from bulk RNA-seq by treating **gene abundance** and **isoform switching** as
separate channels of a multiplex network: starting from transcript-level counts it
builds gene-local switch coordinates from compositional transcript usage and a
standardized per-gene abundance signal, infers sparse gene-module structure via a VAE or
linear backend, classifies each module gene by the channel that drives its membership
(`switch_only`, `abundance_only`, `coupled`, `discordant`), and links modules to traits.

**Backends.** Five interchangeable backends — `baseline`, `latent`, `graph`,
`vae` (default), and `wgcna` — all support multiplex mode and the benchmark suite. The
VAE backend requires a separate PyTorch install (CPU/GPU builds are platform-specific);
the WGCNA backend requires R with the `WGCNA` package and `Rscript` on `PATH`.

**Key source modules** (`src/isograph/`):

- `features/` — channel construction: `switch.py` (compositional switch coordinates),
  `composition.py`, `residualize.py` (covariate residualization), `reliability.py`
  (degradation-aware edge weighting), `graph.py`.
- `models/` — `vae.py` (default backend), `multiplex.py` (abundance+switch fusion and
  α-threshold selection), `wgcna.py`, `latent.py`, `graph.py`, `baseline.py`.
- `explain/` — `isograph explain-module`: gene driver tables, transcript polarity,
  high-vs-low contrasts, PDF plots, and optional VAE decoder / Captum integrated-gradient
  attribution.
- `evaluation/` — metrics, selection, snapshot/regression tracking.
- `workflow/` — the `isograph` CLI (`fit`, `benchmark`, `explain-module`,
  `annotate-structure`, `freeze-real`) built on Hydra configs.

**Install:**

```bash
pip install isograph            # core (Python 3.11–3.14)
pip install torch               # for the VAE backend (match your CPU/GPU/CUDA)
# WGCNA backend: R with the WGCNA package + Rscript on PATH
```

IsoGraph's core runtime dependencies are `numpy`, `pandas`, `pyarrow`, `scipy`,
`scikit-learn`, `patsy`, `networkx`, `matplotlib`, `pydantic`, `hydra-core`, `mpmath`,
and `PyYAML` (optional extras: `captum` for Integrated-Gradients attribution,
`goatools` for GO enrichment, plus PyTorch and R/WGCNA for those backends).

This repository's analysis layer adds a small set of packages on top of IsoGraph,
pinned in `requirements-analysis.txt` at the repo root. The shared packages there track
IsoGraph's supported version ranges (e.g. `pyarrow>=17`, `pandas>=2.2`, `numpy>=1.26`,
`scipy>=1.13`, `scikit-learn>=1.5`, `patsy>=0.5`, `pyyaml>=6`); the only additions are
`statsmodels` (spline aging models), `pyhere` (repo-relative paths), and `pytest`/`ruff`
for development. Install into the same environment as IsoGraph:

```bash
pip install -r requirements-analysis.txt
```

**Fit your own bundle** (providing a `gene_counts` matrix activates the abundance
channel alongside the switch channel):

```bash
isograph fit --dataset-path path/to/bundle --backend vae \
  --output-dir artifacts/fits/my_dataset
```

Outputs: `modules.parquet`, `edges.parquet`, `traits.parquet`,
`feature_scores.parquet`, and `module_gene_roles.parquet`. See the IsoGraph repository
README, the [Read the Docs API](https://isograph.readthedocs.io/en/latest/), and the
[GitHub Wiki](https://github.com/heart-gen/IsoGraph/wiki) for full documentation, and
`CITATION.cff` for the software citation.

## Reproducibility notes

- The synthetic grid is fully specified by `configs/synthetic_grid.yaml` (the single
  source of truth); dataset hashes are stable, so adding scenarios does not perturb
  existing run IDs.
- Result tables are versioned parquet under each stage's `_m/` directory; bulky raw run
  outputs (`_o/`) and dataset bundles are git-ignored.
- Statistics are seed-controlled (base seed 13) with 10,000-resample bootstrap CIs and
  BH-FDR over the full comparison family.
