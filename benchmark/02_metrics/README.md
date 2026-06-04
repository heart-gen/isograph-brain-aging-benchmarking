---
title: "IsoGraph Synthetic Benchmark - Metric Analysis"
keywords:
  - isoform switching
  - coexpression networks
  - benchmarking
  - brain aging
  - IsoGraph
---

## Overview

This directory contains the downstream metric analysis for the IsoGraph synthetic benchmark described in `benchmark/01_synthetic/`.
The benchmark evaluates IsoGraph and WGCNA [@doi:10.1186/1471-2105-9-559] across six parameterized synthetic scenarios designed to stress-test network recovery under conditions relevant to brain-aging transcriptomics.
Figures and tables here are intended for inclusion in the IsoGraph manuscript targeting *Nature Methods*.
The collected result table includes CPU, SLURM, memory, runtime, software, and hardware metadata for every completed run.
Scale and scale_realistic scenarios are complete; `abundance_switch_mixed` and `isograph_vae_multiplex` runs are pending.

## Synthetic Data Generation

Synthetic datasets were generated using a negative binomial count model with ground-truth gene modules and isoform-switching signals.
Each non-scale dataset contains 400 genes (800 transcripts) across 160 samples with simulated covariates (age, sex, diagnosis).
Switching genes are assigned to latent modules correlated with age and disease status; the proportion of switching genes (`switching_fraction`) and expression noise (`noise_sd`) are varied across scenarios.

The `scale` scenario uses 1,000–12,000 genes across 240 samples; `scale_realistic` uses 16,000 genes across 300 samples (matching BrainSEQ dimensions).
Datasets were generated with a deterministic seed sequence (base seed 13) to ensure reproducibility.
Each dataset is identified by a SHA1 hash of its scenario parameters and seed, ensuring stable run IDs across grid changes.

## Benchmark Design

Seven scenarios are evaluated across the applicable method sets.
All non-scale, non-multiplex scenarios are run for five CPU methods.
The scale scenarios are run for IsoGraph VAE and WGCNA only (other methods are memory-prohibitive at 16k genes).
The `abundance_switch_mixed` scenario uses the multiplex method set.

| Scenario | Key parameters | Seeds | Methods |
|---|---|---:|---:|
| **Idealized switching** | `switching_fraction` in {0.10, 0.25, 0.50, 0.75}; `noise_sd` in {0.10, 0.25, 0.40} | 30 | 7 (excl. multiplex) |
| **Noise stress** | `count_dispersion` in {3, 7, 15, 30}; `noise_sd` in {0.05, 0.10, 0.25, 0.50} | 20 | 7 |
| **Feature space interactions** | `interaction_strength` in {0.0, 0.5, 1.0, 2.0, 3.0}; `interaction_fraction` in {0.25, 0.50, 0.75} | 20 | 7 |
| **Non-switching background** | `switching_fraction` in {0.05, 0.10, 0.25, 0.50} | 20 | 7 |
| **Unequal isoform abundance** | `abundance_imbalance` in {1, 4, 10, 25} | 20 | 7 |
| **Abundance switch mixed** | `abundance_fraction` in {0.0, 0.3, 0.5, 0.7, 1.0}; `n_genes` in {200, 500} | 20 | 3 (multiplex set) |
| **Scale** | `n_genes` in {1000, 3000, 6000, 12000}; `switching_fraction` in {0.15, 0.25} | 15 | 2 |
| **Scale realistic** | `n_genes` = 16000; `n_samples` = 300; `switching_fraction` in {0.15, 0.25} | 15 | 2 |

**Methods evaluated:**

| Method | Description | Compute | Scenarios |
|---|---|---|---|
| `isograph_baseline` | IsoGraph baseline network model (`alpha = 0.10`) | CPU | all except scale, multiplex |
| `isograph_latent` | IsoGraph latent-space model (3-fold CV for components, `alpha = 0.10`) | CPU | all except scale, multiplex |
| `isograph_graph` | IsoGraph graph-regularized model (3-fold CV, `alpha = 0.10`) | CPU | all except scale, multiplex |
| `isograph_cpu_latent` | IsoGraph latent-space model (LatentNetworkModel, `alpha = 0.10`) | CPU | all except scale, multiplex |
| `isograph_vae` | Variational autoencoder network model | CPU | all |
| `isograph_vae_multiplex` | IsoGraph VAE with multiplex (abundance + switch) features | CPU | multiplex only |
| `wgcna_gene` | WGCNA gene-level coexpression [@doi:10.1186/1471-2105-9-559] | CPU | all |

Total planned runs: **9,480** across 7 scenarios.

## Run Status

Scale and scale_realistic scenarios are fully complete. The `abundance_switch_mixed` and `isograph_vae_multiplex` runs are pending.

| Scenario | Methods | Planned | Completed |
|---|---|---:|---:|
| Idealized switching | 7 (excl. multiplex) | 2,520 | 2,520 |
| Noise stress | 7 | 2,240 | 2,240 |
| Feature space interactions | 7 | 2,100 | 1,800 |
| Non-switching background | 7 | 560 | 560 |
| Unequal isoform abundance | 7 | 560 | 560 |
| Abundance switch mixed | 3 | 1,200 | 0 |
| Scale | 2 (VAE + WGCNA) | 240 | 240 |
| Scale realistic | 2 (VAE + WGCNA) | 60 | 60 |
| **Total** | | **9,480** | **7,980** |

**Pending runs (2,340):**
- `isograph_vae_multiplex` across all non-scale, non-multiplex scenarios (1,140 runs)
- `abundance_switch_mixed` for all three multiplex methods (1,200 runs)

Current figures and tables use the **7,140 completed runs** in the current grid that cover the five non-multiplex CPU methods and both scale scenarios.

## Results

See generated figures in `figures/`, the main summary table `_m/table1_benchmark_summary.csv`, and the scale compute table `_m/tableS_scale_compute_summary.csv`.
All figures were regenerated with R, ggpubr, ggplot2, and patchwork.

**Primary metrics:**

- **Module recovery** - Area-under-curve score comparing predicted module assignments to ground-truth modules; higher is better (range 0 to 1).
- **Switch gene detection rate** - Recall of true isoform-switching genes among predicted module members; higher is better (range 0 to 1).
- **Non-switching gene module rate** - Fraction of non-switching genes incorrectly assigned to modules; lower is better (range 0 to 1). Main figures show the transformed specificity measure `1 - non-switching gene module rate` so higher is better.
- **Runtime** - Wall-clock time per run in seconds; runtime panels use log scaling where needed for readability.
- **Resource metadata** - CPU/GPU backend, requested CPU and GPU counts, requested memory, observed maximum resident memory, torch CUDA availability, and GPU peak memory where available.

**Key findings (from `figures/fig1_benchmark_overview.pdf`):**

Across the five non-scale scenarios, CPU IsoGraph VAE showed the strongest module recovery while maintaining complete switch-gene detection.
219 of 240 paired Wilcoxon tests (FDR < 0.05) favored IsoGraph over WGCNA.
The graph-regularized IsoGraph model improved module recovery relative to the regular latent variant under several stress conditions.
WGCNA retained high switch-gene detection but showed higher non-switching gene module rates in several scenarios.
See supplementary figures S1-S5 for parameter-resolved analyses and Figure S6 for the scale and scale_realistic compute comparison.

In the scale scenarios, IsoGraph VAE and WGCNA are compared for runtime and peak host memory across 1,000–16,000 genes.

## Methods

### Synthetic benchmark framework

Synthetic transcriptomic datasets were generated by simulation as described in `isograph_benchmark/benchmark/synthetic_data.py`.
Each gene was modeled with two isoforms; a proportion of genes (`switching_fraction`) were designated as switching genes whose predominant isoform varies with a latent factor correlated with age and diagnostic status.
Expression counts were drawn from a negative binomial distribution (mean parameterized by isoform proportions, dispersion = 3 unless otherwise specified).
Isoform inclusion ratios (PSI) were computed from simulated counts.

### Methods evaluated

**IsoGraph** methods use gene-level and transcript-level coexpression signals derived from PSI and count features.
Baseline, Latent, Graph, and cpu_latent models used significance threshold `alpha = 0.10`.
Latent and Graph models selected the number of latent components by 3-fold cross-validation where applicable.
VAE models selected latent dimension from grid {2, 4, 6, 8, 12}, used hidden dimension 128 or 256 depending on gene count, trained for up to 300 epochs with patience = 35, and used significance threshold `alpha = 0.70`.
**WGCNA** was run with signed network topology; soft-thresholding power selected by scale-free fit (R2 target >= 0.85); merge cut height = 0.25.

### Statistical analysis

All metric summaries report the mean and 95% bootstrap confidence interval (10,000 resampling iterations per group).
Groups are defined by scenario x method for the main summaries and by scenario x method x parameter combination in parameter-resolved supplemental panels.
Multiple comparison correction uses the Benjamini-Hochberg procedure [@doi:10.1111/j.2517-6161.1995.tb02031.x] where applicable.
Figure significance annotations use `ggpubr::geom_pwc()` with Wilcoxon tests and Benjamini-Hochberg adjusted p-values.
Main benchmark and parameter-sweep annotations compare IsoGraph VAE with WGCNA; scale compute annotations compare VAE GPU and WGCNA with VAE CPU.

### Compute environment

Jobs were submitted to the Pittsburgh Supercomputing Center (PSC) Bridges-2 cluster [@doi:10.1145/3437359.3465593].
All jobs used the RM-shared partition with 4 CPU cores per task and 2 GB memory per CPU (8 GB total) for non-scale runs.
Scale and scale_realistic jobs used 8 CPU cores (16 GB total) with a 4-hour time limit.

Observed CPU hardware included AMD EPYC 7742, Intel Xeon Gold 6248, and Intel Xeon Platinum 8470 nodes.
The IsoGraph conda environment is at `/ocean/projects/bio260021p/shared/opt/envs/isograph`.
The Python figure-generation environment used by the metric script is the same IsoGraph conda environment.
The R figure-generation environment (alternative) is at `/ocean/projects/bio250020p/shared/opt/env/R_env`.

## Reproduction

All outputs in this directory can be regenerated from the raw run outputs in `benchmark/01_synthetic/_o/runs/`.
Scripts in `_h/` can be run directly on a login node or submitted via `sbatch` for cluster execution.

SLURM opens the output log path before the script body runs, so create the metrics log directory before submitting jobs:

```bash
mkdir -p benchmark/02_metrics/_m/logs
```

```bash
# Step 1 - Collect all telemetry into a single parquet
bash benchmark/02_metrics/_h/step_1_collect.sh
# or: sbatch benchmark/02_metrics/_h/step_1_collect.sh

# Step 2 - Compute bootstrap CI summaries and paired Wilcoxon tests (BH FDR)
bash benchmark/02_metrics/_h/step_2_summarize.sh
# or: sbatch benchmark/02_metrics/_h/step_2_summarize.sh

# Step 3 - Generate figures and Table 1 (Python/matplotlib or R/ggplot2)
bash benchmark/02_metrics/_h/step_3_figures.sh
# or: sbatch benchmark/02_metrics/_h/step_3_figures.sh
```

Step 3 has two implementations:
- `isograph_benchmark/figures/synthetic_benchmark.py` — primary Python/matplotlib figures (used by `step_3_figures.sh`)
- `isograph_benchmark/figures/synthetic_benchmark.R` — alternative R/ggplot2 figures (requires `R_env`)

All three steps complete in under 30 minutes on a login node with at least 4 CPU cores and 8 GB RAM.
When run through SLURM, the scripts write logs as `collect-<jobid>.log`, `summarize-<jobid>.log`, and `figures-<jobid>.log` in `_m/logs/`.

## Files

| File | Description |
|---|---|
| `_h/step_1_collect.sh` | Collects `telemetry.json` files into `benchmark/01_synthetic/_m/synthetic_results.parquet` |
| `_h/step_2_summarize.sh` | Computes bootstrap CI summaries into `_m/synthetic_metric_summary.parquet` |
| `_h/step_3_figures.sh` | Generates all R/ggplot2 figures and summary tables |
| `_m/synthetic_results.parquet` | Raw collected telemetry; 9,617 rows including historical and current runs |
| `_m/synthetic_metric_summary.parquet` | Bootstrap summary by scenario, method, and metric; 282 rows |
| `_m/synthetic_metric_long.parquet` | Per-run long-format metrics; 57,702 rows |
| `_m/synthetic_pairwise_tests.parquet` | Paired Wilcoxon signed-rank tests vs. WGCNA with BH FDR; 240 tests |
| `_m/table1_benchmark_summary.csv` | Manuscript-ready CPU benchmark summary table with significance markers |
| `_m/tableS_scale_compute_summary.csv` | Scale compute table for IsoGraph VAE and WGCNA across gene-count grid |
| `_m/logs/` | Reproducibility logs for collection, summary, and figure/table generation |
| `figures/fig1_benchmark_overview.pdf` | Main CPU accuracy and specificity benchmark overview |
| `figures/figS1_idealized_switching.pdf` | Response dot plot: switching fraction by noise SD |
| `figures/figS2_noise_stress.pdf` | Response dot plot: count dispersion by noise SD |
| `figures/figS3_feature_interactions.pdf` | Response dot plot: interaction strength by interaction fraction |
| `figures/figS4_nonswitching_background.pdf` | Module recovery and non-switching specificity vs. switching fraction |
| `figures/figS5_unequal_abundance.pdf` | Module recovery and switch detection vs. isoform abundance imbalance |
| `figures/figS6_scale_compute.pdf` | Runtime and peak host RAM vs. number of genes (scale scenario) |

PNG versions of all figures are written alongside the PDF files for quick inspection.

## References

<!-- References auto-resolved by manubot from [@doi:...] tags above -->
