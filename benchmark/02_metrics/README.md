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
The analysis was refreshed after completion of all benchmark jobs, and the collected result table now includes CPU, GPU, SLURM, torch, memory, runtime, software, and hardware metadata for every completed run.

## Synthetic Data Generation

Synthetic datasets were generated using a negative binomial count model with ground-truth gene modules and isoform-switching signals.
Each non-scale dataset contains 400 genes (800 transcripts) across 160 samples with simulated covariates (age, sex, diagnosis).
Switching genes are assigned to latent modules correlated with age and disease status; the proportion of switching genes (`switching_fraction`) and expression noise (`noise_sd`) are varied across scenarios.

The scale scenario uses larger synthetic datasets with 1,000, 3,000, 6,000, or 12,000 genes across 240 samples.
Datasets were generated with a deterministic seed sequence (base seed 13) to ensure reproducibility.
All 1,204 unique synthetic datasets are stored in `benchmark/01_synthetic/_m/datasets/`.

## Benchmark Design

Six scenarios were evaluated across the applicable method sets.
All non-scale scenarios were run for all eight methods.
The scale scenario was run for the scalable CPU/GPU IsoGraph variants and WGCNA.

| Scenario | Key parameters | Seeds | Runs per included method |
|---|---|---:|---:|
| **Idealized switching** | `switching_fraction` in {0.10, 0.25, 0.50, 0.75}; `noise_sd` in {0.10, 0.25, 0.40} | 30 | 360 |
| **Noise stress** | `count_dispersion` in {3, 7, 15, 30}; `noise_sd` in {0.05, 0.10, 0.25, 0.50} | 20 | 320 |
| **Feature space interactions** | `interaction_strength` in {0.0, 0.5, 1.0, 2.0, 3.0}; `interaction_fraction` in {0.25, 0.50, 0.75} | 20 | 300 |
| **Non-switching background** | `switching_fraction` in {0.05, 0.10, 0.25, 0.50} | 20 | 80 |
| **Unequal isoform abundance** | `abundance_imbalance` in {1, 4, 10, 25} | 20 | 80 |
| **Scale** | `n_genes` in {1000, 3000, 6000, 12000}; `switching_fraction` in {0.15, 0.25} | 8 | 64 |

**Methods evaluated:**

| Method | Description | Compute |
|---|---|---|
| `isograph_baseline` | IsoGraph baseline network model (`alpha = 0.10`) | CPU |
| `isograph_latent` | IsoGraph latent-space model (3-fold CV for components, `alpha = 0.10`) | CPU |
| `isograph_graph` | IsoGraph graph-regularized model (3-fold CV, `alpha = 0.10`) | CPU |
| `isograph_cpu_latent` | GPU-architecture latent model run on CPU | CPU |
| `isograph_gpu_latent` | GPU-accelerated latent model | GPU |
| `isograph_vae` | Variational autoencoder network model | CPU |
| `isograph_vae_gpu` | GPU-accelerated VAE model | GPU |
| `wgcna_gene` | WGCNA gene-level coexpression [@doi:10.1186/1471-2105-9-559] | CPU |

Total planned and completed runs: **9,440** across 1,204 unique datasets.

## Run Status

All planned runs completed successfully and were included in the refreshed collection step.

| Method | Completed runs | Status |
|---|---:|---|
| isograph_baseline | 1,140 | Complete |
| isograph_latent | 1,140 | Complete |
| isograph_graph | 1,140 | Complete |
| isograph_cpu_latent | 1,204 | Complete |
| isograph_gpu_latent | 1,204 | Complete |
| isograph_vae | 1,204 | Complete |
| isograph_vae_gpu | 1,204 | Complete |
| wgcna_gene | 1,204 | Complete |
| **Total completed** | **9,440** | **Complete** |

**Scenario coverage:**

| Scenario | Methods included | Runs |
|---|---:|---:|
| Idealized switching | 8 | 2,880 |
| Noise stress | 8 | 2,560 |
| Feature space interactions | 8 | 2,400 |
| Non-switching background | 8 | 640 |
| Unequal isoform abundance | 8 | 640 |
| Scale | 5 | 320 |

**Scale scenario:** The scale scenario uses `isograph_cpu_latent`, `isograph_gpu_latent`, `isograph_vae`, `isograph_vae_gpu`, and `wgcna_gene`.
All five scale methods are included in the supplementary scale analysis (Figure S6).

Current analysis uses the **9,440 completed runs** covering all planned CPU and GPU methods.

## Results

See generated figures in `figures/` and the summary table `_m/table1_benchmark_summary.csv`.
All figures were regenerated with R, ggplot2, and patchwork.

**Primary metrics:**

- **Module recovery** - Area-under-curve score comparing predicted module assignments to ground-truth modules; higher is better (range 0 to 1).
- **Switch gene detection rate** - Recall of true isoform-switching genes among predicted module members; higher is better (range 0 to 1).
- **Non-switching gene module rate** - Fraction of non-switching genes incorrectly assigned to modules; lower is better (range 0 to 1).
- **Runtime** - Wall-clock time per run in seconds; runtime panels use log scaling where needed for readability.
- **Resource metadata** - CPU/GPU backend, requested CPU and GPU counts, requested memory, observed maximum resident memory, torch CUDA availability, and GPU peak memory where available.

**Key findings (from `figures/fig1_benchmark_overview.pdf`):**

Across the five non-scale scenarios, IsoGraph VAE and IsoGraph VAE GPU showed the strongest module recovery while maintaining complete switch-gene detection.
The GPU VAE had similar accuracy to the CPU VAE and lower median runtime in the refreshed summaries.
The graph-regularized IsoGraph model improved module recovery relative to the latent and CPU-latent variants under several stress conditions.
WGCNA retained high switch-gene detection but showed higher non-switching gene module rates than the latent and graph IsoGraph variants.
See supplementary figures (S1-S5) for parameter-resolved analyses within each scenario and Figure S6 for scaling behavior.

In the scale scenario, IsoGraph VAE and IsoGraph VAE GPU maintained high module recovery across increasing gene counts.
GPU-accelerated latent and VAE runs are now included in the scale figure and the summary table.

## Methods

### Synthetic benchmark framework

Synthetic transcriptomic datasets were generated by simulation as described in `isograph_benchmark/benchmark/synthetic_data.py`.
Each gene was modeled with two isoforms; a proportion of genes (`switching_fraction`) were designated as switching genes whose predominant isoform varies with a latent factor correlated with age and diagnostic status.
Expression counts were drawn from a negative binomial distribution (mean parameterized by isoform proportions, dispersion = 3 unless otherwise specified).
Isoform inclusion ratios (PSI) were computed from simulated counts.

### Methods evaluated

**IsoGraph** methods use gene-level and transcript-level coexpression signals derived from PSI and count features.
Baseline, Latent, Graph, CPU Latent, and GPU Latent models used significance threshold `alpha = 0.10`.
Latent, Graph, CPU Latent, and GPU Latent models selected the number of latent components by 3-fold cross-validation where applicable.
VAE models selected latent dimension from grid {2, 4, 6, 8, 12}, used hidden dimension 128 or 256 depending on gene count, trained for up to 300 epochs with patience = 35, and used significance threshold `alpha = 0.70`.
**WGCNA** was run with signed network topology; soft-thresholding power selected by scale-free fit (R2 target >= 0.85); merge cut height = 0.25.

### Statistical analysis

All metric summaries report the mean and 95% bootstrap confidence interval (10,000 resampling iterations per group).
Groups are defined by scenario x method for the main summaries and by scenario x method x parameter combination in parameter-resolved supplemental panels.
Multiple comparison correction uses the Benjamini-Hochberg procedure [@doi:10.1111/j.2517-6161.1995.tb02031.x] where applicable.

### Compute environment

Jobs were submitted to the Pittsburgh Supercomputing Center (PSC) Bridges-2 cluster [@doi:10.1145/3437359.3465593].
CPU jobs used the RM-shared partition with 16 to 64 CPU cores per task and 2 GB memory per CPU for most runs.
WGCNA CPU jobs requested 50 CPU cores and 100 GB memory for non-scale runs.
GPU jobs used the GPU-shared partition with 1 GPU, 8 CPU cores per task, and 64 to 128 GB requested memory.

Observed CPU hardware included AMD EPYC 7742, Intel Xeon Gold 6248, and Intel Xeon Platinum 8470 nodes.
Observed accelerators included Tesla V100-SXM2-32GB and NVIDIA H100 80GB HBM3.
The IsoGraph conda environment is at `/ocean/projects/bio260021p/shared/opt/envs/isograph`.
The R figure-generation environment used by the metric script is at `/ocean/projects/bio250020p/shared/opt/env/R_env`.

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

# Step 2 - Compute bootstrap CI summary statistics
bash benchmark/02_metrics/_h/step_2_summarize.sh
# or: sbatch benchmark/02_metrics/_h/step_2_summarize.sh

# Step 3 - Generate R/ggplot2 figures and Table 1
bash benchmark/02_metrics/_h/step_3_figures.sh
# or: sbatch benchmark/02_metrics/_h/step_3_figures.sh
```

The refreshed local reproducibility logs are stored in `benchmark/02_metrics/_m/logs/`:

| Log file | Purpose |
|---|---|
| `collect-local-20260501-122624.log` | Collection of 9,440 telemetry records into `synthetic_results.parquet` |
| `summarize-local-20260501-123505.log` | Bootstrap summaries into `synthetic_metric_summary.parquet` |
| `figures-local-20260501-124549.log` | R/ggplot2 figure and table regeneration |

When run through SLURM, the scripts write equivalent logs as `collect-<jobid>.log`, `summarize-<jobid>.log`, and `figures-<jobid>.log`.

All three steps complete in under 30 minutes on a login node with at least 4 CPU cores and 8 GB RAM.

## Files

| File | Description |
|---|---|
| `_h/step_1_collect.sh` | Collects `telemetry.json` files into `benchmark/01_synthetic/_m/synthetic_results.parquet` |
| `_h/step_2_summarize.sh` | Computes bootstrap CI summaries into `_m/synthetic_metric_summary.parquet` |
| `_h/step_3_figures.sh` | Generates all R/ggplot2 figures and Table 1 |
| `_m/synthetic_metric_summary.parquet` | Bootstrap summary by scenario, method, and metric; 270 rows after refresh |
| `_m/table1_benchmark_summary.csv` | Manuscript-ready benchmark summary table; 40 rows and clean column names |
| `_m/logs/` | Reproducibility logs for collection, summary, and figure/table generation |
| `figures/fig1_benchmark_overview.pdf` | Main multi-panel benchmark overview |
| `figures/figS1_idealized_switching.pdf` | Parameter sweep: switching fraction x noise SD |
| `figures/figS2_noise_stress.pdf` | Parameter sweep: count dispersion x noise SD |
| `figures/figS3_feature_interactions.pdf` | Parameter sweep: interaction strength x interaction fraction |
| `figures/figS4_nonswitching_background.pdf` | Module and FP rate vs. background switching fraction |
| `figures/figS5_unequal_abundance.pdf` | Module recovery vs. isoform abundance imbalance |
| `figures/figS6_scale.pdf` | Module recovery and runtime vs. number of genes (scale scenario) |

PNG versions of all figures are written alongside the PDF files for quick inspection.

## References

<!-- References auto-resolved by manubot from [@doi:...] tags above -->
