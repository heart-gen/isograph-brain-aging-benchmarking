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

This directory contains the downstream metric analysis for the IsoGraph synthetic benchmark described in `01_synthetic_benchmark/01_synthetic/`.
The benchmark evaluates IsoGraph, a Spearman-Leiden baseline, and WGCNA [@doi:10.1186/1471-2105-9-559] across nine parameterized synthetic scenarios designed to stress-test network recovery under conditions relevant to brain-aging transcriptomics.
Figures and tables here are intended for inclusion in the IsoGraph manuscript targeting *Nature Methods*.
The collected result table includes CPU, SLURM, memory, runtime, software, and hardware metadata for every completed run.
Scale, scale_realistic, and abundance_switch_mixed scenarios are complete; the `isograph_vae_multiplex` sweep across the remaining scenarios is still in progress.

## Synthetic Data Generation

Synthetic datasets were generated using a negative binomial count model with ground-truth gene modules and isoform-switching signals.
Each non-scale dataset contains 400 genes (800 transcripts) across 160 samples with simulated covariates (age, sex, diagnosis).
Switching genes are assigned to latent modules correlated with age and disease status; the proportion of switching genes (`switching_fraction`) and expression noise (`noise_sd`) are varied across scenarios.

The `scale` scenario uses 1,000–12,000 genes across 240 samples; `scale_realistic` uses 16,000 genes across 300 samples (matching BrainSEQ dimensions).
Datasets were generated with a deterministic seed sequence (base seed 13) to ensure reproducibility.
Each dataset is identified by a SHA1 hash of its scenario parameters and seed, ensuring stable run IDs across grid changes.

## Benchmark Design

Nine scenarios are evaluated across the applicable method sets.
All non-scale scenarios are run for the full set of eight methods: seven CPU methods
(including `isograph_vae_multiplex`) plus the GPU-only `isograph_vae_gpu`, which is the same
VAE model run on GPU and reported only as a supplementary compute comparison.
The scale scenarios are run for IsoGraph VAE (CPU and GPU) and WGCNA only (other methods are
memory-prohibitive at 16k genes). The `abundance_switch_mixed` scenario uses the three-method
multiplex set.

| Scenario | Key parameters | Seeds | Methods |
|---|---|---:|---:|
| **Idealized switching** | `switching_fraction` in {0.10, 0.25, 0.50, 0.75}; `noise_sd` in {0.10, 0.25, 0.40} | 30 | 8 |
| **Noise stress** | `count_dispersion` in {3, 7, 15, 30}; `noise_sd` in {0.05, 0.10, 0.25, 0.50} | 20 | 8 |
| **Feature space interactions** | `interaction_strength` in {0.0, 0.5, 1.0, 2.0, 3.0}; `interaction_fraction` in {0.25, 0.50, 0.75} | 20 | 8 |
| **Non-switching background** | `switching_fraction` in {0.05, 0.10, 0.25, 0.50} | 20 | 8 |
| **Unequal isoform abundance** | `abundance_imbalance` in {1, 4, 10, 25} | 20 | 8 |
| **Abundance switch mixed** | `abundance_fraction` in {0.0, 0.3, 0.5, 0.7, 1.0}; `n_genes` in {200, 500} | 20 | 3 (multiplex set) |
| **Negative control noise** | `switching_fraction` = 0.0; `noise_sd` in {0.5, 1.0} | 20 | 8 |
| **Scale** | `n_genes` in {1000, 3000, 6000, 12000}; `switching_fraction` in {0.15, 0.25} | 15 | 3 |
| **Scale realistic** | `n_genes` = 16000; `n_samples` = 300; `switching_fraction` in {0.15, 0.25} | 15 | 3 |

**Methods evaluated:**

| Method | Description | Compute | Scenarios |
|---|---|---|---|
| `isograph_baseline` | IsoGraph baseline network model (`alpha = 0.10`) | CPU | all except scale, multiplex |
| `isograph_latent` | IsoGraph latent-space model (3-fold CV for components, `alpha = 0.10`) | CPU | all except scale, multiplex |
| `isograph_graph` | IsoGraph graph-regularized model (3-fold CV, `alpha = 0.10`) | CPU | all except scale, multiplex |
| `isograph_vae` | Variational autoencoder network model | CPU | all |
| `isograph_vae_gpu` | Same VAE model and config as `isograph_vae`, run on GPU; supplementary compute comparison only | GPU | all except multiplex |
| `isograph_vae_residual` | `isograph_vae` that regresses recorded nuisance covariates (RIN, neuron_frac, batch, library_size) out of the features before graph construction | CPU | confound scenarios |
| `isograph_vae_multiplex` | IsoGraph VAE with multiplex (abundance + switch) features; `leiden_resolution = 2.0` to prevent the dense abundance channel fusing modules into one giant component | CPU | multiplex, rna_degradation_coupled |
| `isograph_vae_reliability` | `isograph_vae_multiplex` plus degradation-aware switch reliability: per-gene switch edges are downweighted by alignment with an observed 3' coverage covariate (`median_tin`), so degradation-corrupted genes fall back to the abundance channel (approach #3) | CPU | rna_degradation_coupled |
| `isograph_spearman_leiden` | Spearman r on the shared abundance+switch feature matrix + Leiden clustering (`min_r = 0.30`); same input and feature→gene mapping as WGCNA | CPU | all except scale, multiplex |
| `wgcna_gene` | WGCNA gene-level coexpression [@doi:10.1186/1471-2105-9-559] | CPU | all |

Total planned runs: **12,450**.

### Degradation robustness (rna_degradation vs rna_degradation_coupled)

`rna_degradation` sweeps a one-sided per-gene 3' coverage bias on pure-switch
modules (`abundance_fraction = 0`): with no abundance channel to fall back to it
shows IsoGraph's domain-of-validity limit. `rna_degradation_coupled` adds a
degradation-robust abundance channel (`abundance_fraction = 0.4`,
`dual_signal_fraction = 0.5`) and the generator records `median_tin`, an observed
sample-level 3' coverage metric that is a sharper proxy of the artifact than RIN.
Spot-test (mean module recovery, fast config; full-config sweep pending):

| bias | vae | vae_residual | multiplex | reliability | wgcna |
|---|---|---|---|---|---|
| 0.0 | 0.220 | 0.220 | 0.676 | 0.676 | 0.876 |
| 1.0 | 0.201 | 0.201 | 0.685 | 0.687 | 0.885 |
| 2.0 | 0.166 | 0.166 | 0.678 | 0.682 | 0.780 |

Takeaways: (1) regressing RIN out does nothing (`vae_residual ≡ vae`) — one-sided
per-gene degradation is ~orthogonal to the RIN axis; (2) the abundance channel is
the dominant lever (0.22→0.68), but only after `leiden_resolution` fixes a
giant-module collapse the dense abundance edges otherwise cause under connected
components; (3) reliability via `median_tin` adds a small, consistent edge over
multiplex that grows with degradation and never hurts (identical at bias 0);
(4) IsoGraph holds ~0.68 while WGCNA degrades (0.885→0.780), so the gap is smallest
under heavy degradation. The residual ~0.20 gap at bias 0 is background/grey
rejection, not degradation.

## Run Status

The remaining work is the `isograph_vae_multiplex` sweep across the non-multiplex scenarios,
the new `isograph_vae_gpu` runs for the scenarios added since its telemetry was captured
(negative control and the expanded scale grids), plus a small tail of
`isograph_spearman_leiden` and `isograph_vae` runs.

| Scenario | Methods | Planned | Completed |
|---|---|---:|---:|
| Idealized switching | 8 | 2,880 | 2,520 |
| Noise stress | 8 | 2,560 | 2,240 |
| Feature space interactions | 8 | 2,400 | 2,165 |
| Non-switching background | 8 | 640 | 555 |
| Unequal isoform abundance | 8 | 640 | 480 |
| Abundance switch mixed | 3 | 1,200 | 1,199 |
| Negative control noise | 8 | 320 | 200 |
| Scale | 3 (VAE CPU/GPU + WGCNA) | 360 | 304 |
| Scale realistic | 3 (VAE CPU/GPU + WGCNA) | 90 | 60 |
| **Total** | | **11,090** | **9,723** |

**Pending runs (1,367):**
- `isograph_vae_multiplex` across the non-multiplex scenarios (1,116 runs) — `vae` resource class
- `isograph_vae_gpu` new scenarios (40 `gpu` + 86 `gpu_scale` = 126 runs) — GPU resource classes
- `isograph_spearman_leiden` tail (85 runs) — `cpu_short` resource class
- `isograph_vae` (40, negative-control) — `vae` resource class

Submit the remaining work with `run_batch_vae.sh` (1,156 missing), `run_batch_cpu_short.sh`
(85 missing), `run_batch_gpu.sh` (40 missing), and `run_batch_gpu_scale.sh` (86 missing); all
skip runs that already have a `done.json`. Current figures and tables are regenerated by the
collect → summarize → figures pipeline over whatever runs are complete.

## Results

See generated figures in `figures/`, the supplementary benchmark summary table `_m/tableS_benchmark_summary.csv` (six core accuracy scenarios × six main methods), and the scale compute table `_m/tableS_scale_compute_summary.csv`.
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
Baseline, Latent, and Graph models used significance threshold `alpha = 0.10`.
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
The R figure-generation environment is at `/ocean/projects/bio250020p/shared/opt/env/R_env`.

## Reproduction

All outputs in this directory can be regenerated from the raw run outputs in `01_synthetic_benchmark/01_synthetic/_o/runs/`.
Scripts in `_h/` can be run directly on a login node or submitted via `sbatch` for cluster execution.

SLURM opens the output log path before the script body runs, so create the metrics log directory before submitting jobs:

```bash
mkdir -p 01_synthetic_benchmark/03_metrics/_m/logs
```

```bash
# Step 1 - Collect all telemetry into a single parquet
bash 01_synthetic_benchmark/03_metrics/_h/step_1_collect.sh
# or: sbatch 01_synthetic_benchmark/03_metrics/_h/step_1_collect.sh

# Step 2 - Compute bootstrap CI summaries and paired Wilcoxon tests (BH FDR)
bash 01_synthetic_benchmark/03_metrics/_h/step_2_summarize.sh
# or: sbatch 01_synthetic_benchmark/03_metrics/_h/step_2_summarize.sh

# Step 3 - Generate figures and summary tables (R/ggplot2)
#   (ISOGRAPH_TABLES_ONLY=1 regenerates only the CSV summary tables, no figures)
bash 01_synthetic_benchmark/03_metrics/_h/step_3_figures.sh
# or: sbatch 01_synthetic_benchmark/03_metrics/_h/step_3_figures.sh
```

Step 3 runs a single publication-quality figure script,
`isograph_benchmark/figures/synthetic_benchmark.R`
(ggpubr/ggplot2/patchwork, requires `R_env`). It is data-driven: methods and
scenarios are drawn from whatever completed runs are present, so newly added
methods (e.g. `isograph_spearman_leiden`, `isograph_vae_multiplex`) and
scenarios (`negative_control_noise`, `abundance_switch_mixed`) appear
automatically once the benchmark is re-run; absent ones are skipped.

All three steps complete in under 30 minutes on a login node with at least 4 CPU cores and 8 GB RAM.
When run through SLURM, the scripts write logs as `collect-<jobid>.log`, `summarize-<jobid>.log`, and `figures-<jobid>.log` in `_m/logs/`.

## Files

| File | Description |
|---|---|
| `_h/step_1_collect.sh` | Collects `telemetry.json` files into `01_synthetic_benchmark/01_synthetic/_m/synthetic_results.parquet` |
| `_h/step_2_summarize.sh` | Computes bootstrap CI summaries into `_m/synthetic_metric_summary.parquet` |
| `_h/step_3_figures.sh` | Generates all R/ggplot2 figures and summary tables |
| `_m/synthetic_results.parquet` | Raw collected telemetry; 9,617 rows including historical and current runs |
| `_m/synthetic_metric_summary.parquet` | Bootstrap summary by scenario, method, and metric; 282 rows |
| `_m/synthetic_metric_long.parquet` | Per-run long-format metrics; 57,702 rows |
| `_m/synthetic_pairwise_tests.parquet` | Paired Wilcoxon signed-rank tests vs. WGCNA with BH FDR; 240 tests |
| `_m/tableS_benchmark_summary.csv` | Supplementary CPU benchmark summary table (six core accuracy scenarios × six main methods) |
| `_m/tableS_scale_compute_summary.csv` | Scale compute table for IsoGraph VAE and WGCNA across gene-count grid |
| `_m/logs/` | Reproducibility logs for collection, summary, and figure/table generation |
| `figures/fig1_benchmark_overview.pdf` | Main CPU accuracy and specificity benchmark overview |
| `figures/figS1_idealized_switching.pdf` | Response dot plot: switching fraction by noise SD |
| `figures/figS2_noise_stress.pdf` | Response dot plot: count dispersion by noise SD |
| `figures/figS3_feature_interactions.pdf` | Response dot plot: interaction strength by interaction fraction |
| `figures/figS4_nonswitching_background.pdf` | Module recovery and non-switching specificity vs. switching fraction |
| `figures/figS5_unequal_abundance.pdf` | Module recovery and switch detection vs. isoform abundance imbalance |
| `figures/figS6_scale_compute.pdf` | Runtime and peak host RAM vs. number of genes (scale + scale_realistic) |
| `figures/figS7_abundance_roles.pdf` | Abundance-shift gene detection and isoform-role composition (abundance_switch_mixed; written only when those runs exist) |

PNG versions of all figures are written alongside the PDF files for quick inspection.

## References

<!-- References auto-resolved by manubot from [@doi:...] tags above -->
