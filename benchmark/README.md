# Synthetic Benchmark of IsoGraph

This document summarizes the synthetic benchmark that evaluates IsoGraph against a
gene-level WGCNA baseline for recovering co-regulated transcript programs from bulk
RNA-seq. It is written for direct insertion into the manuscript (Manubot Markdown):
each number below is read from the result tables under `benchmark/03_metrics/_m/`
and `benchmark/02_interpret/_m/`, and the language is scoped to what those tables
support — including the conditions under which IsoGraph does **not** improve on the
baseline.

The benchmark is intentionally framed as a stress test rather than a leaderboard:
the scenarios are hand-designed to probe specific failure modes (abundance-dominated
modules, technical confounds, RNA degradation, the null), so the appropriate readout
is the per-scenario effect, not an aggregate rank.

## Methods

### Design

We generated paired synthetic datasets with known module structure, isoform-switch
ground truth, and (where relevant) abundance-channel ground truth. The grid spans
**15 scenarios** evaluated over **12,530 completed runs** with **15–30 random seeds
per scenario × method cell** (`configs/synthetic_grid.yaml`). All methods are run on
the *same* `dataset_id` values within a scenario, so every method-vs-baseline
comparison is paired on identical data.

Core scenarios use 400 genes × 160 samples; dedicated scale scenarios extend to
1,000–12,000 genes (`scale`) and to a BrainSEQ-matched 16,000 genes × 300 samples
(`scale_realistic`). The scenarios fall into five groups:

| Group | Scenarios | What it probes |
|---|---|---|
| Switch recovery | `idealized_switching`, `noise_stress`, `unequal_isoform_abundance`, `non_switching_background`, `feature_space_interactions` | Module recovery when isoform switching is the defining signal, across switching fraction, count dispersion, isoform imbalance, background contamination, and nonlinear feature-space interactions |
| Abundance vs switch | `abundance_switch_mixed`, `multi_isoform_switch` | Behavior when modules are partly or wholly abundance-driven, and when genes have >2 isoforms |
| Technical confounds | `batch_effects`, `cell_composition`, `library_depth` | Robustness to processing-batch offsets, neuron/glia composition variation, and sequencing-depth spread |
| RNA degradation | `rna_degradation`, `rna_degradation_coupled` | 3′ coverage bias that corrupts the switch signal, without vs with a degradation-robust abundance channel to fall back on |
| Null / specificity | `negative_control_noise`, `non_switching_background` | False module structure under pure noise or non-switching background |

### Methods compared

The primary comparison is IsoGraph's variational backend (`isograph_vae`) against a
gene-level WGCNA baseline (`wgcna_gene`) [@doi:10.1186/1471-2105-9-559]. Additional
IsoGraph configurations are run where the scenario calls for them: the multiplex
abundance+switch model (`isograph_vae_multiplex`), a covariate-residualizing variant
for the confound scenarios (`isograph_vae_residual`), a degradation-reliability
variant (`isograph_vae_reliability`), a GPU build (`isograph_vae_gpu`, numerically
matched to the CPU VAE and used for scale timing), and linear baselines
(`isograph_baseline`, `isograph_latent`, `isograph_graph`, `isograph_spearman_leiden`)
on the core switch-recovery scenarios.

### Metrics

Per run we record the three recovery metrics defined below, the number of predicted
modules, edge count, and wall-clock runtime; the multiplex scenarios additionally
record per-gene role recovery (`switch_only`, `abundance_only`, `coupled`,
`discordant`). Interpretation accuracy is scored separately in
`benchmark/02_interpret/` using top-1 switch-driver transcript accuracy and
switch-magnitude Spearman correlation.

#### Module recovery score

Let $\mathcal{T} = \{T_1, \dots, T_K\}$ be the set of ground-truth modules and
$\mathcal{P} = \{P_1, \dots, P_M\}$ the set of predicted modules. The module recovery
score is the mean best-Jaccard similarity:

$$\text{MRS}(\mathcal{P}, \mathcal{T}) = \frac{1}{K} \sum_{i=1}^{K} \max_{j \in [M]} \frac{|T_i \cap P_j|}{|T_i \cup P_j|}$$

Each ground-truth module $T_i$ is matched to the predicted module $P_j$ that maximises
the Jaccard index. The score ranges from 0 (no overlap between any truth and predicted
module) to 1 (every truth module is exactly recovered by some predicted module).

Properties:

- **Asymmetric in direction** — measures how well predicted modules cover the truth,
  not how well truth covers predictions. A method that predicts one giant module
  covering all genes scores poorly because per-module Jaccard is diluted by the large
  union.
- **Penalises fragmentation** — if a truth module $T_i$ is split across two predicted
  modules $P_a$ and $P_b$, the best Jaccard is below 1 even when
  $T_i \subseteq P_a \cup P_b$.
- **Independent of module count** — the number of predicted modules affects the score
  only through partition quality, not directly.

Implementation: `isograph.evaluation.metrics.module_recovery_score(predicted, truth)`,
where both arguments are DataFrames with columns `gene_id` and `module_id`. The
best-match Jaccard approach for cluster-recovery evaluation follows Lancichinetti &
Fortunato (2009) *Physical Review E* 80:056117 (used by igraph `compare_communities`
under the name "best-match F1 / Jaccard").

#### Switch-gene detection rate

Fraction of ground-truth switching genes assigned to any predicted module:

$$\text{SGDR} = \frac{|\hat{G} \cap G_\text{switch}|}{|G_\text{switch}|}$$

where $\hat{G}$ is the set of all genes assigned to any predicted module and
$G_\text{switch}$ is the set of genes whose PSI varies with the simulated trait.

#### Non-switch gene module rate (false-positive rate)

Fraction of non-switching background genes that are assigned to a predicted module:

$$\text{NGMR} = \frac{|\hat{G} \cap G_\text{background}|}{|G_\text{background}|}$$

Lower is better. A method that assigns every gene to a module scores 1.0 regardless of
specificity. Implementation of SGDR and NGMR:
`isograph_benchmark/benchmark/run_one.py::compute_metrics`.

### Statistics

For each scenario × metric we compare each method against `wgcna_gene` with a paired
Wilcoxon signed-rank test (pairs ≥ 10; pairing on `dataset_id`). Effect sizes are
reported alongside every p-value: the matched-pairs rank-biserial correlation
((T⁺ − T⁻)/(T⁺ + T⁻)) and Cliff's δ [@doi:10.1037/0033-2909.114.3.494] with
conventional negligible/small/medium/large magnitude bins. Multiple testing is
controlled by Benjamini–Hochberg FDR over the **entire** family of
scenario × metric × method comparisons (family size = 468)
[@doi:10.1111/j.2517-6161.1995.tb02031.x]; a within-metric FDR column is also
provided as a sensitivity analysis. Confidence intervals are bootstrap percentile
intervals (10,000 resamples, 95%). We deliberately do **not** report an
across-scenario omnibus ranking: the largest balanced block is too small to separate
methods by a Nemenyi critical difference, the scenarios are designed stress tests
rather than a random sample, and a balanced block would exclude the confound- and
degradation-specific methods. In figures and tables, significance is marked
`**` for FDR < 0.05 and `*` for FDR < 0.10. Implementation:
`isograph_benchmark/stats/hypothesis_tests.py::paired_tests` (tests + effect sizes)
and `isograph_benchmark/stats/summarize.py` (bootstrap CIs).

## Results

### IsoGraph recovers switch-defined modules more accurately when switching is the signal

When isoform switching defines the module structure, the IsoGraph VAE recovers ground
-truth modules substantially better than gene-level WGCNA. Every comparison below is
significant at FDR < 0.05 with a medium or large Cliff's δ (mean module recovery,
IsoGraph VAE vs WGCNA):

| Scenario | IsoGraph VAE | WGCNA | Δ (mean) | Cliff's δ | Magnitude |
|---|---|---|---|---|---|
| `idealized_switching` | 0.968 | 0.572 | +0.396 | 0.99 | large |
| `noise_stress` | 0.974 | 0.629 | +0.346 | 1.00 | large |
| `non_switching_background` | 0.841 | 0.535 | +0.306 | 0.76 | large |
| `feature_space_interactions` | 0.769 | 0.476 | +0.293 | 0.61 | large |
| `unequal_isoform_abundance` | 0.780 | 0.553 | +0.227 | 0.59 | large |
| `scale` (1k–12k genes) | 0.820 | 0.640 | +0.180 | 0.45 | medium |

Switch-gene detection rate is at or near ceiling (≥ 0.99) for the VAE across these
scenarios, so the recovery gain reflects better *module assignment* of switch genes,
not merely detecting more of them.

### The advantage narrows to parity at BrainSEQ scale

At the BrainSEQ-matched dimension (16,000 genes × 300 samples, `scale_realistic`),
IsoGraph VAE and WGCNA are statistically indistinguishable on module recovery
(0.625 vs 0.625; not significant). The benchmark therefore supports IsoGraph as the
stronger method at the 400–12,000-gene scale and as *comparable* — not superior — at
full transcriptome scale on this metric.

### IsoGraph is more specific under the null

Under pure noise with no true module signal (`negative_control_noise`), IsoGraph VAE
returns near-zero module recovery (0.054) and ~9 predicted modules, whereas WGCNA
returns a spuriously high recovery score (0.405) by fragmenting the data into ~282
modules. IsoGraph's lower null recovery is the desired behavior: it does not
manufacture module structure where none exists.

### WGCNA remains stronger when modules are abundance-dominated

When module identity is carried mainly by gene-level abundance rather than isoform
switching (`abundance_switch_mixed`), gene-level WGCNA is the most accurate method
(0.675), the single-channel VAE is weakest (0.331), and the multiplex model recovers
much of the gap (0.492) by reading the abundance channel. We report this as a genuine
limitation: IsoGraph's strength is switch-defined structure, and a gene-level method
is the right tool when abundance is the dominant axis.

### Covariate residualization restores robustness under technical confounds

Technical confounds degrade the single-channel VAE, but the residualizing variant
(`isograph_vae_residual`) recovers or exceeds WGCNA (mean module recovery):

| Scenario | IsoGraph VAE | IsoGraph VAE (residualized) | WGCNA |
|---|---|---|---|
| `cell_composition` | 0.372 | **0.626** | 0.505 |
| `library_depth` | 0.362 | **0.616** | 0.457 |
| `batch_effects` | 0.520 | **0.622** | 0.565 |

The practical recommendation is therefore to residualize recorded nuisance covariates
before fitting; the gain over WGCNA under confounds is significant (FDR < 0.05) with
small-to-large Cliff's δ.

### RNA degradation defines a domain-of-validity limit

3′ coverage bias corrupts the isoform-ratio (switch) signal. With no abundance channel
to fall back on (`rna_degradation`), both VAE variants fall below WGCNA
(0.466 / 0.472 vs 0.613) — an honest limit of a switch-only method. When a
degradation-robust abundance channel is available and reliability weighting downweights
corrupted switch edges (`rna_degradation_coupled`), the multiplex and reliability
variants recover most of the lost signal (0.679 / 0.680) relative to the single-channel
VAE (0.196), though still below the gene-level WGCNA baseline (0.755). Degradation
robustness thus requires the abundance channel and does not, on its own, beat a
gene-level method.

### Module interpretation correctly identifies the switch-driver transcript

On genes with four isoforms (`multi_isoform_switch`, chance = 0.25), IsoGraph's built-in
interpretation identifies the correct switch-driver transcript with top-1 accuracy of
0.999 (`isograph_vae`) and 1.000 (`isograph_vae_multiplex`). Recovery of the *signed
magnitude* of switching is weaker — switch-magnitude Spearman correlations are modest
and scenario-dependent (e.g. 0.23 under `cell_composition`, 0.05 under
`multi_isoform_switch`) — so the interpretation layer is reliable for *which* transcript
drives a module but should be read cautiously for *how much*.

### Runtime is competitive and scales with GPU acceleration

On the core 400-gene grid, IsoGraph VAE and WGCNA have comparable median runtimes
(~4–5 s). At BrainSEQ scale (`scale_realistic`), the GPU VAE is the fastest method
(24.6 s) versus the CPU VAE (226.3 s) and WGCNA (383.5 s); the GPU and CPU VAE backends
are numerically matched.

## Summary

Across 12,530 paired synthetic runs, IsoGraph's VAE backend recovers switch-defined
co-expression modules more accurately than gene-level WGCNA whenever isoform switching
is the dominant signal (large effect sizes, FDR < 0.05), is more specific under the
null, and reaches parity — not superiority — at full transcriptome scale. WGCNA remains
the stronger method when module structure is abundance-dominated or when 3′ degradation
corrupts the switch channel without an abundance fallback. Covariate residualization is
required for robustness to technical confounds, and IsoGraph's interpretation layer
reliably identifies the switch-driver transcript while only modestly recovering switch
magnitude. IsoGraph is best positioned as a complementary, switch-aware layer rather
than a universal replacement for gene-level co-expression analysis.

## Reproducing these results

```bash
# 1. Generate the synthetic grid and run all method × scenario × seed cells
sbatch benchmark/00_design/_h/*.sh
sbatch benchmark/01_synthetic/_h/*.sh

# 2. Score module interpretation (runs before metrics; metrics figS11 reads its summary)
sbatch benchmark/02_interpret/_h/step_1.sh
sbatch benchmark/02_interpret/_h/step_2.sh

# 3. Collect, summarize (bootstrap CIs + paired tests + effect sizes), and plot
sbatch benchmark/03_metrics/_h/step_1_collect.sh
python -m isograph_benchmark.stats.summarize        # step_2_summarize.sh
sbatch benchmark/03_metrics/_h/step_3_figures.sh
```

Result tables: `benchmark/03_metrics/_m/synthetic_metric_summary.parquet` (bootstrap
CIs), `synthetic_metric_long.parquet` (per-run tidy table),
`synthetic_pairwise_tests.parquet` (paired tests + effect sizes + FDR), and
`benchmark/02_interpret/_m/synthetic_interpret_summary.parquet` (interpretation
accuracy). Figures and Table 1 are written under `benchmark/03_metrics/figures/`.

*Citation note:* `[@...]` entries are Manubot citation keys; the IsoGraph software
itself should be cited from its repository metadata
(`https://github.com/heart-gen/IsoGraph`, `CITATION.cff`) until a preprint DOI exists.
