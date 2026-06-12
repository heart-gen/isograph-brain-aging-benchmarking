# Synthetic Benchmark of IsoGraph

This document summarizes the synthetic benchmark that evaluates IsoGraph against a
WGCNA baseline for recovering **gene-level modules of coordinated isoform switching**
from bulk RNA-seq. It is written for direct insertion into the manuscript (Manubot
Markdown): each number below is read from the result tables under
`benchmark/03_metrics/_m/` and `benchmark/02_interpret/_m/`, and the language is scoped
to what those tables support — including the conditions under which IsoGraph does
**not** improve on the baseline.

The benchmark is intentionally framed as a stress test rather than a leaderboard:
the scenarios are hand-designed to probe specific failure modes (abundance-dominated
modules, technical confounds, RNA degradation, the null), so the appropriate readout
is the per-scenario effect, not an aggregate rank.

### What this benchmark tests (and what it does not)

The comparison is **not** "switch representation versus gene abundance." Every method —
IsoGraph and WGCNA alike — is fit on the *identical* per-gene feature matrix (see
[Shared input](#shared-input-all-methods-see-the-same-features)). The question is
therefore strictly about *network inference*:

> Given the same per-gene feature representation (a switch coordinate and an abundance
> coordinate per gene), does IsoGraph's network inference recover planted co-switching
> modules better than a standard WGCNA correlation network?

WGCNA is **not** run on a conventional total-gene-abundance matrix, and it is **not**
given a different or privileged representation. It receives the same gene-indexed
switch + abundance feature rows that IsoGraph receives; the methods differ only in how
they turn those features into a gene-gene network and modules.

## Methods

### Design

We generated paired synthetic datasets with known module structure, isoform-switch
ground truth, and (where relevant) abundance-channel ground truth. The grid spans
**15 scenarios** evaluated over **12,530 completed runs** with **15–30 random seeds
per scenario × method cell** (`configs/synthetic_grid.yaml`). Here a "completed run"
is one **scenario × method × dataset-seed cell** — i.e. one fitted model on one dataset
(each run has a unique `run_id`). All methods are run on the *same* `dataset_id` values
within a scenario, so every method-vs-baseline comparison is paired on identical data.

*(Earlier manuscript drafts quoted "9,440 runs across 1,204 synthetic datasets"; that
figure is stale. The current grid is 12,530 runs across 15 scenarios — the Introduction,
Methods, Results, and figure captions should all use 12,530 / 15.)*

Core scenarios use 400 genes × 160 samples; dedicated scale scenarios extend to
1,000–12,000 genes (`scale`) and to a BrainSEQ-matched 16,000 genes × 300 samples
(`scale_realistic`). The scenarios fall into five groups:

| Group | Scenarios | What it probes |
|---|---|---|
| Switch recovery | `idealized_switching`, `noise_stress`, `unequal_isoform_abundance`, `non_switching_background`, `feature_space_interactions` | Module recovery when isoform switching is the defining signal, across switching fraction, count dispersion, isoform imbalance, background contamination, and nonlinear feature-space interactions |
| Abundance vs switch | `abundance_switch_mixed`, `multi_isoform_switch` | Behavior when modules are partly or wholly abundance-driven, and when genes have >2 isoforms |
| Technical confounds | `batch_effects`, `cell_composition`, `library_depth` | Robustness to processing-batch offsets, neuron/glia composition variation, and sequencing-depth spread |
| RNA degradation | `rna_degradation`, `rna_degradation_coupled` | 3′ coverage bias that corrupts the switch signal, without vs with a degradation-robust abundance channel to fall back on |
| Null / specificity | `negative_control_noise`, `non_switching_background` | False module structure under pure noise or under non-switching background |

`non_switching_background` is listed twice on purpose. It plants true switching modules
inside a large pool of non-switching background genes, so it simultaneously measures
**sensitivity** (recovery of the planted switching modules) and **specificity** (not
sweeping background genes into modules). It contributes to both the switch-recovery and
the null/specificity readouts.

### Shared input: all methods see the same features

Every backend is fed the same three tables (`transcript_counts`, `transcript_table`,
`sample_table`) and internally builds the **same** per-gene feature matrix via
`isograph.features.channels.gene_feature_channels`. That matrix contains, for each
gene:

- a **switch coordinate** — PC1 of the gene's centered-log-ratio (CLR) isoform
  composition, present for every gene with ≥ 2 retained transcripts; and
- an **abundance coordinate** — the gene's standardized log-CPM (the summed transcript
  counts of the gene), present for every gene.

WGCNA (`wgcna_gene`) is a thin subprocess wrapper around `wgcna_runner.R`
(`isograph.models.wgcna.WgcnaNetworkModel`): it builds a standard **signed**
correlation network over **all** of these feature rows — abundance and switch alike —
then maps the resulting feature modules back to genes. So WGCNA does have access to the
abundance channel; it is not switch-blind and it is not abundance-only. What differs
between methods is the network step:

| Method (grid id) | Display name | Edges formed between features |
|---|---|---|
| `wgcna_gene` | WGCNA | Signed correlation over all rows → switch–switch, abundance–abundance, and cross edges |
| `isograph_vae` | IsoGraph VAE (single-channel) | **Switch–switch only** (abundance–abundance and switch↔abundance cross edges suppressed) |
| `isograph_vae_multiplex` | IsoGraph VAE (multiplex) | Switch–switch **and** abundance–abundance |
| `isograph_vae_residual` | IsoGraph VAE (residualized) | Switch–switch only, after regressing recorded covariates out of the features |
| `isograph_vae_reliability` | IsoGraph VAE (reliability) | Multiplex, with degradation-corrupted switch edges downweighted toward the abundance channel |

The label `wgcna_gene` denotes that WGCNA produces **gene-level** modules from this
gene-indexed feature matrix (every feature is mapped to a `gene_id`); it does **not**
mean WGCNA was run on raw gene abundance. (The grid still uses the identifier
`wgcna_gene` in `configs/synthetic_grid.yaml` and in the result tables; prose and
figures use the display name "WGCNA".)

This design has a direct consequence that recurs in the Results: the single-channel
`isograph_vae` deliberately ignores abundance–abundance structure, whereas WGCNA always
sees it. When a module's identity is carried by abundance rather than switching, WGCNA
can read it directly from the shared input while single-channel IsoGraph cannot — which
is exactly why the multiplex variant exists.

### Methods compared

The primary comparison is IsoGraph's variational backend (`isograph_vae`) against the
WGCNA baseline (`wgcna_gene`) [@doi:10.1186/1471-2105-9-559]. Additional IsoGraph
configurations are run where the scenario calls for them: the multiplex abundance+switch
model (`isograph_vae_multiplex`), a covariate-residualizing variant for the confound
scenarios (`isograph_vae_residual`), a degradation-reliability variant
(`isograph_vae_reliability`), a GPU build (`isograph_vae_gpu`), and linear baselines
(`isograph_baseline`, `isograph_latent`, `isograph_graph`, `isograph_spearman_leiden`)
on the core switch-recovery scenarios.

`isograph_vae_gpu` is **numerically matched** to the CPU `isograph_vae` (same model, same
features, same graph construction) and exists to report scale **timing**. Its accuracy
rows are therefore redundant with the CPU VAE: the accuracy tables and headline figures
report the CPU `isograph_vae`, and the GPU build is reported only in the runtime/compute
tables. (For completeness, `isograph_vae_gpu` rows are still present in the paired-test
table and the FDR family; because they duplicate the CPU VAE numerically, they should be
read as the same method, not as an independent eighth comparator.)

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
- **Does not explicitly penalise module *count*** — MRS rewards the single best-matching
  predicted module per truth module and does not subtract for producing many extra
  modules. A method that shatters the data into very many small modules can therefore
  earn an inflated apparent MRS by chance overlap, even when those modules carry no real
  structure (see `negative_control_noise` below). MRS should accordingly be read
  **alongside** the predicted module count, module size, and the non-switch gene module
  rate (NGMR), not on its own.

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
Wilcoxon signed-rank test (pairs ≥ 10; pairing on `dataset_id`). Because the design is
paired by `dataset_id`, the primary effect size is the **matched-pairs rank-biserial
correlation** $r = (T^{+} - T^{-})/(T^{+} + T^{-})$; **Cliff's δ**
[@doi:10.1037/0033-2909.114.3.494] (with conventional negligible/small/medium/large
bins) is reported as a secondary, unpaired effect size.

Multiple testing is controlled by Benjamini–Hochberg FDR over the **entire** family of
scenario × metric × method comparisons [@doi:10.1111/j.2517-6161.1995.tb02031.x]. The
family size is **468**: of the 477 scenario × metric × method cells, 9 are degenerate
(zero-variance/identical paired vectors yield no test statistic) and are excluded, so
468 comparisons carry a valid p-value and enter the correction. Every "FDR < 0.05"
statement in this document and in the figures refers to this **full-family** adjusted
p-value (`p_adj`); a within-metric FDR column (`p_adj_within_metric`) is provided only
as a sensitivity analysis.

Confidence intervals reported with each per-arm mean are bootstrap percentile intervals
(10,000 resamples, 95%; `ci_low`/`ci_high` in `synthetic_metric_summary.parquet`).
Note that these are per-arm CIs; the paired tables report the mean paired difference
$\Delta$ as a point estimate (a bootstrap CI on $\Delta$ itself is not currently emitted
and would be a useful supplemental addition).

In the **main figure** (Fig. 1), the per-scenario IsoGraph-VAE-vs-WGCNA brackets are
drawn directly from the precomputed **full-family** adjusted p-value
(`p_adj` in `synthetic_pairwise_tests.parquet`) via `ggpubr::stat_pvalue_manual`, so the
stars correspond *exactly* to the reported full-family FDR. The standard tiers are shown
(`*` FDR < 0.05, `**` < 0.01, `***` < 0.001, `****` < 0.0001) and any comparison with
FDR ≥ 0.05 is left unbracketed, so **nothing below FDR < 0.05 is ever starred**. The
0.05 ≤ FDR < 0.10 band is treated as **suggestive** and lives only in the table's
`significant_10` / `p_adj_within_metric` columns as a sensitivity analysis — it is never
starred in a figure. The dose-response and compute *supplements* mark finer-grained
per-parameter-bin / per-gene-count Wilcoxon comparisons (with local per-panel BH) that
are deliberately **not** part of the 468-comparison family; their stars are visual
guides, and the full-family table remains the reference for every significance claim in
the text. We deliberately do **not** report an
across-scenario omnibus ranking: the largest balanced block is too small to separate
methods by a Nemenyi critical difference, the scenarios are designed stress tests rather
than a random sample, and a balanced block would exclude the confound- and
degradation-specific methods. Implementation:
`isograph_benchmark/stats/hypothesis_tests.py::paired_tests` (tests + effect sizes)
and `isograph_benchmark/stats/summarize.py` (bootstrap CIs).

## Results

### IsoGraph improves module recovery over WGCNA on the same switch + abundance input when switching is the signal

When isoform switching defines the module structure, the IsoGraph VAE recovers ground
-truth modules substantially better than WGCNA from the **same** per-gene feature matrix.
Every comparison below is significant at full-family FDR < 0.05 with a near-perfect
matched-pairs rank-biserial and a medium-to-large Cliff's δ (mean module recovery,
95% bootstrap CI in brackets):

| Scenario | IsoGraph VAE | WGCNA | Δ (mean) | Rank-biserial | Cliff's δ | Magnitude |
|---|---|---|---|---|---|---|
| `idealized_switching` | 0.968 [0.961, 0.974] | 0.572 [0.557, 0.586] | +0.396 | +1.00 | 0.99 | large |
| `noise_stress` | 0.974 [0.971, 0.977] | 0.629 [0.619, 0.638] | +0.346 | +1.00 | 1.00 | large |
| `non_switching_background` | 0.841 [0.778, 0.898] | 0.535 [0.499, 0.570] | +0.306 | +0.85 | 0.76 | large |
| `feature_space_interactions` | 0.769 [0.739, 0.799] | 0.476 [0.456, 0.496] | +0.293 | +1.00 | 0.61 | large |
| `unequal_isoform_abundance` | 0.780 [0.707, 0.849] | 0.553 [0.524, 0.580] | +0.227 | +0.76 | 0.59 | large |
| `scale` (1k–12k genes) | 0.820 [0.787, 0.853] | 0.640 [0.634, 0.645] | +0.180 | +0.94 | 0.45 | medium |

Switch-gene detection rate is at or near ceiling (≥ 0.99) for the VAE across these
scenarios, so the recovery gain reflects better *module assignment* of switch genes,
not merely detecting more of them.

### The advantage narrows to parity at BrainSEQ scale

At the BrainSEQ-matched dimension (16,000 genes × 300 samples, `scale_realistic`),
IsoGraph VAE and WGCNA are statistically indistinguishable on module recovery
(0.625 [0.622, 0.628] vs 0.625 [0.622, 0.628]; full-family FDR not significant). The
benchmark therefore supports IsoGraph as the stronger method at the 400–12,000-gene
scale and as *comparable* — not superior — at full transcriptome scale on this metric.

### IsoGraph is more specific under the null

`negative_control_noise` sets `switching_fraction = 0`, so the data-generating model
expresses **no** module signal. (The truth table retains a single degenerate planted
module gene purely so MRS is mathematically defined — `n_module_genes = max(1, 0) = 1` —
but the simulated counts carry no coordinated structure.) Under this null, IsoGraph VAE
returns near-zero module recovery (0.054 [0.017, 0.100]) and ~9 predicted modules,
whereas WGCNA returns a spuriously high recovery score (0.405 [0.359, 0.446]) by
fragmenting the noise into ~282 tiny modules. This is the inflation that MRS does not
penalise directly (see the metric note above): WGCNA's many small modules raise the
chance that one of them happens to overlap the lone truth gene. Read alongside the
predicted module count, IsoGraph's behaviour is the desired one — it does not
manufacture module structure where none exists.

### WGCNA remains stronger when modules are abundance-dominated, despite matched input

When module identity is carried mainly by gene-level abundance rather than isoform
switching (`abundance_switch_mixed`), WGCNA is the most accurate method
(0.675 [0.653, 0.697]), the single-channel VAE is weakest
(0.331 [0.311, 0.352]; vs WGCNA rank-biserial −0.96, Cliff's δ −0.72, large), and the
multiplex model recovers much of the gap (0.492 [0.472, 0.513]). This is **not** an
input advantage: WGCNA reads the abundance signal straight from the shared feature
matrix, because its correlation network always includes the abundance–abundance edges.
The single-channel `isograph_vae` *suppresses* abundance–abundance edges by design and
so cannot see an abundance-defined module; `isograph_vae_multiplex` re-enables those
edges and closes most of the gap. We report this as a genuine limitation of the
switch-only configuration, and as the motivation for the multiplex model when abundance
is the dominant axis.

### Covariate residualization restores robustness under technical confounds

Technical confounds degrade the single-channel VAE, but the residualizing variant
(`isograph_vae_residual`) recovers or exceeds WGCNA (mean module recovery, 95% CI):

| Scenario | IsoGraph VAE | IsoGraph VAE (residualized) | WGCNA |
|---|---|---|---|
| `cell_composition` | 0.372 [0.324, 0.422] | **0.626 [0.614, 0.637]** | 0.505 [0.478, 0.530] |
| `library_depth` | 0.362 [0.314, 0.408] | **0.616 [0.598, 0.632]** | 0.457 [0.425, 0.489] |
| `batch_effects` | 0.520 [0.478, 0.560] | **0.622 [0.609, 0.636]** | 0.565 [0.536, 0.593] |

**What each method residualizes.** This is important for interpreting the contrast.
The simulator records the active confound as a sample covariate (`neuron_frac`,
`library_size`, or `batch`). `isograph_vae_residual` regresses these recorded confound
covariates out of the features before fitting. `isograph_vae` and **WGCNA** use the
default integrity-covariate list, of which only `RIN` exists in the simulated
`sample_table` (the others are silently skipped); they therefore do **not** residualize
the scenario-specific confound — they receive input that is effectively non-residualized
with respect to it.

Consequently the clean, like-for-like WITH/WITHOUT-residualization ablation is
`isograph_vae` vs `isograph_vae_residual` (same architecture, differing only in the
covariate set). The `isograph_vae_residual`-vs-WGCNA gap reflects residualization **plus**
the network-inference difference and should be read as such; it is significant at
full-family FDR < 0.05 with rank-biserial +0.93/+0.94/+0.53 (small-to-large Cliff's δ)
across the three scenarios. The practical recommendation is to residualize recorded
nuisance covariates before fitting either method.

### RNA degradation defines a domain-of-validity limit

3′ coverage bias corrupts the isoform-ratio (switch) signal. Which methods have an
abundance fallback is the crux here, so it is worth stating explicitly: `isograph_vae`
and `isograph_vae_residual` form **switch–switch edges only** (no abundance fallback),
`isograph_vae_multiplex` and `isograph_vae_reliability` enable abundance–abundance edges
(a fallback path), and WGCNA always includes both channels.

With no abundance channel in the data-generating model to fall back on
(`rna_degradation`, `abundance_fraction = 0`), both switch-only VAE variants fall below
WGCNA (0.466 [0.424, 0.508] / 0.472 [0.428, 0.514] vs 0.613 [0.588, 0.635]) — an honest
limit of a switch-only method. When a degradation-robust abundance channel **is** present
(`rna_degradation_coupled`) and reliability weighting downweights corrupted switch edges
toward it, the multiplex and reliability variants recover most of the lost signal
(0.679 [0.672, 0.686] / 0.680 [0.673, 0.686]) relative to the switch-only VAE
(0.196 [0.190, 0.202]), though still below WGCNA (0.755 [0.703, 0.802]) — which reads the
degradation-robust abundance edges directly. Degradation robustness thus requires the
abundance channel and does not, on its own, beat WGCNA.

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
are numerically matched, so the GPU build is a drop-in accelerator, not a different model
(see `tableS_scale_compute_summary.csv` for the full genes × method timing/memory grid).

## Summary

Across 12,530 paired synthetic runs — every method fit on the **same** per-gene
switch + abundance feature matrix — IsoGraph's VAE backend recovers switch-defined
modules more accurately than WGCNA whenever isoform switching is the dominant signal
(large rank-biserial effects, full-family FDR < 0.05), is more specific under the null,
and reaches parity — not superiority — at full transcriptome scale. WGCNA remains the
stronger method when module structure is abundance-dominated (it reads the shared
abundance channel directly, which single-channel IsoGraph suppresses) or when 3′
degradation corrupts the switch channel without an abundance fallback. Covariate
residualization is required for robustness to technical confounds, and IsoGraph's
interpretation layer reliably identifies the switch-driver transcript while only modestly
recovering switch magnitude. IsoGraph is best positioned as a complementary, switch-aware
layer rather than a universal replacement for gene-level co-expression analysis.

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

Result tables: `benchmark/03_metrics/_m/synthetic_metric_summary.parquet` (per-arm means
+ bootstrap CIs), `synthetic_metric_long.parquet` (per-run tidy table),
`synthetic_pairwise_tests.parquet` (paired tests + rank-biserial + Cliff's δ + full-family
FDR), and `benchmark/02_interpret/_m/synthetic_interpret_summary.parquet` (interpretation
accuracy). Figures and Table 1 are written under `benchmark/03_metrics/figures/`.

*Citation note:* `[@...]` entries are Manubot citation keys; the IsoGraph software
itself should be cited from its repository metadata
(`https://github.com/heart-gen/IsoGraph`, `CITATION.cff`) until a preprint DOI exists.
