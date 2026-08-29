# Residual-cost ablation: is residualization free on unconfounded data?

`isograph_vae` vs `isograph_vae_residual` on **100 fresh datasets** across 5 unconfounded scenarios, paired within dataset (both arms see byte-identical input; only the residualization flag differs).

**Scope.** In unconfounded scenarios `RIN`, `neuron_frac` and `batch` are constants that `build_design_matrix` drops, so this measures the cost of residualizing **library depth only** (`library_size`). It is not a test of residualization in general and must not be reported as one.

**Why fresh datasets.** The archived unconfounded datasets predate the covariate columns, and 1,204 of them no longer regenerate from the current generator (`benchmark/01_synthetic/_m/SAMPLE_TABLE_REFRESH.md`), so their sample tables could not be repaired. A paired within-dataset contrast does not need to be commensurable with the archive, so these were generated fresh and pinned.

Generator pin: commit `012b1a37001c15a8b4f7debc87472c2777f87bf1`, `synthetic_data.py` sha256 `1771763f63581b2f`, dirty=False.

## Pooled across scenarios

`delta` is residual minus baseline, so a **negative** delta is a cost. The CI is a paired bootstrap (2,000 resamples, seed 13) on the mean delta.

| metric | n pairs | baseline | residual | delta | 95% CI | worse/better | Wilcoxon p | q (BH) |
|---|---|---|---|---|---|---|---|---|
| `metrics_ari_planted` | 100 | 0.8735 | 0.8774 | +0.0039 | +0.0002 to +0.0079 | 8/23 | 0.0255 | 0.0764 |
| `metrics_ari_assigned` | 100 | 0.9303 | 0.9325 | +0.0022 | -0.0009 to +0.0058 | 9/15 | 0.331 | 0.398 |
| `metrics_module_recovery` | 100 | 0.8579 | 0.8613 | +0.0034 | +0.0005 to +0.0066 | 12/23 | 0.0629 | 0.126 |
| `metrics_homogeneity_planted` | 100 | 0.9824 | 0.9832 | +0.0008 | -0.0005 to +0.0023 | 8/10 | 0.444 | 0.444 |
| `metrics_completeness_planted` | 100 | 0.9101 | 0.9126 | +0.0025 | +0.0005 to +0.0048 | 9/23 | 0.0151 | 0.0764 |
| `metrics_frac_planted_assigned` | 100 | 0.9355 | 0.9388 | +0.0033 | +0.0003 to +0.0066 | 9/20 | 0.0868 | 0.13 |

BH is applied across the five pooled metrics, which are five tests on the same 100 pairs. The per-scenario table below breaks those same tests down and is descriptive; it is not separately corrected.

## Per scenario (`ari_planted`)

| scenario | n pairs | baseline | residual | delta | 95% CI | Wilcoxon p |
|---|---|---|---|---|---|---|
| feature_space_interactions | 20 | 0.5896 | 0.5984 | +0.0089 | -0.0007 to +0.0182 | 0.084 |
| idealized_switching | 20 | 0.9905 | 0.9906 | +0.0001 | +0.0000 to +0.0004 | 0.317 |
| noise_stress | 20 | 0.9905 | 0.9906 | +0.0001 | +0.0000 to +0.0004 | 0.317 |
| non_switching_background | 20 | 0.8008 | 0.8102 | +0.0093 | -0.0056 to +0.0249 | 0.214 |
| unequal_isoform_abundance | 20 | 0.9961 | 0.9971 | +0.0011 | -0.0011 to +0.0038 | 0.593 |

## Verdict

The mean ARI change is +0.0039, 95% CI +0.0002 to +0.0079, BH q=0.0764: **no cost detectable** after correcting across the five metrics. The point estimate is slightly positive, which is consistent with library depth being a real nuisance axis in every scenario, but nothing here survives multiple testing and it must not be reported as a benefit.

The supportable claim is an **equivalence bounded by the CI**: residualizing library depth changes ARI by at most 0.0079 on unconfounded data. That licenses the at-no-cost clause alongside the confound-repair result -- stated as a bound, not as a proof of exactly zero effect.
