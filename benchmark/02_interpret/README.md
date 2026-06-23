# Synthetic Module Interpretation

This stage evaluates IsoGraph module interpretation on the synthetic benchmark
using ground-truth switching genes and module assignments.

Run (sharded array, then collect):

```bash
JOB=$(sbatch --parsable benchmark/02_interpret/_h/step_1.sh)
sbatch --dependency=afterok:$JOB benchmark/02_interpret/_h/step_2.sh
```

`step_1.sh` is a 16-way array: each task evaluates `runs[shard::16]` and writes a
per-shard partial under `_m/_shards/`. `step_2.sh` merges the partials and runs the
bootstrap summary once. The expensive per-run `explain-module` output is
checkpointed on disk, so reruns reuse it.

Main outputs:

- `benchmark/02_interpret/_m/synthetic_interpret_results.parquet`
- `benchmark/02_interpret/_m/synthetic_interpret_module_metrics.parquet`
- `benchmark/02_interpret/_m/synthetic_interpret_summary.parquet`

Large per-run `explain-module` outputs are written under
`benchmark/02_interpret/_o/runs/` and are intentionally ignored by git.

Synthetic transcript IDs are simulated, so this stage does not use a human GTF.
GTF-based structural annotation is enabled by default for the real-data
interpretation runners.

## Metrics

Gene-level driver / role attribution (all scenarios, vs `truth_switch`/
`truth_modules`/`truth_abundance`): `gene_driver_switch_auroc`,
`gene_driver_truth_module_auroc`, `top10_switch_precision`,
`top10_truth_module_precision`, `switch_strength_auroc`, `opposite_polarity_rate`,
`switch_driver_recall`, `abundance_driver_recall`.

Transcript-level switch-event accuracy (vs `truth_switch_event`, only meaningful
in the `multi_isoform_switch` scenario where genes have ≥3 isoforms — at 2 isoforms
the switch transcript is degenerate):

- `switch_transcript_top1_accuracy` — fraction of switching genes whose highest-|r|
  transcript (the interpretation's predicted switch driver) equals the ground-truth
  driver transcript T1. Chance ≈ 1/n_transcripts (0.25 at n_tx=4).
- `switch_strength_auroc` — AUROC of the interpretation's per-gene `switch_strength`
  for separating switching from non-switching genes within a module (chance 0.5).
  Undefined (NaN) for modules that are all-switching or all-background, so it is scored
  on the subset of modules that contain both classes.
- `switch_magnitude_spearman` — Spearman correlation between the interpretation's
  per-gene `switch_strength` and the true switch magnitude (`true_delta_psi`). **Computed
  but not plotted** — see the construction caveat below.

### Switch-event construction, and why magnitude is not a reported readout

Each switching gene's PSI signal is `p1 = sigmoid(signal)`, where `signal` is driven by a
**unit-variance** module latent (`synthetic_data.py`, `module_latent ~ N(0, 1)`), and the
ground-truth magnitude is the across-sample range
`true_delta_psi = p1.max() − p1.min()` (`synthetic_data.py:282`). Because the latent has
the same unit scale for every gene, over ~160 samples its extremes reach roughly ±3, so
**every** switching gene saturates to `Δψ ≈ 0.93` (empirically mean 0.93, **SD 0.046**,
~11 % pinned at the 0.9998 clip). The generator has **no per-gene switch-magnitude
parameter** — switches are near-complete by construction.

Consequently the ground truth has **no magnitude gradient to recover**, and
`switch_magnitude_spearman` is ≈ 0 by construction for *any* method and *any* reasonable
predicted-strength proxy — it reflects sampling noise, not calibration (its per-module
estimate's SD ≈ 0.29 is exactly the n≈12 sampling SE of a correlation, centred at ~0.05).
We therefore keep the metric in the tables for completeness but report
`switch_transcript_top1_accuracy` (the switch *driver-identity* claim) and
`switch_strength_auroc` (does predicted strength rank switching genes above background) as
the two interpretation-accuracy panels (figS11). Making magnitude calibration testable
would require a per-gene amplitude parameter in the generator (e.g. scale `module_latent`
per gene so `Δψ` spans a wide range); that is a scenario-scoped change and is not currently
wired.

The `multi_isoform_switch` scenario (`configs/synthetic_grid.yaml`,
`n_transcripts_per_gene=4`) exists to exercise these; it is hashed on
`n_transcripts_per_gene` so its datasets are distinct from the 2-isoform grid.
