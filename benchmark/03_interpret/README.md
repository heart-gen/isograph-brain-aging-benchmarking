# Synthetic Module Interpretation

This stage evaluates IsoGraph module interpretation on the synthetic benchmark
using ground-truth switching genes and module assignments.

Run:

```bash
sbatch benchmark/03_interpret/_h/step_1_interpret.sh
```

Main outputs:

- `benchmark/03_interpret/_m/synthetic_interpret_results.parquet`
- `benchmark/03_interpret/_m/synthetic_interpret_module_metrics.parquet`
- `benchmark/03_interpret/_m/synthetic_interpret_summary.parquet`

Large per-run `explain-module` outputs are written under
`benchmark/03_interpret/_o/runs/` and are intentionally ignored by git.

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
- `switch_magnitude_spearman` — Spearman correlation between the interpretation's
  per-gene `switch_strength` and the true switch magnitude (`true_delta_psi`).

The `multi_isoform_switch` scenario (`configs/synthetic_grid.yaml`,
`n_transcripts_per_gene=4`) exists to exercise these; it is hashed on
`n_transcripts_per_gene` so its datasets are distinct from the 2-isoform grid.
