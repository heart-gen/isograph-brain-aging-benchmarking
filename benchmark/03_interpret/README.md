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
