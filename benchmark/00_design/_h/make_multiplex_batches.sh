#!/usr/bin/env bash
# Regenerate synthetic grid and create multiplex-only batch files for Bridges-2.
# Run this locally before submitting run_batch_multiplex.sh.
# Outputs: benchmark/01_synthetic/_m/multiplex_batches.parquet
#          benchmark/01_synthetic/_m/multiplex_batches_vae.tsv  (array task list)
set -euo pipefail

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${PWD}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: run from the isograph-brain-aging-benchmarking repo root."
    exit 1
fi
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

echo "=== Step 1: Regenerate full synthetic grid ==="
<<<<<<< HEAD
python -m isograph_benchmark.benchmark.run_synthetic
=======
python -m isograph_benchmark.benchmark.build_synthetic_grid
>>>>>>> a0fe94a (Add abundance+switch (multiplex) mode to benchmark and real data pipeline.)

echo ""
echo "=== Step 2: Create multiplex-only batch files ==="
python -m isograph_benchmark.benchmark.make_batches \
  --grid benchmark/00_design/_m/synthetic_run_grid.parquet \
  --out benchmark/01_synthetic/_m/multiplex_batches.parquet \
  --method-filter isograph_vae_multiplex \
  --scenario-filter abundance_switch_mixed \
  --scenario-filter idealized_switching \
  --scenario-filter unequal_isoform_abundance \
  --scenario-filter non_switching_background \
  --scenario-filter noise_stress \
  --scenario-filter feature_space_interactions

echo ""
echo "=== Done ==="
echo "Update --array upper bound in run_batch_multiplex.sh to match the vae batch count above."
echo "Then submit: sbatch benchmark/01_synthetic/_h/run_batch_multiplex.sh"
