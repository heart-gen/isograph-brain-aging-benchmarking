#!/usr/bin/env bash
# Run benchmark batches locally using GNU parallel (no SLURM required).
#
# Usage:
#   ./benchmark/01_synthetic/_h/run_local.sh <resource_class> [n_parallel]
#
# resource_class: cpu_short | vae | wgcna_cpu | scale
# n_parallel:     number of concurrent batches (default: 2)
#
# Examples:
#   ./benchmark/01_synthetic/_h/run_local.sh cpu_short 4
#   ./benchmark/01_synthetic/_h/run_local.sh vae 2
#   ./benchmark/01_synthetic/_h/run_local.sh wgcna_cpu 1
#   ./benchmark/01_synthetic/_h/run_local.sh scale 1

set -euo pipefail

RESOURCE_CLASS="${1:?Usage: $0 <resource_class> [n_parallel]}"
N_PARALLEL="${2:-2}"

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BATCH_TSV="${REPO_ROOT}/benchmark/01_synthetic/_m/batches_${RESOURCE_CLASS}.tsv"
BATCH_PARQUET="${REPO_ROOT}/benchmark/01_synthetic/_m/synthetic_batches.parquet"
LOG_DIR="${REPO_ROOT}/benchmark/01_synthetic/_m/logs"

mkdir -p "${LOG_DIR}"

if [[ ! -f "${BATCH_TSV}" ]]; then
    echo "ERROR: batch file not found: ${BATCH_TSV}" >&2
    exit 1
fi

echo "Running ${RESOURCE_CLASS} batches with N_PARALLEL=${N_PARALLEL}"
echo "Batch file: ${BATCH_TSV}"
echo "Logs: ${LOG_DIR}"

tail -n +2 "${BATCH_TSV}" | awk '{print $1}' | \
    parallel --jobs "${N_PARALLEL}" \
             --joblog "${LOG_DIR}/parallel_${RESOURCE_CLASS}.log" \
    "python -m isograph_benchmark.benchmark.run_batch \
        --batch-file '${BATCH_PARQUET}' \
        --batch-id {} \
        > '${LOG_DIR}/batch_${RESOURCE_CLASS}_{}.out' \
        2> '${LOG_DIR}/batch_${RESOURCE_CLASS}_{}.err'"

echo "Done. Check ${LOG_DIR}/parallel_${RESOURCE_CLASS}.log for per-batch status."
