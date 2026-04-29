#!/usr/bin/env bash
#SBATCH --job-name=isograph-synth
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=32
#SBATCH --array=1-1%50
#SBATCH --output=benchmark/01_synthetic/_m/logs/slurm-%A_%a.out
#SBATCH --error=benchmark/01_synthetic/_m/logs/slurm-%A_%a.err

set -euo pipefail

source ~/.venvs/isograph/bin/activate

BATCH_FILE="${BATCH_FILE:-benchmark/01_synthetic/_m/synthetic_batches.parquet}"
BATCH_INDEX="${SLURM_ARRAY_TASK_ID:-1}"

python -m isograph_benchmark.benchmark.run_batch \
  --batch-file "${BATCH_FILE}" \
  --batch-index "${BATCH_INDEX}"
