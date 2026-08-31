#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=isograph-partition-backfill
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --array=0-7
#SBATCH --output=01_synthetic_benchmark/03_metrics/_m/logs/partition-backfill-%A_%a.log

## Step 0: backfill the fragmentation-sensitive partition metrics (reviewer items 1 and 4)
## onto every completed synthetic run.  Recomputes ARI / AMI / homogeneity / completeness /
## V-measure and the module-count-preserving null for best-match Jaccard from each run's
## stored modules.parquet plus its dataset's truth tables — NO model is re-fitted.
##
## New keys are merged into telemetry.json / done.json, so step_1_collect picks them up with
## no change.  Each shard asserts that the recomputed best-match Jaccard equals the stored
## metrics.module_recovery, and fails the job if any run disagrees.
##
## Run BEFORE step_1_collect.sh.  Usage:
##   sbatch 01_synthetic_benchmark/03_metrics/_h/step_0_backfill_partition_metrics.sh
set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 01_synthetic_benchmark/03_metrics/_m/logs

N_SHARDS="${SLURM_ARRAY_TASK_COUNT:-8}"
SHARD="${SLURM_ARRAY_TASK_ID:-0}"

log_message "**** Step 0: partition-metric backfill (shard ${SHARD}/${N_SHARDS}) ****"
echo "User: ${USER}"
echo "Host: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1

log_message "Activating isograph environment"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

python -m isograph_benchmark.benchmark.backfill_metrics \
  --run-root 01_synthetic_benchmark/01_synthetic/_o/runs \
  --out-dir 01_synthetic_benchmark/01_synthetic/_m/partition_backfill \
  --shard "${SHARD}" \
  --n-shards "${N_SHARDS}" \
  --n-perm 200 \
  --seed 13 \
  "$@"

conda deactivate
log_message "**** Step 0 complete (shard ${SHARD}) ****"
