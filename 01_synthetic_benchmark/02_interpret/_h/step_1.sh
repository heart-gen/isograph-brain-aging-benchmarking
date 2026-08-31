#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=isograph-interpret
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=0-15
#SBATCH --time=02:00:00
#SBATCH --output=01_synthetic_benchmark/02_interpret/_m/logs/interpret-%A_%a.log

# Synthetic module interpretation benchmark, sharded across array tasks. Each
# task processes runs[shard::N_SHARDS] and writes a per-shard partial under
# _m/_shards/; step_2.sh (afterok) merges them and runs the bootstrap
# summary once. The expensive explain_artifact_modules output is checkpointed on
# disk per run, so reruns reuse it. Pass extra args (e.g. --method) after "$@".

set -euo pipefail

N_SHARDS=16

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Synthetic interpretation: shard ${SLURM_ARRAY_TASK_ID}/${N_SHARDS} ****"
mkdir -p 01_synthetic_benchmark/02_interpret/_m/logs

module purge
module load anaconda3/2024.10-1

log_message "Activating isograph environment"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

python -m isograph_benchmark.benchmark.interpret_modules \
    --shard "${SLURM_ARRAY_TASK_ID}" --n-shards "${N_SHARDS}" "$@"

conda deactivate
log_message "**** Complete ****"
