#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=refresh-sample-tables
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --array=0-7
#SBATCH --time=04:00:00
#SBATCH --output=01_synthetic_benchmark/01_synthetic/_m/logs/%x-%A_%a.log

## Rewrites the cached synthetic `samples.parquet` files so that the covariates
## isograph_vae_residual regresses out (RIN, neuron_frac, batch, library_size) are actually
## present. 1730 of the cached datasets predate those columns, which makes residualization a
## SILENT NO-OP -- the residual arm reduces to plain isograph_vae, and the open "is
## residualization free on unconfounded data?" ablation would answer itself.
##
## Only samples.parquet is written. Each dataset is fully regenerated so its count matrices,
## feature tables and truth tables can be fingerprinted against the cached copies, and the
## regenerated arrays are then DISCARDED: 13,060 completed runs key off the cached expression
## data, so a single non-reproducing dataset must fail loudly rather than overwrite results
## already in the figures.
##
## Each task re-runs the single-dataset probe first and exits without writing if it fails.
## Exit 1 = at least one dataset did not reproduce (a reproducibility finding -- escalate,
## do not re-run). Exit 2 = the probe itself failed.
##
## Usage:
##   python -m isograph_benchmark.benchmark.refresh_sample_tables --probe    # seconds, no writes
##   sbatch 01_synthetic_benchmark/01_synthetic/_h/refresh_sample_tables.sh
##   python -m isograph_benchmark.benchmark.refresh_sample_tables --collect  # after all shards

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 01_synthetic_benchmark/01_synthetic/_m/logs

N_SHARDS=8
SHARD="${SLURM_ARRAY_TASK_ID:-0}"

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** refresh sample tables: shard ${SHARD}/${N_SHARDS} ****"
python -u -m isograph_benchmark.benchmark.refresh_sample_tables \
    --shard "${SHARD}" --n-shards "${N_SHARDS}" "$@"
status=$?
conda deactivate
log_message "**** Complete (shard ${SHARD}, exit ${status}) ****"
exit "${status}"
