#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=extract-tin
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4  # 4 x 2000M = 8G; one TIN file in memory at a time + ~0.2G matrix
#SBATCH --time=02:00:00
#SBATCH --output=real_data/stability/_m/logs/%x-%j.log

# Parse per-sample RSeQC TIN tables into a bundle-aligned transcript x sample matrix
# (cached parquet) + per-sample median-TIN sidecar. Heavy parse (hundreds of files x
# ~388k rows) -> SLURM only. Override the (cohort, region) via COHORT/REGION env.

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/stability/_m/logs

COHORT="${COHORT:-brainseq}"
REGION="${REGION:-caudate}"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** TIN extract: ${COHORT}/${REGION} starts ****"
python -m isograph_benchmark.real_data.tin extract --cohort "${COHORT}" --region "${REGION}"
conda deactivate
log_message "**** Complete ****"
