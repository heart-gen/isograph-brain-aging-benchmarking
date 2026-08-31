#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=module-meta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16  # 16 x 2000M = 32G; per-gene SVD reconstruction over the half-fits
#SBATCH --time=02:00:00
#SBATCH --output=03_module_trust/_m/logs/%x-%j.log

# Step 0 (post-hoc, no VAE re-fit): reconstruct per-half module eigengene Age effect +
# driver transcripts for the split-half partitions. SVD-heavy -> SLURM only. Override
# (cohort, region, method) via COHORT/REGION/METHOD env.

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
mkdir -p 03_module_trust/_m/logs

COHORT="${COHORT:-brainseq}"
REGION="${REGION:-caudate}"
METHOD="${METHOD:-isograph}"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** module meta: ${COHORT}/${REGION}/${METHOD} starts ****"
python -u -m isograph_benchmark.real_data.module_trust meta \
    --cohort "${COHORT}" --region "${REGION}" --method "${METHOD}" --k 5
conda deactivate
log_message "**** Complete ****"
