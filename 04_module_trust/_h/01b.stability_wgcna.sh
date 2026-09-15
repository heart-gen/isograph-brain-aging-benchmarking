#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=stability-wgcna
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --array=1-6
#SBATCH --output=04_module_trust/_m/logs/%x-%A_%a.log

# WGCNA within-cohort split-half stability (abundance-network reference ceiling) for
# the same 6 cohort-regions. Writes per-half partitions to the shared partitions/
# dir for the uniform Python aggregator.

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
mkdir -p 04_module_trust/_m/logs

SPECS=(
    "brainseq caudate"
    "brainseq hippocampus"
    "brainseq dlpfc"
    "gtex caudate_basal_ganglia"
    "gtex hippocampus"
    "gtex frontal_cortex_ba9"
)
spec="${SPECS[$((SLURM_ARRAY_TASK_ID - 1))]}"
read -r COHORT REGION <<< "${spec}"

SEEDS="${STABILITY_SEEDS:-5}"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** WGCNA split-half: ${COHORT}/${REGION} (${SEEDS} seeds) starts ****"
Rscript 04_module_trust/_h/stability_wgcna.R "${COHORT}" "${REGION}" "${SEEDS}"
conda deactivate
log_message "**** Complete ****"
