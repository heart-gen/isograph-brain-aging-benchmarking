#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-iso-aging
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=32
#SBATCH --array=1-3
#SBATCH --time=04:00:00
#SBATCH --output=real_data/brainseq/_m/logs/%x-%A_%a.log
# Run IsoGraph VAE on BrainSEQ adult control aging bundles:
# caudate, hippocampus, and DLPFC.
set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"
mkdir -p real_data/brainseq/_m/logs

log_message "**** BrainSEQ aging IsoGraph job starts ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID:-local}"
echo "Job name: ${SLURM_JOB_NAME:-brainseq-iso-aging}"
echo "Node name: ${SLURM_NODENAME:-local}"
echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "Activating IsoGraph environment"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

REGIONS=(caudate hippocampus dlpfc)
RUN_ARGS=("$@")
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    REGION="${REGIONS[$((SLURM_ARRAY_TASK_ID - 1))]}"
    RUN_ARGS=(--region "${REGION}" "$@")
    log_message "Running region ${REGION}"
fi

python -m isograph_benchmark.real_data.run_models brainseq-aging "${RUN_ARGS[@]}"

conda deactivate
log_message "**** BrainSEQ aging IsoGraph job ends ****"
