#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-wgcna-sczd
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH --output=real_data/brainseq/_m/logs/%x-%j.log

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"
mkdir -p real_data/brainseq/_m/logs

log_message "**** BrainSEQ SCZD WGCNA job starts ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID:-local}"
echo "Job name: ${SLURM_JOB_NAME:-brainseq-wgcna-sczd}"
echo "Node name: ${SLURM_NODENAME:-local}"
echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "Activating R environment"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
export ISOGRAPH_PYTHON="${ISOGRAPH_PYTHON:-/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python}"

Rscript real_data/brainseq/_h/02.run_wgcna_sczd.R "$@"

conda deactivate
log_message "**** BrainSEQ SCZD WGCNA job ends ****"
