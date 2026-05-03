#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-drd2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=real_data/brainseq/_m/logs/%x-%j.log

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "ERROR: project root does not exist: ${PROJECT_ROOT}"
    exit 1
fi
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the isograph-brain-aging-benchmarking repo root or set ISOGRAPH_BENCHMARK_ROOT."
    echo "Current project root candidate: ${PROJECT_ROOT}"
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/brainseq/_m/logs

log_message "**** BrainSEQ DRD2 case study job starts ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID:-local}"
echo "Job name: ${SLURM_JOB_NAME:-brainseq-drd2}"
echo "Node name: ${SLURM_NODENAME:-local}"
echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "Activating R environment"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

Rscript real_data/brainseq/_h/03.drd2_case_study.R "$@"

conda deactivate
log_message "**** BrainSEQ DRD2 case study job ends ****"
