#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=GPU-shared
#SBATCH --gres=gpu:v100-16:1
#SBATCH --job-name=brainseq-iso-aging
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --array=1-3
#SBATCH --time=02:00:00
#SBATCH --output=real_data/brainseq/_m/logs/%x-%A_%a.log
# Run IsoGraph VAE on BrainSEQ adult control aging bundles:
# caudate, hippocampus, and DLPFC.
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

log_message "Checking Python analysis dependencies"
python -c "import isograph_benchmark, numpy, pandas, scipy, patsy"

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
