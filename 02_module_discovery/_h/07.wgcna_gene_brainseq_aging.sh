#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=bs-wgcna-aging
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16
#SBATCH --array=1-3
#SBATCH --time=04:00:00
#SBATCH --output=02_module_discovery/_m/logs/%x-%A_%a.log

# Gene-level WGCNA aging baseline on the three BrainSEQ control regions
# (caudate, hippocampus, dlpfc). Provides the BrainSEQ WGCNA arm for the
# cross-cohort replication against GTEx.

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "ERROR: project root does not exist: ${PROJECT_ROOT}"
    exit 1
fi
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the isograph-brain-aging-benchmarking repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 02_module_discovery/_m/logs

log_message "**** BrainSEQ WGCNA aging job starts ****"
echo "Job id: ${SLURM_JOBID:-local} | Node: ${SLURM_NODENAME:-local} | Host: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
log_message "Activating R environment"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

REGIONS=(caudate hippocampus dlpfc)
RUN_ARGS=("$@")
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    REGION="${REGIONS[$((SLURM_ARRAY_TASK_ID - 1))]}"
    RUN_ARGS=("${REGION}" "$@")
    log_message "Running region ${REGION}"
fi

Rscript 02_module_discovery/_h/07.wgcna_gene_brainseq_aging.R "${RUN_ARGS[@]}"

conda deactivate
log_message "**** BrainSEQ WGCNA aging job ends ****"
