#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-module-enrich
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-4
#SBATCH --time=03:00:00
#SBATCH --output=02_module_discovery/brainseq/_m/logs/module-enrich-%A_%a.log

# Module-level GO:BP enrichment + network metrics for BOTH IsoGraph and WGCNA
# partitions -- the module-network narrative. Writes
# 02_module_discovery/brainseq/<region>/_m/module_enrichment/. WGCNA is only present for
# caudate_sczd (task 1); aging tasks characterize IsoGraph modules only until an
# aging WGCNA exists.
#
#   task 1 -> brainseq-sczd        (IsoGraph + WGCNA)
#   task 2 -> brainseq-aging caudate
#   task 3 -> brainseq-aging hippocampus
#   task 4 -> brainseq-aging dlpfc
#
# Requires modules.parquet (+ edges.parquet for IsoGraph network metrics) and the
# saved trait/diagnosis tables. Flags forwarded after the script name, e.g.
# --variant with-abundance / --method isograph / --no-go.

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 02_module_discovery/brainseq/_m/logs
log_message "**** Module enrichment + network (task ${SLURM_ARRAY_TASK_ID:-1}) ****"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
case "${TASK_ID}" in
    1) python -m isograph_benchmark.real_data.module_enrichment brainseq-sczd "$@" ;;
    2) python -m isograph_benchmark.real_data.module_enrichment brainseq-aging --region caudate "$@" ;;
    3) python -m isograph_benchmark.real_data.module_enrichment brainseq-aging --region hippocampus "$@" ;;
    4) python -m isograph_benchmark.real_data.module_enrichment brainseq-aging --region dlpfc "$@" ;;
    *) echo "ERROR: unknown array task id ${TASK_ID} (expected 1-4)"; exit 1 ;;
esac

conda deactivate
log_message "**** Complete ****"
