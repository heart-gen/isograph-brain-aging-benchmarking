#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-celltype-composition
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-4
#SBATCH --time=02:00:00
#SBATCH --output=02_module_discovery/brainseq/_m/logs/celltype-composition-%A_%a.log

# Cell-type composition confound test (BrainSEQ; reuses the committed MuSiC
# deconvolution from ../sex_context_brain).
#
# Per region:
#   1. celltype_composition  -> writes celltype_fractions.parquet (BrNum join) +
#      celltype_composition/marker_enrichment.parquet (module marker-depletion cut).
#   2. incremental_association --composition  -> repeats the de-confounded switch-vs-
#      abundance test with cell-type fractions added as inference covariates, writing to
#      isograph_vae/incremental_association_composition/. Compare against the canonical
#      isograph_vae/incremental_association/ for the with-vs-without contrast.
#
#   task 1 -> brainseq-sczd (caudate_sczd, diagnosis)
#   task 2 -> brainseq-aging caudate      (df=3 spline age)
#   task 3 -> brainseq-aging hippocampus
#   task 4 -> brainseq-aging dlpfc
#
# Requires modules.parquet + feature_scores.parquet (canonical fit) and the canonical
# incremental_association/ already written (08.incremental_association.sh).

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 02_module_discovery/brainseq/_m/logs
log_message "**** Cell-type composition (task ${SLURM_ARRAY_TASK_ID:-1}) ****"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

comp() { python -m isograph_benchmark.real_data.celltype_composition fractions "$@"; }
incr() { python -m isograph_benchmark.real_data.incremental_association "$@" --composition; }

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
case "${TASK_ID}" in
    1) comp brainseq-sczd;                       incr brainseq-sczd ;;
    2) comp brainseq-aging --region caudate;     incr brainseq-aging --region caudate ;;
    3) comp brainseq-aging --region hippocampus; incr brainseq-aging --region hippocampus ;;
    4) comp brainseq-aging --region dlpfc;       incr brainseq-aging --region dlpfc ;;
    *) echo "ERROR: unknown array task id ${TASK_ID} (expected 1-4)"; exit 1 ;;
esac
# After all array tasks finish, roll up the with-vs-without contrast:
#   python -m isograph_benchmark.real_data.celltype_composition meta

conda deactivate
log_message "**** Complete ****"
