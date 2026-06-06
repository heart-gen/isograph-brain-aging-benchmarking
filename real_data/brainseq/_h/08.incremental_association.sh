#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-incremental-assoc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-4
#SBATCH --time=02:00:00
#SBATCH --output=real_data/brainseq/_m/logs/incremental-assoc-%A_%a.log

# Switch-vs-abundance incremental association (strengths/limitations of IsoGraph).
# Runs the gene-level DE-CONFOUNDED test + the module-level incremental contrast
# on the saved IsoGraph features, writing to
# real_data/brainseq/<region>/_m/isograph_vae/incremental_association/.
#
#   task 1 -> brainseq-sczd (caudate_sczd, diagnosis)
#   task 2 -> brainseq-aging caudate      (df=3 spline age)
#   task 3 -> brainseq-aging hippocampus
#   task 4 -> brainseq-aging dlpfc
#
# Requires feature_scores.parquet + modules.parquet (01/02 runs + Part 1 sweep).
# Flags forwarded after the script name, e.g. --variant with-abundance.

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/brainseq/_m/logs
log_message "**** Incremental association (task ${SLURM_ARRAY_TASK_ID:-1}) ****"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
case "${TASK_ID}" in
    1) python -m isograph_benchmark.real_data.incremental_association brainseq-sczd "$@" ;;
    2) python -m isograph_benchmark.real_data.incremental_association brainseq-aging --region caudate "$@" ;;
    3) python -m isograph_benchmark.real_data.incremental_association brainseq-aging --region hippocampus "$@" ;;
    4) python -m isograph_benchmark.real_data.incremental_association brainseq-aging --region dlpfc "$@" ;;
    *) echo "ERROR: unknown array task id ${TASK_ID} (expected 1-4)"; exit 1 ;;
esac

conda deactivate
log_message "**** Complete ****"
