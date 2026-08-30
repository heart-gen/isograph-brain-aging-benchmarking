#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-leiden-sweep-abund
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-4
#SBATCH --time=06:00:00
#SBATCH --output=02_module_discovery/brainseq/_m/logs/leiden-sweep-abund-%A_%a.log

# Same Leiden-resolution sweep as 05.sweep_leiden.sh, but operating on the
# WITH-ABUNDANCE refit (isograph_vae_with_abundance/) rather than the switch-only
# baseline. This selects each with-abundance graph's resolution by the SAME GO
# n_enriched criterion used for the baseline, so the two IsoGraph variants are
# compared at consistently-selected resolutions (the with-abundance graph is
# denser, so equal resolution != equal module count).
#
# Runs as a 4-task ARRAY (one analysis/region per task):
#   task 1 -> brainseq-sczd (caudate_sczd)
#   task 2 -> brainseq-aging caudate
#   task 3 -> brainseq-aging hippocampus
#   task 4 -> brainseq-aging dlpfc
#
# Requires the with-abundance runs (06.run_isograph_with_abundance.sh) to have
# written edges.parquet + feature_scores.parquet to
# 02_module_discovery/brainseq/<region>/_m/isograph_vae_with_abundance/.
#
# Flags forwarded to the sweep (pass after the script name), e.g. --write-best.
# --variant with-abundance is supplied automatically below.

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 02_module_discovery/brainseq/_m/logs
log_message "**** BrainSeq Leiden sweep [with-abundance] (task ${SLURM_ARRAY_TASK_ID:-1}) ****"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
case "${TASK_ID}" in
    1) log_message "Sweeping SCZD caudate [with-abundance]"
       python -m isograph_benchmark.real_data.sweep_leiden brainseq-sczd --variant with-abundance "$@" ;;
    2) log_message "Sweeping aging caudate [with-abundance]"
       python -m isograph_benchmark.real_data.sweep_leiden brainseq-aging --region caudate --variant with-abundance "$@" ;;
    3) log_message "Sweeping aging hippocampus [with-abundance]"
       python -m isograph_benchmark.real_data.sweep_leiden brainseq-aging --region hippocampus --variant with-abundance "$@" ;;
    4) log_message "Sweeping aging dlpfc [with-abundance]"
       python -m isograph_benchmark.real_data.sweep_leiden brainseq-aging --region dlpfc --variant with-abundance "$@" ;;
    *) echo "ERROR: unknown array task id ${TASK_ID} (expected 1-4)"; exit 1 ;;
esac

conda deactivate
log_message "**** Complete ****"
