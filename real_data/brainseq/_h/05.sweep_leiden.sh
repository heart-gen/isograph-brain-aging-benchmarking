#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-leiden-sweep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-4
#SBATCH --time=06:00:00
#SBATCH --output=real_data/brainseq/_m/logs/leiden-sweep-%A_%a.log

# Re-clusters saved IsoGraph edges at a grid of Leiden resolutions WITHOUT
# refitting the VAE (reuses edges.parquet + feature_scores.parquet). Reports
# n_modules / giant fraction / GO-enrichment / trait-association counts at each
# resolution so we can pick a biologically interpretable module count.
#
# Runs as a 4-task ARRAY (one analysis/region per task) so each finishes well
# under the wall. GO enrichment costs ~8-13 min per resolution, so a single
# serial pass over all four analyses times out at 2h -- hence the split.
#   task 1 -> brainseq-sczd (caudate_sczd)
#   task 2 -> brainseq-aging caudate
#   task 3 -> brainseq-aging hippocampus
#   task 4 -> brainseq-aging dlpfc
#
# Requires the base runs (01.run_isograph_aging.sh / 02.run_isograph_sczd.sh)
# to have already written artifacts to real_data/brainseq/<region>/_m/.
#
# Flags forwarded to the sweep (pass after the script name):
#   --dry-run      compute + print only, write nothing
#   --write-best   commit the best-resolution artifacts (selected by GO density)
#   --no-go        skip GO enrichment (use trait-association fallback for --write-best)
#   --resolutions  override the default per-analysis resolution grid

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/brainseq/_m/logs
log_message "**** BrainSeq Leiden resolution sweep (task ${SLURM_ARRAY_TASK_ID:-1}) ****"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
case "${TASK_ID}" in
    1) log_message "Sweeping SCZD caudate"
       python -m isograph_benchmark.real_data.sweep_leiden brainseq-sczd "$@" ;;
    2) log_message "Sweeping aging caudate"
       python -m isograph_benchmark.real_data.sweep_leiden brainseq-aging --region caudate "$@" ;;
    3) log_message "Sweeping aging hippocampus"
       python -m isograph_benchmark.real_data.sweep_leiden brainseq-aging --region hippocampus "$@" ;;
    4) log_message "Sweeping aging dlpfc"
       python -m isograph_benchmark.real_data.sweep_leiden brainseq-aging --region dlpfc "$@" ;;
    *) echo "ERROR: unknown array task id ${TASK_ID} (expected 1-4)"; exit 1 ;;
esac

conda deactivate
log_message "**** Complete ****"
