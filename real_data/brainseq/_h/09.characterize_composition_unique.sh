#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-char-comp-unique
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-4
#SBATCH --time=03:00:00
#SBATCH --output=real_data/brainseq/_m/logs/char-comp-unique-%A_%a.log

# Characterize the composition-unique gene sets (GO:BP enrichment + IsoGraph
# module concentration) from the de-confounded gene-level test. Writes to
# real_data/brainseq/<region>/_m/isograph_vae/composition_unique/.
#
# DEPENDENCY: this consumes incremental_association/gene_level.parquet, so it
# MUST afterok on 08.incremental_association.sh -- NOT 04.interpret_modules.sh.
# (Mis-wiring to interpret is what cancelled the S3 cascade on 2026-06-25: the
# array launched before incremental finished writing. The consumer now also
# bounded-waits for gene_level.parquet to absorb afterok skew, but wire the
# correct parent regardless.)
#
#   task 1 -> brainseq-sczd
#   task 2 -> brainseq-aging caudate
#   task 3 -> brainseq-aging hippocampus
#   task 4 -> brainseq-aging dlpfc
#
# After the array, build the cross-analysis overlap with:
#   python -m isograph_benchmark.real_data.characterize_composition_unique --combine
# Flags forwarded after the script name, e.g. --variant with-abundance / --no-go.

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/brainseq/_m/logs
log_message "**** Characterize composition-unique (task ${SLURM_ARRAY_TASK_ID:-1}) ****"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
case "${TASK_ID}" in
    1) python -m isograph_benchmark.real_data.characterize_composition_unique brainseq-sczd "$@" ;;
    2) python -m isograph_benchmark.real_data.characterize_composition_unique brainseq-aging --region caudate "$@" ;;
    3) python -m isograph_benchmark.real_data.characterize_composition_unique brainseq-aging --region hippocampus "$@" ;;
    4) python -m isograph_benchmark.real_data.characterize_composition_unique brainseq-aging --region dlpfc "$@" ;;
    *) echo "ERROR: unknown array task id ${TASK_ID} (expected 1-4)"; exit 1 ;;
esac

conda deactivate
log_message "**** Complete ****"
