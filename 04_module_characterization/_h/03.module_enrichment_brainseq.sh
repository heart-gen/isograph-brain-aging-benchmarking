#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-module-enrich
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-4
#SBATCH --time=03:00:00
#SBATCH --output=04_module_characterization/_m/logs/module-enrich-%A_%a.log

# Module-level GO:BP enrichment + network metrics for every partition present in the
# store (IsoGraph, classical WGCNA, and the matched-feature WGCNA baselines). Writes
# 02_module_discovery/brainseq/<region>/_m/module_enrichment/.
#
#   task 1 -> brainseq-sczd        (IsoGraph + classical WGCNA; no matched baselines)
#   task 2 -> brainseq-aging caudate
#   task 3 -> brainseq-aging hippocampus
#   task 4 -> brainseq-aging dlpfc
#
# Defaults to --method all. all_modules.parquet and summary.json are rewritten from the
# methods passed on THIS run, so the CLI's own default (both = isograph + wgcna) would
# silently drop the matched-baseline rows that stage 04's baseline comparison reads.
# Methods whose modules.parquet is absent are skipped, so `all` is safe for caudate_sczd.
#
# Requires modules.parquet (+ edges.parquet for IsoGraph network metrics) and the
# saved trait/diagnosis tables. Flags forwarded after the script name, e.g.
# --variant with-abundance / --method isograph / --no-go.
#
# Calls the env interpreter directly rather than `module load` + `conda activate`, so an
# array task cannot die on "module: command not found".

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 04_module_characterization/_m/logs
log_message "**** Module enrichment + network (task ${SLURM_ARRAY_TASK_ID:-1}) ****"

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

ARGS=("$@")
[[ " $* " == *" --method "* ]] || ARGS+=(--method all)

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
case "${TASK_ID}" in
    1) "${PY}" -m isograph_benchmark.real_data.module_enrichment brainseq-sczd "${ARGS[@]}" ;;
    2) "${PY}" -m isograph_benchmark.real_data.module_enrichment brainseq-aging --region caudate "${ARGS[@]}" ;;
    3) "${PY}" -m isograph_benchmark.real_data.module_enrichment brainseq-aging --region hippocampus "${ARGS[@]}" ;;
    4) "${PY}" -m isograph_benchmark.real_data.module_enrichment brainseq-aging --region dlpfc "${ARGS[@]}" ;;
    *) echo "ERROR: unknown array task id ${TASK_ID} (expected 1-4)"; exit 1 ;;
esac

log_message "**** Complete ****"
