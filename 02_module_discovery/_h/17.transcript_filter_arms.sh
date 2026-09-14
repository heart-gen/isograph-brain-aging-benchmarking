#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=transcript-filter-arms
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=24
#SBATCH --array=1-9
#SBATCH --time=03:00:00
#SBATCH --output=02_module_discovery/_m/logs/transcript-filter-arms-%A_%a.log

## BrainSEQ transcript-filter arms: the production IsoGraph fit re-run with a switching (usage)
## filter, beside production-setting refits (noise floors), for caudate / hippocampus / DLPFC
## aging and the SCZD caudate fit. Writes ONLY under 02_module_discovery/_m/transcript_filter_arms/;
## the production stores are never written. One fit per array task; list with `--list`. Then:
##
##   a=$(sbatch --parsable 02_module_discovery/_h/17.transcript_filter_arms.sh)
##   sbatch --dependency=afterany:${a} --array=1 --export=ALL,FILTER_ARMS_COMPARE=1 \
##       02_module_discovery/_h/17.transcript_filter_arms.sh
##
## Seed floor (the same 9 arms at seeds 14 and 15; compare picks them up automatically):
##
##   s=$(sbatch --parsable --array=1-18 --export=ALL,FILTER_ARMS_GRID=seed \
##       02_module_discovery/_h/17.transcript_filter_arms.sh)
##   sbatch --dependency=afterany:${s} --array=1 --export=ALL,FILTER_ARMS_COMPARE=1 \
##       02_module_discovery/_h/17.transcript_filter_arms.sh
##
## A real-data fit peaks near 30 GB (DLPFC 48 GB): 24 x 2000MB = 48G; do NOT pass --mem.
## Calls the env interpreter directly, so a task cannot die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 02_module_discovery/_m/logs
# As the production wrappers: cap BLAS threads (per-thread workspace over the SVD calls).
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
if [[ "${FILTER_ARMS_COMPARE:-0}" == "1" ]]; then
    log "transcript_filter_arms compare"
    "${PY}" -u -m isograph_benchmark.real_data.transcript_filter_arms compare "$@"
else
    log "transcript_filter_arms fit grid ${FILTER_ARMS_GRID:-filter} index ${SLURM_ARRAY_TASK_ID:-1}"
    "${PY}" -u -m isograph_benchmark.real_data.transcript_filter_arms fit \
        --grid "${FILTER_ARMS_GRID:-filter}" --index "${SLURM_ARRAY_TASK_ID:-1}" "$@"
fi
log "**** complete ****"
