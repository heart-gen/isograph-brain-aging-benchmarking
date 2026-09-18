#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=dtu-added-value-analyze
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --array=1-6
#SBATCH --time=02:00:00
#SBATCH --output=04_module_trust/_m/logs/%x-%A_%a.log

## Module context and held-out DTU evidence (PI review item 12a), step 2 of 3: the held-out test for one
## region -- 5 seeds x 2 directions, each with the stratified and plain permutation nulls
## (1,000 each) and the abundance-neighbour joint model. Reads the 02e satuRn halves and the
## 01a split-half partitions and graphs.
##
## Array: the 6 split-half regions. Bridges memory is --cpus-per-task x 2000MB (8G); do NOT
## pass --mem. Extra arguments are forwarded (e.g. --n-perm 200 for a quick look).
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 04_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

REGIONS=(brainseq:caudate brainseq:hippocampus brainseq:dlpfc
         gtex:caudate_basal_ganglia gtex:hippocampus gtex:frontal_cortex_ba9)
CR="${REGIONS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}"
log "module_dtu_added_value analyze --cohort ${CR%%:*} --region ${CR##*:}"
"${PY}" -u -m isograph_benchmark.real_data.module_dtu_added_value analyze \
    --cohort "${CR%%:*}" --region "${CR##*:}" "$@"
log "**** complete ****"
