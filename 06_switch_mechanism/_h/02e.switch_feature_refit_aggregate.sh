#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=switch-feature-refit-aggregate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/%x-%j.log

## Roll the per-setting refits (01j.switch_feature_refit.sh) into the refit-sensitivity report,
## against the published setting's noise floor. Same COHORT/REGION as the array (default
## brainseq/caudate):
##   sbatch --dependency=afterok:<01j> 06_switch_mechanism/_h/02e.switch_feature_refit_aggregate.sh
## Bridges memory is --cpus-per-task x 2000MB; do NOT pass --mem. Calls the env interpreter directly.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 06_switch_mechanism/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
COHORT="${COHORT:-brainseq}"
REGION="${REGION:-caudate}"

log "switch_feature_refit --aggregate ${COHORT}/${REGION}"
"${PY}" -u -m isograph_benchmark.real_data.switch_feature_refit \
    --cohort "${COHORT}" --region "${REGION}" --aggregate "$@"
log "**** complete ****"
