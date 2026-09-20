#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=switch-feature-refit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=24
#SBATCH --array=1-14
#SBATCH --time=04:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/switch-feature-refit-%A_%a.log

## Full IsoGraph refit per preprocessing setting (the follow-on the fixed-partition harness,
## _h/01h and _h/01i, names as its limitation). One production fit per array task: the
## published setting (task 1, the noise floor) and the 13 perturbed settings of the
## pseudocount, expression-filter and minor-isoform axes. The rollup is its own step,
##
## **Keep --array in step with the grid.** It said 1-12 while the grid held 14 settings, so
## the two minor-isoform settings were never re-run on the switching filter and silently kept
## results computed against the legacy committed fit. `--list` prints the grid and its size;
## the aggregate step now refuses to mix rows compared against different committed fits.
## 06_switch_mechanism/_h/02e.switch_feature_refit_aggregate.sh, after this array.
##
## Defaults to brainseq/caudate, the region the harness was defined on; override with
## COHORT/REGION env and list the grid (and its size) with `--list`. A real-data fit peaks
## near 30 GB, so 24 x 2000MB = 48G; do NOT pass --mem. Calls the env interpreter directly
## rather than `module load` + `conda activate`, so an array task cannot die on
## "module: command not found".
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

log "switch_feature_refit ${COHORT}/${REGION} index ${SLURM_ARRAY_TASK_ID:-1}"
"${PY}" -u -m isograph_benchmark.real_data.switch_feature_refit \
    --cohort "${COHORT}" --region "${REGION}" --index "${SLURM_ARRAY_TASK_ID:-1}" "$@"
log "**** complete ****"
