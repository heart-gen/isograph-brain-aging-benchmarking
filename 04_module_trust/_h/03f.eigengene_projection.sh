#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=eigengene-projection
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16
#SBATCH --array=1-6
#SBATCH --time=06:00:00
#SBATCH --output=04_module_trust/_m/logs/%x-%A_%a.log

## Cross-cohort eigengene projection (PI review 03-4): freeze each trusted module's eigengene
## weights in one cohort, project them onto the other, and test preservation (signed kME vs a
## type- and size-matched null) and age-sign concordance, in both directions. IsoGraph switch
## axes are oriented across cohorts by shared-transcript loadings first.
##
## Array: {isograph, wgcna} x {caudate, hippocampus, dlpfc_ba9}. The cross-pair rollup is its own
## step, 04_module_trust/_h/04c.eigengene_projection_aggregate.sh, after this array.
##
## Needs the Q1 `stability` tables (trusted modules, 02c) and the production fits. Bridges memory
## is --cpus-per-task x 2000MB (32G); do NOT pass --mem. Calls the env interpreter directly.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 04_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

METHODS=(isograph wgcna)
PAIRS=(caudate hippocampus dlpfc_ba9)
IDX=$((${SLURM_ARRAY_TASK_ID:-1} - 1))
METHOD="${METHODS[$((IDX / 3))]}"
PAIR="${PAIRS[$((IDX % 3))]}"
log "eigengene_projection run --pair ${PAIR} --method ${METHOD}"
"${PY}" -u -m isograph_benchmark.real_data.eigengene_projection run \
    --pair "${PAIR}" --method "${METHOD}" "$@"
log "**** complete ****"
