#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=eigengene-projection-aggregate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=04_module_trust/_m/logs/%x-%j.log

## Roll the six eigengene-projection runs (03f: {isograph, wgcna} x 3 region pairs) into the
## cross-pair summary. Reads what 03f wrote; run after that array:
##   sbatch --dependency=afterok:<03f> 04_module_trust/_h/04c.eigengene_projection_aggregate.sh
## Bridges memory is --cpus-per-task x 2000MB; do NOT pass --mem. Calls the env interpreter directly.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 04_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "eigengene_projection aggregate"
"${PY}" -u -m isograph_benchmark.real_data.eigengene_projection aggregate "$@"
log "**** complete ****"
