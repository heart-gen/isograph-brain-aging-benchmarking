#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=module-anchoring-meta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%j.log
## Roll the per-analysis module-level genetic anchoring (04) into
## _m/module_genetic_anchoring_meta/ (Table S16). Gate on the 04 array:
##   sbatch --dependency=afterok:<04> 05_genetic_anchoring/_h/02c.module_anchoring_meta.sh
## Calls the env interpreter directly, so the job cannot die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 05_genetic_anchoring/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "**** module_genetic_anchoring --meta ****"
"${PY}" -u -m isograph_benchmark.real_data.module_genetic_anchoring --meta "$@"
log "**** complete ****"
