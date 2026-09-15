#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=age-curvature
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --output=04_module_trust/_m/logs/age-curvature-%j.log

## Does the module Age trajectory actually bend?
##
## Decides which age model the cross-cohort replication count is quoted under. The linear
## and df=3 spline arms disagree (25/130 vs 15/130 at matched covariates); this test says
## whether that gap is a correction (real curvature the line misses) or a power tax (two
## df spent on noise). It reads the committed age_spline/age_linear tables and refits
## nothing, so it cannot drift from the numbers those tables carry.
##
## Cheap -- six small parquet reads. Memory on Bridges is --cpus-per-task x 2000MB, so
## 2 cpus = 4G; do NOT pass --mem.
##
## Extra arguments are forwarded, e.g.
##   sbatch 04_module_trust/_h/01d.age_model_curvature.sh --fdr 0.10
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 04_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "age_model_curvature $*"
"${PY}" -m isograph_benchmark.real_data.age_model_curvature "$@"
log "**** complete ****"
