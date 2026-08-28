#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=scz-confound
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=real_data/_m/logs/scz-confound-%j.log

## Clinical/technical confound sensitivity for the SCZD caudate module findings.
##   audit  -- availability table only (what BrainSEQ actually releases)
##   run    -- audit + the full measured/proxy/combined ladder (default)
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p real_data/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
log "scz_confound_sensitivity ${*:-run}"
"${PY}" -m isograph_benchmark.real_data.scz_confound_sensitivity "${@:-run}"
log "**** complete ****"
