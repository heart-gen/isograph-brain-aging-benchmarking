#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=rep-perm-report
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=04_module_trust/_m/logs/%x-%j.log

## Collect every replication-permutation stats json that 03d wrote (all three --covariates modes)
## into REPLICATION_PERMUTATION.md. Pure reads, no permutation. Run after every 03d array:
##   sbatch --dependency=afterok:<03d full>:<03d complement>:<03d none> \
##       04_module_trust/_h/04b.replication_permutation_report.sh
## Calls the env interpreter directly, so the job cannot die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 04_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "**** replication permutation report ****"
"${PY}" -u -m isograph_benchmark.real_data.replication_permutation --report "$@"
log "**** complete ****"
