#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=module-sizes
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=02_module_discovery/_m/logs/%x-%j.log

## Module-size tables for every method in every cohort x region store (IsoGraph and the three
## WGCNA baselines): 02_module_discovery/_m/module_sizes/{module_sizes,module_size_summary}.tsv,
## each fit stamped with its partition hash. Reports the giant modules rather than fixing them.
## Reads every fit, so run after all of 01a-01i:
##   sbatch --dependency=afterok:<01a..01i> 02_module_discovery/_h/02c.module_sizes.sh
## Calls the env interpreter directly, so the job cannot die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 02_module_discovery/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "**** module_sizes ****"
"${PY}" -u -m isograph_benchmark.real_data.module_sizes "$@"
log "**** complete ****"
