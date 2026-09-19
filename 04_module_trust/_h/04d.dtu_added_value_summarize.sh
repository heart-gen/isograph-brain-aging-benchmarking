#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=dtu-added-value-summarize
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:45:00
#SBATCH --output=04_module_trust/_m/logs/%x-%j.log

## Module context and held-out DTU evidence (PI review item 12a), step 3 of 3: the post-hoc
## giant-module supplement (re-reads the 02e halves and 01a partitions), then the region summaries
## (split-direction distributions, BH across the six regions, the pre-registered decision rule) and
## DTU_ADDED_VALUE.md. Fails loudly if any region's 03i output is missing. 2 cpus = 4G.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 04_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
log "module_dtu_added_value giant-sensitivity"
"${PY}" -u -m isograph_benchmark.real_data.module_dtu_added_value giant-sensitivity
log "module_dtu_added_value summarize"
"${PY}" -u -m isograph_benchmark.real_data.module_dtu_added_value summarize
log "**** complete ****"
