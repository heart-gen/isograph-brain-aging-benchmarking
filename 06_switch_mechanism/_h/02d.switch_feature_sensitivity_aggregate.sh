#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=switch-sensitivity-aggregate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/%x-%j.log

## Cross-region feature-sensitivity report (switch_feature_sensitivity --aggregate): collects the
## per-region runs of 01h (BrainSEQ caudate) and 01i (the 13 GTEx regions) into
## SWITCH_FEATURE_SENSITIVITY.md and the region worst-case table, and runs the quantification-
## pipeline axis once. Run after both:
##   sbatch --dependency=afterok:<01h>:<01i> 06_switch_mechanism/_h/02d.switch_feature_sensitivity_aggregate.sh
## Bridges memory is --cpus-per-task x 2000MB; do NOT pass --mem. Calls the env interpreter directly.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 06_switch_mechanism/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
log "switch_feature_sensitivity --aggregate"
"${PY}" -m isograph_benchmark.real_data.switch_feature_sensitivity --aggregate "$@"
log "**** complete ****"
