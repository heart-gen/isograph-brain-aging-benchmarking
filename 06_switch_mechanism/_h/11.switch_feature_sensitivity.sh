#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=switch-sensitivity
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/switch-sensitivity-%j.log

## Preprocessing sensitivity of the switch representation: pseudocount, transcript
## expression filter, minor-isoform threshold, identifiability, quantification pipeline.
##   sbatch 06_switch_mechanism/_h/11.switch_feature_sensitivity.sh --cohort brainseq --region caudate
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 06_switch_mechanism/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
log "switch_feature_sensitivity $*"
"${PY}" -m isograph_benchmark.real_data.switch_feature_sensitivity "$@"
log "**** complete ****"
