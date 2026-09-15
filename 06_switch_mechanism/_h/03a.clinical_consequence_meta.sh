#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=clinical-consequence-meta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=06_switch_mechanism/_m/logs/%x-%j.log

## Cross-region clinical-consequence rollup: clinical_consequence_meta.parquet +
## CLINICAL_CONSEQUENCE_META.md from the 17 per-region runs of 02c.clinical_consequence.sh.
## Previously a comment at the bottom of that wrapper and run by hand. Run after the array:
##   sbatch --dependency=afterok:<02c> 06_switch_mechanism/_h/03a.clinical_consequence_meta.sh
## Calls the env interpreter directly, so the job cannot die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 06_switch_mechanism/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "**** clinical_consequence_meta ****"
"${PY}" -u -m isograph_benchmark.real_data.clinical_consequence_meta "$@"
log "**** complete ****"
