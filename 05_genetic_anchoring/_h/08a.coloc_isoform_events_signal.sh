#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-events-signal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%j.log

## Resolved isoform events for the signal-level (coloc.susie) all-introns nominations, written
## under coloc_signal_susie/all_introns/coloc_isoform_events.parquet. The CLPP-layer events and the
## signed direction are 04b.coloc_direction.sh; direction is not re-run here.
##
## Run after 07b.coloc_signal_susie_meta.sh --stage meta --sqtl all (the nominations) and 04b:
##   sbatch 05_genetic_anchoring/_h/08a.coloc_isoform_events_signal.sh
## Stage 06 reads the table (06_switch_mechanism/_h/02b.switch_orthogonal_confirm.sh --events signal).
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

log "**** resolve isoform events, signal layer (all-introns nominations) ****"
"${PY}" -u -m isograph_benchmark.real_data.coloc_isoform_events --layer signal "$@"
log "**** complete ****"
