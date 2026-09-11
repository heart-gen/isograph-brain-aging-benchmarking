#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=switch-orthogonal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/switch-orthogonal-%j.log

## Orthogonal ONT long-read confirmation of the genetically anchored switch pairs.
##   --mode anchored     (default) splicing-led coloc genes vs abundance-matched background
##   --mode global-null  matched null for the overall long-read switch-like rate
## Extra arguments are forwarded, e.g.
##   sbatch 06_switch_mechanism/_h/05.switch_orthogonal_confirm.sh --mode global-null
##
## --events picks the colocalization layer the anchored pairs come from:
##   clpp    (default) eCAVIAR events; writes _m/switch_orthogonal_confirm/ (S-real-8, S14)
##   signal  coloc.susie all-introns nominations; writes _m/switch_orthogonal_confirm/signal_coloc/
##           Cross-tissue exceptions (UNC13A; coloc_isoform_events.CROSS_TISSUE_EXCEPTIONS) are
##           scored by the same code into exception_pair_confirmation.parquet, kept out of the
##           set-level comparison and its background -- they change the interpretation.
## The signal layer needs its event table first (login-node safe):
##   python -m isograph_benchmark.real_data.coloc_isoform_events --layer signal
##   sbatch 06_switch_mechanism/_h/05.switch_orthogonal_confirm.sh --events signal
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 06_switch_mechanism/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

## The Bambu matrix is git-LFS tracked; a fresh clone leaves a pointer behind.
COUNTS=inputs/_m/longread_aged_dlpfc/counts_transcript.txt
if [[ ! -s "${COUNTS}" ]] || head -c 40 "${COUNTS}" | grep -q 'git-lfs'; then
    log "restoring ${COUNTS} from git-lfs"
    git lfs pull --include="${COUNTS}" || \
        "${PY}" -m isograph_benchmark.real_data.longread_switch_confirm fetch
fi

log "switch_orthogonal_confirm $*"
"${PY}" -m isograph_benchmark.real_data.switch_orthogonal_confirm "$@"
log "**** complete ****"
