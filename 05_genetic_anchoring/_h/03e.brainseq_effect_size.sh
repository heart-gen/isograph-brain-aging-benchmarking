#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-effect-size
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16
#SBATCH --array=1-2
#SBATCH --time=03:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%A_%a.log

## Switch vs abundance QTL effect sizes without selection asymmetry (brainseq_switch_qtl
## --stage effect_size). The meta stage compares |slope| only on doubly significant genes,
## which selects the noisier switch axis on larger effects. This stage uses every gene tested
## on both axes: own-lead |z|/|slope|, each axis at the other's lead, and both at a common
## variant, overall and within power deciles. Reads the paired tables from --stage meta and the
## nominal cis_qtl_pairs files from --stage map; recomputes nothing upstream.
##
##   task 1 -> all_samples   task 2 -> ea_only
##
## Bridges memory is --cpus-per-task x 2000MB (32G); do NOT pass --mem. Calls the env
## interpreter directly, so an array task cannot die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 05_genetic_anchoring/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
ARMS=(all_samples ea_only)
ARM="${ARMS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}"

log "brainseq_switch_qtl --stage effect_size --arm ${ARM}"
"${PY}" -u -m isograph_benchmark.real_data.brainseq_switch_qtl --stage effect_size --arm "${ARM}" "$@"
log "**** complete ****"
