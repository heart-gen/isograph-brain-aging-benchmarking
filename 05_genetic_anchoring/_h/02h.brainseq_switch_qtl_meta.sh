#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-swqtl-meta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-2
#SBATCH --time=02:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%A_%a.log
## BrainSEQ switch-QTL vs abundance-QTL comparison (brainseq_switch_qtl --stage meta): the paired
## per-region tables, modality_contrast.parquet, yield_by_power_bin.parquet and
## BRAINSEQ_SWITCH_QTL.md. 01i.brainseq_switch_qtl.sh runs phenotypes + map only, so this is the
## step that turns a mapped arm into the tables 03e.brainseq_effect_size.sh reads.
##
##   task 1 -> all_samples   task 2 -> ea_only
##
##   sbatch --dependency=afterok:<25 jobs> 05_genetic_anchoring/_h/02h.brainseq_switch_qtl_meta.sh
##   sbatch --dependency=afterok:<this>    05_genetic_anchoring/_h/03e.brainseq_effect_size.sh
##
## Bridges memory is --cpus-per-task x 2000MB (16G); do NOT pass --mem. Calls the env
## interpreter directly, so an array task cannot die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 05_genetic_anchoring/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
ARMS=(all_samples ea_only)
ARM="${ARMS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}"

log "brainseq_switch_qtl --stage meta --arm ${ARM}"
"${PY}" -u -m isograph_benchmark.real_data.brainseq_switch_qtl --stage meta --arm "${ARM}" "$@"
log "**** complete ****"
