#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-mod-compare
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%j.log
## Per-gene sQTL-vs-eQTL colocalization contrast, step 4: the four gene-pool arms side by side
## (switch, background, wgcna_switch, wgcna_multiplex) -> arm_comparison.parquet. Reads each arm's
## contrast.parquet, so run it after every 06a.coloc_modality_meta.sh:
##   sbatch --dependency=afterok:<06a x4> 05_genetic_anchoring/_h/07a.coloc_modality_compare.sh
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

log "**** coloc modality contrast: compare arms ****"
"${PY}" -u -m isograph_benchmark.real_data.coloc_modality_contrast --stage compare "$@"
log "**** complete ****"
