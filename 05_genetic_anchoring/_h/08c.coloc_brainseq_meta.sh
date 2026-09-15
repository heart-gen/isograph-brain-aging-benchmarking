#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-brainseq-meta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%j.log
## BrainSEQ signal-level coloc, assembly: cells, estimator hierarchy, S_g-vs-A_g contrast and
## BRAINSEQ_COLOC.md from the 07c.coloc_brainseq_susie.sh shards. The report cross-references the
## GTEx all-introns nominations, so it also waits on 07b.coloc_signal_susie_meta.sh --sqtl all:
##   sbatch 05_genetic_anchoring/_h/08c.coloc_brainseq_meta.sh
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

log "**** coloc_brainseq --stage meta $* ****"
"${PY}" -u -m isograph_benchmark.real_data.coloc_brainseq --stage meta "$@"
log "**** complete ****"
