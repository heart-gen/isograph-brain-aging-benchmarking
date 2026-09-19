#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-signal-meta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%j.log
## Signal-level coloc, assembly: turn the shards 06b.coloc_signal_susie.sh wrote into the cell
## tables, the estimator hierarchy and COLOC_SIGNAL_SUSIE.md. One call writes both the
## gtex_matched primary and the all_signals arm.
##
## The flags must match the env the shards were produced under (COLOC_SIGNAL_SQTL,
## COLOC_GWAS_MAX_SNPS); they pick the directory, so a mismatch reads an empty tree:
##   sbatch 05_genetic_anchoring/_h/07b.coloc_signal_susie_meta.sh                                # representative
##   sbatch 05_genetic_anchoring/_h/07b.coloc_signal_susie_meta.sh --sqtl all                     # all introns (primary nominations)
##   sbatch 05_genetic_anchoring/_h/07b.coloc_signal_susie_meta.sh --max-snps 30000               # SNP-guard sensitivity
##   sbatch 05_genetic_anchoring/_h/07b.coloc_signal_susie_meta.sh --sqtl all --max-snps 30000
##
## Memory on PSC is --cpus-per-task x 2000MB; 8 cpus = 16 GB. Do NOT pass --mem. Calls the env
## interpreter directly, so the job cannot die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 05_genetic_anchoring/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "**** coloc_signal_susie --stage meta ${*:-(representative)} ****"
"${PY}" -u -m isograph_benchmark.real_data.coloc_signal_susie --stage meta "$@"
log "**** complete ****"
