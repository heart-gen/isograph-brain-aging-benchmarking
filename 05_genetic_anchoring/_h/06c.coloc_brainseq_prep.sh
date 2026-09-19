#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-brainseq-prep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%j.log
## BrainSEQ signal-level coloc, prep: targets and the work_list.tsv that
## 07c.coloc_brainseq_susie.sh indexes. EA-only by construction -- the CLI refuses any other arm,
## and refuses ea_only unless both 02g.brainseq_qtl_checks.sh checks passed for it.
##
## Run after 03c.coloc_gwas_susie.sh (GWAS SuSiE cache), 05b.coloc_signal_susie_prep.sh (target
## grid) and 02g.brainseq_qtl_checks.sh with SWQTL_ARM=ea_only:
##   sbatch 05_genetic_anchoring/_h/06c.coloc_brainseq_prep.sh
## The assembly step is 08c.coloc_brainseq_meta.sh.
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

log "**** coloc_brainseq --stage prep $* ****"
"${PY}" -u -m isograph_benchmark.real_data.coloc_brainseq --stage prep "$@"
log "**** complete ****"
