#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-signal-prep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%j.log
## Signal-level coloc, prep: freeze the target grid from the abf arm and write the GTEx credible
## sets + work_list.tsv that 06b.coloc_signal_susie.sh indexes. Reads the targets that
## 04d.coloc_modality_prep.sh wrote for the same arm, so both estimators share one grid; the
## assembly step is 07b.coloc_signal_susie_meta.sh.
##   sbatch 05_genetic_anchoring/_h/05b.coloc_signal_susie_prep.sh            # switch arm
##   sbatch 05_genetic_anchoring/_h/05b.coloc_signal_susie_prep.sh --arm background
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

log "**** coloc_signal_susie --stage prep $* ****"
"${PY}" -u -m isograph_benchmark.real_data.coloc_signal_susie --stage prep "$@"
log "**** complete ****"
