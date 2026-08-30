#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=mtrust-rep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --array=0-5
#SBATCH --output=03_module_trust/_m/logs/mtrust-rep-%A_%a.log

## Q3 cross-cohort aging replication (module_trust `replication` subcommand).
##
## This driver was referenced by MODULE_TRUST_SUMMARY.md but never committed — the arm that
## produces the "N of M modules replicate" per-pair tables had no reproducible launcher.
## Array covers {isograph, wgcna} x {caudate, hippocampus, dlpfc_ba9}.
##
## Writes 03_module_trust/_m/stability/module_trust/module_aging_replication__{pair}__{method}.parquet
## Requires Q1 `stability` to have run for both cohorts of each pair.
##
## Usage: sbatch 03_module_trust/_h/09.module_trust_replication.sh
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 03_module_trust/_m/logs

METHODS=(isograph wgcna)
PAIRS=(caudate hippocampus dlpfc_ba9)
IDX="${SLURM_ARRAY_TASK_ID:-0}"
METHOD="${METHODS[$((IDX / 3))]}"
PAIR="${PAIRS[$((IDX % 3))]}"

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log "**** Q3 replication: method=${METHOD} pair=${PAIR} ****"
python -m isograph_benchmark.real_data.module_trust replication \
  --pair "${PAIR}" --method "${METHOD}" "$@"

conda deactivate
log "**** done (${METHOD}/${PAIR}) ****"
