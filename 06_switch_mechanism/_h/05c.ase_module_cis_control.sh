#!/usr/bin/env bash
#SBATCH --job-name=ase-module-cis
#SBATCH --partition=RM-shared
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/ase-module-cis-%j.log

## BRIDGES-2 (06a Priority 2, item 2). The allelic test reports cis control per GENE; this
## aggregates it to the co-switching MODULE and asks whether the modules with proportionally
## more cis-controlled genes are the ones carrying the aging association. It reads
## allelic_test.parquet and refits nothing -- no model, and no module eigengene is mapped
## (module-level eigen-QTL was rejected on power).
##
##   sbatch 06_switch_mechanism/_h/05c.ase_module_cis_control.sh [--regions caudate hippocampus]
##
## Writes module_cis_control.parquet per region, plus a combined table,
## module_cis_control_summary.json and MODULE_CIS_CONTROL.md at the arm root.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 06_switch_mechanism/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge 2>/dev/null || true
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

log "**** Job starts ****"
python -m isograph_benchmark.real_data.ase_junction_allelic --stage module "$@"
log "**** Complete -> 06_switch_mechanism/_m/ase_junction_switch/ ****"
