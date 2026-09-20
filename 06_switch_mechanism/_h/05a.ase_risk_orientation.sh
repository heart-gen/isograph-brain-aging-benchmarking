#!/usr/bin/env bash
#SBATCH --job-name=ase-risk-orient
#SBATCH --partition=RM-shared
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/ase-risk-orient-%j.log

## BRIDGES-2 (PI item 10a, step 5). Step 4 (Quest) fits beta against the switch-QTL lead's ALT
## allele, which carries no disease meaning. This re-signs it to the GWAS risk allele of the
## locus each gene colocalizes in, using signed LD between the two variants in the BrainSEQ
## genotypes the QTL mapping used. It reads the merged allelic tables and refits nothing, so
## it runs here, where the GWAS sumstats and the genotype panel are.
##
##   sbatch 06_switch_mechanism/_h/05a.ase_risk_orientation.sh [--region caudate]
##
## Writes risk_orientation.parquet + ASE_RISK_ORIENTATION.md beside the allelic outputs.
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
python -m isograph_benchmark.real_data.ase_risk_orientation "$@"
log "**** Complete -> 06_switch_mechanism/_m/ase_junction_switch/<region>/ ****"
