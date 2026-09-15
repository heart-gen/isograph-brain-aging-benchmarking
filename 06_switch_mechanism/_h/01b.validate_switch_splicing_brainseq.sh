#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=validate-switch-splicing
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=12
#SBATCH --time=03:00:00
#SBATCH --array=0-5
#SBATCH --output=06_switch_mechanism/_m/logs/validate-switch-splicing-%A_%a.log

## Orthogonal validation of IsoGraph switches against junction-derived splicing
## (validate_switch_splicing.py): are genes/junctions IsoGraph calls switching also
## supported by the INDEPENDENT LIBD PSI quantification (split-read based, does not use
## the transcript quantifier that feeds IsoGraph)? The 200k-event spline sweep is too
## heavy for the login node -- run here. One (region, trait) per array task.
## Usage: sbatch 06_switch_mechanism/_h/01b.validate_switch_splicing_brainseq.sh [--alpha 0.05]
set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 06_switch_mechanism/_m/logs

# array index -> (trait, region, switch_set, analysis). age in the control cohorts
# (caudate/hippocampus/dlpfc); dx = SCZD-vs-Control in the caudate SCZD cohort. Two switch-sets:
# module = members of trait-associated co-switch modules (cross-region story; only caudate has
# age modules in BrainSEQ -- GTEx has more); marginal = per-gene switch calls (strong caudate
# corroboration, OR~100). event analysis is switch-set-independent so it rides the module tasks.
SPECS=(
    "age caudate     module   both"
    "age caudate     marginal gene"
    "age hippocampus module   both"
    "age dlpfc       module   both"
    "dx  caudate     module   both"
    "dx  caudate     marginal gene"
)
read -r TRAIT REGION SWITCHSET ANALYSIS <<< "${SPECS[${SLURM_ARRAY_TASK_ID:-0}]}"

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** validate switch splicing: trait=${TRAIT} region=${REGION} set=${SWITCHSET} analysis=${ANALYSIS} ****"
python -m isograph_benchmark.real_data.validate_switch_splicing \
    --trait "${TRAIT}" --region "${REGION}" --switch-set "${SWITCHSET}" --analysis "${ANALYSIS}" "$@"
conda deactivate
log_message "**** Complete ****"
