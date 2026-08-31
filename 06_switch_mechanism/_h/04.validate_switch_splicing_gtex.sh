#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=gtex-validate-switch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=12
#SBATCH --time=03:00:00
#SBATCH --array=0-12
#SBATCH --output=06_switch_mechanism/_m/logs/gtex-validate-switch-%A_%a.log

## Orthogonal validation of IsoGraph switches vs GTEx within-gene junction usage
## (validate_switch_splicing.py --cohort gtex). Are members of age-associated co-switch
## modules enriched for an INDEPENDENT junction-level age-splicing signal (split-read based,
## not the RSEM transcript quantifier that feeds IsoGraph)? GTEx has 115 age-associated
## modules across 13 regions (vs BrainSEQ's 8, caudate-only), so this is the powered
## cross-region replication. Requires build_gtex_junction_usage first. One region per task.
## Usage: sbatch 06_switch_mechanism/_h/04.validate_switch_splicing_gtex.sh [--switch-set marginal]
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

REGIONS=(
    amygdala anterior_cingulate_cortex_ba24 caudate_basal_ganglia
    cerebellar_hemisphere cerebellum cortex frontal_cortex_ba9
    hippocampus hypothalamus nucleus_accumbens_basal_ganglia
    putamen_basal_ganglia spinal_cord_cervical_c_1 substantia_nigra
)
REGION="${REGIONS[${SLURM_ARRAY_TASK_ID:-0}]}"

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** GTEx switch validation: ${REGION} (age, module) ****"
python -m isograph_benchmark.real_data.validate_switch_splicing \
    --cohort gtex --trait age --region "${REGION}" --analysis gene --switch-set module "$@"
conda deactivate
log_message "**** Complete ****"
