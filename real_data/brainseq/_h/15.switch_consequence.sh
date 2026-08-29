#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=switch-consequence
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --array=0-16
#SBATCH --output=real_data/brainseq/_m/logs/switch-consequence-%A_%a.log

## Switch coding-consequence enrichment (switch_consequence.py) across all 17 switch-layer
## regions: does the IsoGraph switch axis preferentially select coding/UTR-consequential
## isoform pairs vs a within-gene random-pair null? Deterministic (--seed 13). One region per
## array task. Also emits a gene-level block-bootstrap SE / CI on the log enrichment
## (--n-boot, default 2000), which the cross-region random-effects meta pools.
## Usage: sbatch real_data/brainseq/_h/15.switch_consequence.sh [--n-perm N] [--n-boot N]
set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/brainseq/_m/logs

# array index -> (tree, region): SCZD + 3 BrainSeq aging + 13 GTEx aging. Region names
# collide across trees (hippocampus), so the artifact dir + label are tree-qualified.
SPECS=(
    "brainseq caudate_sczd"
    "brainseq caudate"
    "brainseq hippocampus"
    "brainseq dlpfc"
    "gtex amygdala"
    "gtex anterior_cingulate_cortex_ba24"
    "gtex caudate_basal_ganglia"
    "gtex cerebellar_hemisphere"
    "gtex cerebellum"
    "gtex cortex"
    "gtex frontal_cortex_ba9"
    "gtex hippocampus"
    "gtex hypothalamus"
    "gtex nucleus_accumbens_basal_ganglia"
    "gtex putamen_basal_ganglia"
    "gtex spinal_cord_cervical_c_1"
    "gtex substantia_nigra"
)
read -r TREE REGION <<< "${SPECS[${SLURM_ARRAY_TASK_ID:-0}]}"
ARTIFACT="real_data/${TREE}/${REGION}/_m/isograph_vae"

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** switch consequence: ${TREE}/${REGION} ****"
python -m isograph_benchmark.real_data.switch_consequence \
    --region "${TREE}_${REGION}" --artifact-dir "${ARTIFACT}" "$@"
conda deactivate
log_message "**** Complete ****"
