#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=lr-validation
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=24  # 24 x 2000M = 48G; a full-data dlpfc fit needs ~48G
#SBATCH --time=04:00:00
#SBATCH --array=1-7
#SBATCH --output=real_data/stability/_m/logs/%x-%A_%a.log

# Validation gate B.2 — single-LR/optimizer config. One FULL-data IsoGraph fit per region
# at a single fixed learning rate (NO per-region tuning) with gradient clipping, recording
# reconstruction RMSE and whether training diverged. Tests whether the merged grad_clip_norm
# + divergence guard remove the need for the hand-tuned GTEx lr=3e-4. The 7th region is the
# GTEx region that diverged at lr=1e-3 (nucleus_accumbens) — the decisive case. Heavy
# compute -> SLURM only.

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/stability/_m/logs

# array index -> (cohort, region): 6 trust-funnel regions + the diverging GTEx region.
SPECS=(
    "brainseq caudate"
    "brainseq hippocampus"
    "brainseq dlpfc"
    "gtex caudate_basal_ganglia"
    "gtex hippocampus"
    "gtex frontal_cortex_ba9"
    "gtex nucleus_accumbens_basal_ganglia"
)
spec="${SPECS[$((SLURM_ARRAY_TASK_ID - 1))]}"
read -r COHORT REGION <<< "${spec}"

# One LR for ALL regions (default 1e-3 = the BrainSEQ default that diverges on some GTEx
# without clipping); grad clip on. Override via:
#   sbatch --export=ALL,STABILITY_LR=5e-4,STABILITY_GRAD_CLIP=1.0 ...
LR="${STABILITY_LR:-1e-3}"
GRAD_CLIP="${STABILITY_GRAD_CLIP:-1.0}"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** single-LR fit: ${COHORT}/${REGION} (lr=${LR}, grad_clip=${GRAD_CLIP}) starts ****"
python -m isograph_benchmark.real_data.stability fit-rmse \
    --cohort "${COHORT}" --region "${REGION}" --lr "${LR}" --grad-clip-norm "${GRAD_CLIP}"
conda deactivate
log_message "**** Complete ****"
