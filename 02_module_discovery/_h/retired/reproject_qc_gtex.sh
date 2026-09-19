#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=reproj-qc-gtex
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16   # 16 x 2000M = 32G; projection only (no VAE re-fit)
#SBATCH --time=02:00:00
#SBATCH --array=1-13
#SBATCH --output=02_module_discovery/_m/logs/%x-%A_%a.log

# Re-project the 3 IsoGraph tiers ONLY (no VAE re-fit) to pick up the collision-safe
# GTEx QC merge fix in project_tiers._associate. The prior fan-out (02_module_discovery/_h/retired/refit_qc_gtex.sh)
# left the GTEx tier age_spline_qc_adjusted.parquet byte-identical to baseline because
# the native SMEXNCRT/SM3PB75P columns collided in the merge and silently dropped out.
# Fits + feature_reconstruction.parquet already exist; this only recomputes projections.
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
mkdir -p 02_module_discovery/_m/logs

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

REGIONS=(
    amygdala
    anterior_cingulate_cortex_ba24
    caudate_basal_ganglia
    cerebellar_hemisphere
    cerebellum
    cortex
    frontal_cortex_ba9
    hippocampus
    hypothalamus
    nucleus_accumbens_basal_ganglia
    putamen_basal_ganglia
    spinal_cord_cervical_c_1
    substantia_nigra
)
REGION="${REGIONS[$((SLURM_ARRAY_TASK_ID - 1))]}"
if [[ -z "${REGION:-}" ]]; then
    echo "ERROR: bad array idx ${SLURM_ARRAY_TASK_ID}"; exit 1
fi

log_message "**** [gtex-aging/${REGION}] re-project tiers (collision-safe QC merge) ****"
python -m isograph_benchmark.real_data.project_tiers gtex-aging --region "${REGION}"

conda deactivate
log_message "**** [gtex-aging/${REGION}] QC re-projection complete ****"
