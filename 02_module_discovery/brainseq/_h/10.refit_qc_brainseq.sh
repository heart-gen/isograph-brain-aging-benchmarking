#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=refit-qc-bs
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=24   # 24 x 2000M = 48G (dlpfc fit needs 48G; caudate/hippo fit comfortably)
#SBATCH --time=04:00:00
#SBATCH --array=1-3
#SBATCH --output=02_module_discovery/brainseq/_m/logs/%x-%A_%a.log

# Materialize the DE-aligned QC-adjusted aging associations across the 3 BrainSEQ aging
# regions. Per region: (1) re-fit the VAE (run_models now emits age_spline_qc_adjusted.parquet
# via _rnaseqc_covariate_table + feature_reconstruction.parquet), then (2) re-project the
# 3 IsoGraph tiers (project_tiers now emits per-tier age_spline_qc_adjusted.parquet).
# Heavy compute -> SLURM only.
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
mkdir -p 02_module_discovery/brainseq/_m/logs

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

case "${SLURM_ARRAY_TASK_ID}" in
    1) REGION="caudate" ;;
    2) REGION="hippocampus" ;;
    3) REGION="dlpfc" ;;
    *) echo "ERROR: bad array idx ${SLURM_ARRAY_TASK_ID}"; exit 1 ;;
esac

log_message "**** [brainseq-aging/${REGION}] step 1: re-fit VAE (emits age_spline_qc_adjusted + reconstruction) ****"
python -m isograph_benchmark.real_data.run_models brainseq-aging --region "${REGION}"

log_message "**** [brainseq-aging/${REGION}] step 2: re-project tiers (per-tier QC-adjusted spline) ****"
python -m isograph_benchmark.real_data.project_tiers brainseq-aging --region "${REGION}"

conda deactivate
log_message "**** [brainseq-aging/${REGION}] QC re-fit complete ****"
