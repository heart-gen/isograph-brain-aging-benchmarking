#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=tiers-fanout
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=24   # 24 x 2000M = 48G (dlpfc fit needs 48G; others fit comfortably)
#SBATCH --time=04:00:00
#SBATCH --array=1-5
#SBATCH --output=02_module_discovery/brainseq/_m/logs/%x-%A_%a.log

# Fan out the validated caudate pilot to the remaining 5 trust-funnel aging regions.
# Per region: (1) re-fit the VAE ONCE emitting feature_reconstruction.parquet, then
# (2) project that fit into the 3 IsoGraph tiers, then (3) run the 6 tier checks.
# Caudate (the pilot) is already done. Heavy compute -> SLURM only.
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

# array idx -> "analysis region"
case "${SLURM_ARRAY_TASK_ID}" in
    1) ANALYSIS="brainseq-aging" REGION="hippocampus" ;;
    2) ANALYSIS="brainseq-aging" REGION="dlpfc" ;;
    3) ANALYSIS="gtex-aging"     REGION="caudate_basal_ganglia" ;;
    4) ANALYSIS="gtex-aging"     REGION="hippocampus" ;;
    5) ANALYSIS="gtex-aging"     REGION="frontal_cortex_ba9" ;;
    *) echo "ERROR: bad array idx ${SLURM_ARRAY_TASK_ID}"; exit 1 ;;
esac

log_message "**** [${ANALYSIS}/${REGION}] step 1: re-fit VAE (emits feature_reconstruction) ****"
python -m isograph_benchmark.real_data.run_models "${ANALYSIS}" --region "${REGION}"

log_message "**** [${ANALYSIS}/${REGION}] step 2: project into 3 IsoGraph tiers ****"
python -m isograph_benchmark.real_data.project_tiers "${ANALYSIS}" --region "${REGION}"

log_message "**** [${ANALYSIS}/${REGION}] step 3: six tier checks ****"
python -m isograph_benchmark.real_data.tier_checks "${ANALYSIS}" --region "${REGION}"

conda deactivate
log_message "**** [${ANALYSIS}/${REGION}] fan-out complete ****"
