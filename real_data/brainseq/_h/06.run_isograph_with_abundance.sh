#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-iso-withab
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=32
#SBATCH --array=1-4
#SBATCH --time=03:00:00
#SBATCH --output=real_data/brainseq/_m/logs/%x-%A_%a.log
# Part 2 fan-out: re-run IsoGraph with abundance-abundance edges ENABLED
# (alpha_abundance_grid calibration) for every BrainSEQ dataset, not just
# caudate aging. Each region uses its data-driven best Leiden resolution
# (BEST_LEIDEN_RESOLUTION in run_models.py, selected by the Part 1 sweep).
# Artifacts go to <region>/_m/isograph_vae_with_abundance/ so the baseline
# isograph_vae/ results are preserved for side-by-side comparison.
#
# Array tasks:
#   1 -> brainseq-aging caudate      (resolution 2.25)
#   2 -> brainseq-aging hippocampus  (resolution 2.0)
#   3 -> brainseq-aging dlpfc        (resolution 3.0)
#   4 -> brainseq-sczd  caudate_sczd (resolution 2.0)
#
# Requires the standard runs (01/02) to have produced the input bundles; this
# script refits the VAE, it does not reuse saved edges.
set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "ERROR: project root does not exist: ${PROJECT_ROOT}"
    exit 1
fi
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the isograph-brain-aging-benchmarking repo root or set ISOGRAPH_BENCHMARK_ROOT."
    echo "Current project root candidate: ${PROJECT_ROOT}"
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/brainseq/_m/logs

log_message "**** BrainSEQ with-abundance IsoGraph job starts ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID:-local}"
echo "Array task: ${SLURM_ARRAY_TASK_ID:-none}"
echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "Activating IsoGraph environment"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "Checking Python analysis dependencies"
python -c "import isograph_benchmark, numpy, pandas, scipy, patsy"

# Limit BLAS/OpenMP threads to bound per-thread SVD workspace (17k+ genes).
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
case "${TASK_ID}" in
    1) log_message "with-abundance: brainseq-aging caudate"
       python -m isograph_benchmark.real_data.run_models brainseq-aging \
           --variant with-abundance --region caudate ;;
    2) log_message "with-abundance: brainseq-aging hippocampus"
       python -m isograph_benchmark.real_data.run_models brainseq-aging \
           --variant with-abundance --region hippocampus ;;
    3) log_message "with-abundance: brainseq-aging dlpfc"
       python -m isograph_benchmark.real_data.run_models brainseq-aging \
           --variant with-abundance --region dlpfc ;;
    4) log_message "with-abundance: brainseq-sczd caudate_sczd"
       python -m isograph_benchmark.real_data.run_models brainseq-sczd \
           --variant with-abundance ;;
    *) echo "ERROR: unknown array task id ${TASK_ID} (expected 1-4)"; exit 1 ;;
esac

conda deactivate
log_message "**** BrainSEQ with-abundance IsoGraph job ends ****"
