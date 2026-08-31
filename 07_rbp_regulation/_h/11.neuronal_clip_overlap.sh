#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=neuronal-clip-overlap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=08:00:00
#SBATCH --output=07_rbp_regulation/_m/logs/neuronal-clip-overlap-%j.log

# Build assay-callable window calls and matched effects for eligible public
# human TDP-43 and PTBP2 contexts. Run stage 23 successfully first.
#
# Submit from the repository root:
#   sbatch 07_rbp_regulation/_h/11.neuronal_clip_overlap.sh

set -euo pipefail
umask 0027

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repository root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi

CONFIG="${NEURONAL_CLIP_CONFIG:-configs/neuronal_clip.yaml}"
PYTHON="/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python"
ELIGIBILITY="07_rbp_regulation/_m/neuronal_clip/motif_opportunity/candidate_eligibility.parquet"
if [[ ! -f "${CONFIG}" || ! -x "${PYTHON}" ]]; then
    echo "ERROR: missing configuration or designated Python interpreter."
    exit 1
fi
if [[ ! -f "${ELIGIBILITY}" ]]; then
    echo "ERROR: missing ${ELIGIBILITY}; complete stage 23 first."
    exit 1
fi

export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p 07_rbp_regulation/_m/logs

log_message "**** neuronal CLIP callable-window overlap ****"
sha256sum "${CONFIG}" "${ELIGIBILITY}"
"${PYTHON}" -m isograph_benchmark.real_data.neuronal_clip_overlap \
    --config "${CONFIG}" "$@"
log_message "outputs: 07_rbp_regulation/_m/neuronal_clip/overlap"
log_message "**** complete ****"
