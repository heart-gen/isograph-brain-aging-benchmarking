#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=neuronal-clip-freeze
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=04:00:00
#SBATCH --output=07_rbp_regulation/_m/logs/neuronal-clip-freeze-%j.log

# Freeze the module-RBP candidate set that gates the neuronal CLIP suite
# (stages 22-27) and refresh the dataset / context QC manifests.
#
# Guards: the module-RBP nomination COUNT and the candidate IDENTITY hash.
# The identity hash is the load-bearing one -- a count guard alone passes when
# the nomination set changes membership without changing size.
#
# Submit from the repository root:
#   sbatch 07_rbp_regulation/_h/04b.neuronal_clip_freeze.sh

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
if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: missing configuration: ${CONFIG}"
    exit 1
fi
if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: missing designated Python interpreter: ${PYTHON}"
    exit 1
fi

export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p 07_rbp_regulation/_m/logs

log_message "**** neuronal CLIP candidate freeze ****"
log_message "project: ${PROJECT_ROOT}"
log_message "config: ${CONFIG}"
sha256sum "${CONFIG}"

"${PYTHON}" -m isograph_benchmark.real_data.neuronal_clip_fetch \
    --config "${CONFIG}" manifest "$@"

log_message "outputs: 07_rbp_regulation/_m/neuronal_clip_manifests/"
log_message "**** complete ****"
