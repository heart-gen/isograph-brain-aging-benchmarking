#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=nova-family-renomination
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=08:00:00
#SBATCH --output=02_module_discovery/brainseq/_m/logs/nova-family-renomination-%j.log

# Rebuild NOVA-family nominations from the complete switch-pair universe with
# two-sided intronic-opportunity eligibility and adjusted module enrichment.
#
# Submit from the repository root:
#   sbatch 02_module_discovery/brainseq/_h/25.nova_family_renomination.sh

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
if [[ ! -f "${CONFIG}" || ! -x "${PYTHON}" ]]; then
    echo "ERROR: missing configuration or designated Python interpreter."
    exit 1
fi

export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p 02_module_discovery/brainseq/_m/logs

log_message "**** NOVA-family opportunity-controlled re-nomination ****"
sha256sum "${CONFIG}"
"${PYTHON}" -m isograph_benchmark.real_data.nova_family_renomination \
    --config "${CONFIG}" "$@"
log_message "outputs: real_data/_m/neuronal_clip/nova_family_renomination"
log_message "candidate freeze: reports/neuronal_clip/nova_family_candidate_manifest.parquet"
log_message "**** complete ****"
