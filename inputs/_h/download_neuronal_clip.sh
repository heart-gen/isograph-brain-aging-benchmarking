#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=dl-neuronal-clip
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=02:00:00
#SBATCH --output=inputs/_m/download-neuronal-clip-%j.log

# Reproducibly acquire the public Tier-1 neuronal CLIP inputs.
#
# The configuration pins source URLs, expected byte counts, and SHA-256 values.
# The downloader freezes the motif-derived candidate manifest before reading peaks,
# resumes interrupted transfers through atomic `.part` files, verifies archives,
# safely extracts configured archives, and records controlled-access EGA datasets
# without attempting to download them.
#
# Submit from the repository root:
#   sbatch inputs/_h/download_neuronal_clip.sh
#
# Re-running is safe: existing files and extracted members are checksum-verified.

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
mkdir -p inputs/_m

log_message "**** neuronal CLIP acquisition ****"
log_message "project: ${PROJECT_ROOT}"
log_message "config: ${CONFIG}"
sha256sum "${CONFIG}"

"${PYTHON}" -m isograph_benchmark.real_data.neuronal_clip_fetch \
    --config "${CONFIG}" fetch --extract "$@"

log_message "dataset manifest: 07_rbp_regulation/_m/neuronal_clip_manifests/dataset_manifest.tsv"
log_message "candidate manifest: 07_rbp_regulation/_m/neuronal_clip_manifests/candidate_manifest.parquet"
log_message "context QC: 07_rbp_regulation/_m/neuronal_clip_manifests/context_qc.tsv"
log_message "**** complete ****"
