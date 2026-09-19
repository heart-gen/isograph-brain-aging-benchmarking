#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=abundance-structure-sep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=03_module_characterization/_m/logs/abundance-structure-%j.log

# Abundance-vs-isoform-structure separation inputs for the figSeparation panel.
# Derives the per-gene abundance/switch axis orthogonality, pools the de-confounded
# incremental-test category counts across cohorts/regions, and extracts one
# composition-unique example gene (stable total abundance, real isoform switch),
# writing to
# 02_module_discovery/brainseq/caudate_sczd/_m/isograph_vae/abundance_structure/.
#
# Requires feature_scores.parquet + the incremental_association/ + composition_unique/
# outputs (07/08 runs). Flags forwarded after the script name, e.g. --gene MAP1A.

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 03_module_characterization/_m/logs
log_message "**** Abundance-structure separation ****"

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

python -m isograph_benchmark.real_data.abundance_structure_separation brainseq-sczd "$@"

conda deactivate
log_message "**** Complete ****"
