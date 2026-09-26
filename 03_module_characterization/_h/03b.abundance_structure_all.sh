#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=abundance-structure-all
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-17
#SBATCH --time=01:00:00
#SBATCH --output=03_module_characterization/_m/logs/abundance-structure-all-%A_%a.log

# Abundance-vs-switch axis orthogonality for EVERY analysis, one array task per store.
# The separability statement in the manuscript ("the two channels capture partly distinct
# within-gene variation") must hold across the 16 aging analyses too, not only the BrainSEQ
# SCZD store that _h/03a happens to run for the example gene. Writes
# 02_module_discovery/<cohort>/<region>/_m/isograph_vae/abundance_structure/axis_orthogonality.parquet
# per task; pool them with _h/03c (local), which also rebuilds figSeparation.
#
#   task 1     -> brainseq-sczd (caudate_sczd)
#   tasks 2-4  -> brainseq-aging caudate, hippocampus, dlpfc
#   tasks 5-17 -> gtex-aging, the 13 regions in run_models.GTEX_REGIONS order
#
# Requires feature_scores.parquet (stage 02) and the input bundle. Flags forwarded after
# the script name, e.g. --variant with-abundance.

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the isograph-brain-aging-benchmarking repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 03_module_characterization/_m/logs

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

TARGETS=(
    "brainseq-sczd"
    "brainseq-aging --region caudate"
    "brainseq-aging --region hippocampus"
    "brainseq-aging --region dlpfc"
    "gtex-aging --region amygdala"
    "gtex-aging --region anterior_cingulate_cortex_ba24"
    "gtex-aging --region caudate_basal_ganglia"
    "gtex-aging --region cerebellar_hemisphere"
    "gtex-aging --region cerebellum"
    "gtex-aging --region cortex"
    "gtex-aging --region frontal_cortex_ba9"
    "gtex-aging --region hippocampus"
    "gtex-aging --region hypothalamus"
    "gtex-aging --region nucleus_accumbens_basal_ganglia"
    "gtex-aging --region putamen_basal_ganglia"
    "gtex-aging --region spinal_cord_cervical_c_1"
    "gtex-aging --region substantia_nigra"
)
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
if (( TASK_ID < 1 || TASK_ID > ${#TARGETS[@]} )); then
    echo "ERROR: unknown array task id ${TASK_ID} (expected 1-${#TARGETS[@]})"; exit 1
fi
# shellcheck disable=SC2206  # intentional word split of the target spec
TARGET=(${TARGETS[$((TASK_ID - 1))]})
log_message "**** Axis orthogonality: ${TARGET[*]} ****"

python -u -m isograph_benchmark.real_data.abundance_structure_separation \
    "${TARGET[@]}" --orthogonality-only "$@"

conda deactivate
log_message "**** Complete ****"
