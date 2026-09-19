#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=gtex-incremental-assoc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-13
#SBATCH --time=03:00:00
#SBATCH --output=03_module_characterization/_m/logs/incremental-assoc-%A_%a.log

# Switch-vs-abundance incremental association on the 13 GTEx v11 brain regions.
# Runs the gene-level DE-CONFOUNDED test (switch|abundance and abundance|switch,
# df=3 AGE spline) + the module-level incremental contrast on the saved IsoGraph
# features, writing to 02_module_discovery/gtex/<region>/_m/isograph_vae/incremental_association/.
#
# Requires feature_scores.parquet + modules.parquet from 02_module_discovery/_h/01c.run_isograph_gtex.sh.
# Flags forwarded after the script name, e.g. --variant with-abundance.

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "ERROR: project root does not exist: ${PROJECT_ROOT}"
    exit 1
fi
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the isograph-brain-aging-benchmarking repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 03_module_characterization/_m/logs

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

REGIONS=(
    amygdala anterior_cingulate_cortex_ba24 caudate_basal_ganglia
    cerebellar_hemisphere cerebellum cortex frontal_cortex_ba9
    hippocampus hypothalamus nucleus_accumbens_basal_ganglia
    putamen_basal_ganglia spinal_cord_cervical_c_1 substantia_nigra
)
REGION="${REGIONS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}"
log_message "**** GTEx incremental association: ${REGION} ****"

python -m isograph_benchmark.real_data.incremental_association gtex-aging --region "${REGION}" "$@"

conda deactivate
log_message "**** Complete ****"
