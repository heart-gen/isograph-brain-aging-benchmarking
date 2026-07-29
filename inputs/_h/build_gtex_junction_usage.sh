#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=gtex-junction-usage
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=10
#SBATCH --time=02:00:00
#SBATCH --array=0-12
#SBATCH --output=real_data/gtex/_m/logs/gtex-junction-usage-%A_%a.log

## Ingest GTEx v11 STAR junction counts -> per-region within-gene junction usage
## (build_gtex_junction_usage.py). Streams the 523,817 x 19,788 junction GCT once per
## region with usecols limited to that region's bundle samples, restricts to switch-scored
## genes, and writes inputs/processed/gtex/<region>/junction_usage.parquet -- the GTEx
## analogue of BrainSEQ PSI for the orthogonal switch validation. One region per array task.
## Usage: sbatch inputs/_h/build_gtex_junction_usage.sh
set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/gtex/_m/logs

REGIONS=(
    amygdala anterior_cingulate_cortex_ba24 caudate_basal_ganglia
    cerebellar_hemisphere cerebellum cortex frontal_cortex_ba9
    hippocampus hypothalamus nucleus_accumbens_basal_ganglia
    putamen_basal_ganglia spinal_cord_cervical_c_1 substantia_nigra
)
REGION="${REGIONS[${SLURM_ARRAY_TASK_ID:-0}]}"

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** GTEx junction usage: ${REGION} ****"
python -m isograph_benchmark.inputs.build_gtex_junction_usage --region "${REGION}" "$@"
conda deactivate
log_message "**** Complete ****"
