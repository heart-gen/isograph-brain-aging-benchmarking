#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=qtl-anchoring-matched
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --array=0-25
#SBATCH --output=05_genetic_anchoring/_m/logs/qtl-anchoring-matched-%A_%a.log

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 05_genetic_anchoring/_m/logs

# Anchor the matched WGCNA baselines (same switch features as IsoGraph) across the
# 13 GTEx aging tissues; array index -> (method, region) over 2 methods x 13 tissues.
METHODS=(wgcna_switch_only wgcna_multiplex)
REGIONS=(
    amygdala anterior_cingulate_cortex_ba24 caudate_basal_ganglia
    cerebellar_hemisphere cerebellum cortex frontal_cortex_ba9 hippocampus
    hypothalamus nucleus_accumbens_basal_ganglia putamen_basal_ganglia
    spinal_cord_cervical_c_1 substantia_nigra
)
IDX="${SLURM_ARRAY_TASK_ID:-0}"
METHOD="${METHODS[$((IDX / ${#REGIONS[@]}))]}"
REGION="${REGIONS[$((IDX % ${#REGIONS[@]}))]}"

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** xQTL anchoring (matched): ${METHOD} gtex-aging ${REGION} ****"
python -m isograph_benchmark.real_data.qtl_anchoring \
    --analysis gtex-aging --region "${REGION}" --method "${METHOD}" "$@"
conda deactivate
log_message "**** Complete ****"
