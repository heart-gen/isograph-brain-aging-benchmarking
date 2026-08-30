#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=sqtl-concordance
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --array=0-16
#SBATCH --output=02_module_discovery/brainseq/_m/logs/sqtl-concordance-%A_%a.log

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 02_module_discovery/brainseq/_m/logs

# array index -> (analysis, region); SCZD + 3 BrainSEQ aging + 13 GTEx aging.
# Same cohort grid as 13.qtl_anchoring.sh so the two anchoring layers are paired.
SPECS=(
    "brainseq-sczd -"
    "brainseq-aging caudate"
    "brainseq-aging hippocampus"
    "brainseq-aging dlpfc"
    "gtex-aging amygdala"
    "gtex-aging anterior_cingulate_cortex_ba24"
    "gtex-aging caudate_basal_ganglia"
    "gtex-aging cerebellar_hemisphere"
    "gtex-aging cerebellum"
    "gtex-aging cortex"
    "gtex-aging frontal_cortex_ba9"
    "gtex-aging hippocampus"
    "gtex-aging hypothalamus"
    "gtex-aging nucleus_accumbens_basal_ganglia"
    "gtex-aging putamen_basal_ganglia"
    "gtex-aging spinal_cord_cervical_c_1"
    "gtex-aging substantia_nigra"
)
read -r ANALYSIS REGION <<< "${SPECS[${SLURM_ARRAY_TASK_ID:-0}]}"
REGION_FLAG=""
if [[ "${REGION}" != "-" ]]; then REGION_FLAG="--region ${REGION}"; fi

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** sQTL direction concordance: ${ANALYSIS} ${REGION} ****"
python -m isograph_benchmark.real_data.sqtl_concordance --analysis "${ANALYSIS}" ${REGION_FLAG} "$@"
conda deactivate
log_message "**** Complete ****"
