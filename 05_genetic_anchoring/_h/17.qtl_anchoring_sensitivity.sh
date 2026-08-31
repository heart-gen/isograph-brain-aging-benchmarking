#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=qtl-anchor-sens
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --array=0-50
#SBATCH --output=05_genetic_anchoring/_m/logs/qtl-anchor-sens-%A_%a.log
#
# Pre-specified sensitivity arms for the xQTL anchoring result. The PRIMARY analysis
# (binary sGene/eGene outcome, standard QTL-detectability covariates) is produced by
# 01.qtl_anchoring.sh and is NOT touched here — these arms write to their own
# `qtl_anchoring_<arm>.parquet` files.
#
#   constraint  binary outcome + gnomAD v4.1 LOEUF, missense z and log tissue-matched
#               expression. Tests the leading alternative explanation: co-switch module
#               genes are cis-QTL DEPLETED in BOTH modalities, which is what selective
#               constraint on network-central genes would produce. Emits BOTH covariate
#               sets fitted on the identical constraint-complete gene subset, so the
#               comparison is a nested-model test rather than a change of universe.
#   continuous  rank-INT of -log10(pval_beta) — the permutation statistic the sGene call
#               thresholds away. Threshold-free; adds precision without adding data.
#   dose        Poisson on the SuSiE credible-set count — how many INDEPENDENT cis
#               signals a gene carries, LD-resolved.
#
# 17 analyses x 3 arms = 51 tasks. ~15 min each; the dose arm reads the SuSiE parquets.

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
# arm index -> "--outcome X --covariate-set Y"
ARMS=(
    "--outcome binary --covariate-set constraint"
    "--outcome continuous --covariate-set standard"
    "--outcome dose --covariate-set standard"
)

TASK="${SLURM_ARRAY_TASK_ID:-0}"
N_SPECS=${#SPECS[@]}
read -r ANALYSIS REGION <<< "${SPECS[$(( TASK % N_SPECS ))]}"
ARM="${ARMS[$(( TASK / N_SPECS ))]}"
REGION_FLAG=""
if [[ "${REGION}" != "-" ]]; then REGION_FLAG="--region ${REGION}"; fi

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** xQTL anchoring sensitivity: ${ANALYSIS} ${REGION} [${ARM}] ****"
python -m isograph_benchmark.real_data.qtl_anchoring \
    --analysis "${ANALYSIS}" ${REGION_FLAG} ${ARM} "$@"
conda deactivate
log_message "**** Complete ****"
