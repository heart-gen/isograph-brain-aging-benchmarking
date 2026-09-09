#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=stability-isograph
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=24  # 24 x 2000M = 48G; dlpfc (idx 3) peaks ~46G, OOMs at 32G
#SBATCH --time=04:00:00
#SBATCH --array=1-6
#SBATCH --output=03_module_trust/_m/logs/%x-%A_%a.log

# IsoGraph within-cohort split-half stability for the 3 matched regions in each
# cohort (6 array tasks). Each task refits IsoGraph on both 50/50 sample halves for
# every seed and writes per-half partitions for the aggregator. Heavy compute ->
# SLURM only.

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 03_module_trust/_m/logs

# array index -> (cohort, region)
SPECS=(
    "brainseq caudate"
    "brainseq hippocampus"
    "brainseq dlpfc"
    "gtex caudate_basal_ganglia"
    "gtex hippocampus"
    "gtex frontal_cortex_ba9"
)
spec="${SPECS[$((SLURM_ARRAY_TASK_ID - 1))]}"
read -r COHORT REGION <<< "${spec}"

SEEDS="${STABILITY_SEEDS:-5}"
# Consensus Leiden over N seeded, edge-weighted runs (>=2 enables; partitions are tagged
# 'isograph_consensus' for an A/B against the baseline 'isograph' partitions). Set via
#   sbatch --export=ALL,STABILITY_CONSENSUS=10 ...
CONSENSUS="${STABILITY_CONSENSUS:-1}"
# Covariate-free isoform-estimability switch-edge downweighting (1 enables; partitions
# tagged 'isograph_reliability' for an A/B against the baseline 'isograph' partitions).
# Set via:  sbatch --export=ALL,STABILITY_RELIABILITY=1 ...
RELIABILITY="${STABILITY_RELIABILITY:-0}"
RELIABILITY_FLAG=""
if [[ "${RELIABILITY}" == "1" ]]; then RELIABILITY_FLAG="--reliability"; fi
# Reliability-weight floor (reliability in [F,1]); caps downweighting / guards n_common.
FLOOR_FLAG=""
if [[ -n "${STABILITY_RELIABILITY_FLOOR:-}" ]]; then FLOOR_FLAG="--reliability-floor ${STABILITY_RELIABILITY_FLOOR}"; fi
# Giant-module cap (post-hoc split, negative lever): recursively re-cluster any module over
# this fraction of assigned genes (tags 'isograph_capNN'). Set via:
#   sbatch --export=ALL,STABILITY_MAX_MODULE_FRAC=0.15 ...
CAP_FLAG=""
if [[ -n "${STABILITY_MAX_MODULE_FRAC:-}" ]]; then CAP_FLAG="--max-module-frac ${STABILITY_MAX_MODULE_FRAC}"; fi
# Collapse fix C — resolution-sweep giant cap: pick the smallest Leiden resolution whose
# largest module is <= this fraction of genes (tags 'isograph_gcapNN' for an A/B vs baseline).
# Recommended candidate 0.15. Set via:
#   sbatch --export=ALL,STABILITY_LEIDEN_GIANT_FRAC=0.15 ...
GCAP_FLAG=""
if [[ -n "${STABILITY_LEIDEN_GIANT_FRAC:-}" ]]; then GCAP_FLAG="--leiden-giant-frac ${STABILITY_LEIDEN_GIANT_FRAC}"; fi
# Persist each split half's gene-gene graph, so 16.stability_resolution_sweep.sh can
# re-cluster it across a resolution grid instead of refitting the VAE per resolution.
# Costs ~4 MB per half and nothing in compute. Set via:
#   sbatch --export=ALL,STABILITY_SAVE_EDGES=1 ...
SAVE_EDGES_FLAG=""
if [[ "${STABILITY_SAVE_EDGES:-0}" == "1" ]]; then SAVE_EDGES_FLAG="--save-edges"; fi

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** IsoGraph split-half: ${COHORT}/${REGION} (${SEEDS} seeds, consensus=${CONSENSUS}, reliability=${RELIABILITY}) starts ****"
# One process per FIT (seed x half): each VAE fit's allocations are not fully reclaimed
# in-process, and a single fit on the largest region (~18k genes) peaks near 16G, so
# even two fits in one process OOM-kill at 32G. Isolating to one fit per process keeps
# peak to a single fit; process exit reclaims everything. Bundle is re-read per fit
# (cheap relative to the fit) — the cost of the safe path.
for ((k = 0; k < SEEDS; k++)); do
    for half in A B; do
        log_message "  seed ${k} half ${half} ..."
        python -m isograph_benchmark.real_data.stability fit-isograph \
            --cohort "${COHORT}" --region "${REGION}" --seed "${k}" --half "${half}" \
            --consensus "${CONSENSUS}" ${RELIABILITY_FLAG} ${FLOOR_FLAG} ${CAP_FLAG} ${GCAP_FLAG} \
            ${SAVE_EDGES_FLAG}
    done
done
conda deactivate
log_message "**** Complete ****"
