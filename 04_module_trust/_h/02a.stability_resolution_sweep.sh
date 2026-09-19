#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=stability-res-sweep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8   # 8 x 2000M = 16G; Leiden only, no VAE fit
#SBATCH --time=02:00:00
#SBATCH --array=1-6
#SBATCH --output=04_module_trust/_m/logs/%x-%A_%a.log

## Phenotype-blind Leiden resolution sweep on the split halves.
##
## Answers the standing objection that the canonical resolution (2.0) was chosen while
## looking at a trait: this picks it on a split-half stability curve instead, and demotes
## the GWAS giant-module argument to post-hoc confirmation.
##
## RE-CLUSTERS SAVED GRAPHS -- it does not refit anything. The gene-gene graph does not
## depend on the resolution; only the Leiden step does. Refitting per resolution would cost
## len(GRID) x the split-half fits at ~30-48G each, which is why 01a.stability_isograph.sh
## must first have been run with --save-edges:
##
##   sbatch --export=ALL,STABILITY_SAVE_EDGES=1 04_module_trust/_h/01a.stability_isograph.sh
##   sbatch 04_module_trust/_h/02a.stability_resolution_sweep.sh
##   04_module_trust/_h/03a.stability_aggregate.sh    # -> one row per resolution
##
## The canonical resolution is skipped: its partitions are the committed baseline and the
## sweep must not rewrite them. Every other resolution is tagged 'isograph_resXpY', so
## `aggregate` reports it as its own method and the curve falls out of stability_summary.
## Clustering is production's edge-weighted Leiden (sweep_leiden._build_module_table calls
## IsoGraph's NetworkModel._module_table). 5.0 is skipped when its split-half FITS exist
## (isograph_res5__*, the retired canonical arm), so that point is fitted, not re-clustered.
##
## Bridges memory is --cpus-per-task x 2000MB; do NOT pass --mem. Calls the env interpreter
## directly rather than `module load` + `conda activate`, so an array task cannot die on
## "module: command not found".
##
## Override the grid with:  sbatch --export=ALL,SWEEP_RESOLUTIONS="1 2 4 8" ...
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
mkdir -p 04_module_trust/_m/logs

SPECS=(
    "brainseq caudate"
    "brainseq hippocampus"
    "brainseq dlpfc"
    "gtex caudate_basal_ganglia"
    "gtex hippocampus"
    "gtex frontal_cortex_ba9"
)
spec="${SPECS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}"
read -r COHORT REGION <<< "${spec}"

# Geometric-ish grid bracketing the canonical 2.0 on both sides, so the curve can show an
# optimum away from it rather than only confirming it.
RESOLUTIONS="${SWEEP_RESOLUTIONS:-0.5 1 2 3 5 8 12 20}"

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "**** resolution sweep: ${COHORT}/${REGION} over [${RESOLUTIONS}] ****"
"${PY}" -m isograph_benchmark.real_data.stability sweep \
    --cohort "${COHORT}" --region "${REGION}" --resolutions ${RESOLUTIONS}
log "**** complete -- now run 03a.stability_aggregate.sh ****"
