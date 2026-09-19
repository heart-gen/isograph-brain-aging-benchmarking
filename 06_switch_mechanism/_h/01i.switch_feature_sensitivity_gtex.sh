#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=switch-sensitivity-gtex
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16
#SBATCH --array=1-13
#SBATCH --time=02:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/switch-sensitivity-gtex-%A_%a.log

## Preprocessing sensitivity of the switch representation on every GTEx v11 brain region
## (pseudocount, transcript expression filter, minor-isoform threshold, identifiability).
## The quantification axis is cross-region and runs once in the aggregate step:
##
##   a=$(sbatch --parsable 06_switch_mechanism/_h/01i.switch_feature_sensitivity_gtex.sh)
##   sbatch --dependency=afterok:${a} 06_switch_mechanism/_h/01h.switch_feature_sensitivity.sh --aggregate
##
## GTEx is fit on the unfiltered transcript matrix (run_gtex_region applies no expression
## filter), so its published baseline differs from BrainSEQ's; the harness gates on
## reproducing each region's committed feature_scores before perturbing anything.
##
## Bridges memory is --cpus-per-task x 2000MB (32G here); do NOT pass --mem.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 06_switch_mechanism/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
REGIONS=(
    amygdala anterior_cingulate_cortex_ba24 caudate_basal_ganglia
    cerebellar_hemisphere cerebellum cortex frontal_cortex_ba9
    hippocampus hypothalamus nucleus_accumbens_basal_ganglia
    putamen_basal_ganglia spinal_cord_cervical_c_1 substantia_nigra
)
REGION="${REGIONS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}"

log "switch_feature_sensitivity --cohort gtex --region ${REGION} $*"
"${PY}" -m isograph_benchmark.real_data.switch_feature_sensitivity \
    --cohort gtex --region "${REGION}" "$@"
log "**** complete ****"
