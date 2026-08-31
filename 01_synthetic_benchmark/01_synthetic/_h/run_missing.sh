#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=synth-missing
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --array=1-323%60
#SBATCH --time=02:00:00
#SBATCH --output=01_synthetic_benchmark/01_synthetic/_m/logs/%x-%A_%a.log

## Completes the main synthetic grid: runs every declared row that has no done.json,
## EXCEPT the 1,180 isograph_vae_residual rows on the unconfounded switch scenarios, which
## are deliberately not run (they would be silent no-ops on datasets whose sample tables
## cannot be repaired -- see the comment at the top of configs/synthetic_grid.yaml).
##
## The 323 rows here are:
##   312  isograph_vae_residual x abundance_switch_mixed  (datasets refreshed + reproducing,
##        so residualization is a real operation; a genuine coverage gap in the multiplex
##        feature space)
##    11  genetic_anchoring stragglers across spearman_leiden/vae/multiplex, which are why
##        that scenario currently has unequal n across methods
##
## Every dataset these rows need already exists, so ensure_dataset short-circuits and the
## non-atomic save_dataset_bundle race cannot fire. If that ever stops being true,
## materialise serially first (cf. residual_cost.py materialize).
##
## Prerequisite: 01_synthetic_benchmark/01_synthetic/_m/missing_runs.tsv (regenerate with
##   python -m isograph_benchmark.benchmark.list_missing_runs)
## Follow-up: re-run the metrics chain (step_1_collect.sh -> step_2_summarize.sh ->
##   step_3_figures.sh) so the new runs reach the tables and figures.

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
mkdir -p 01_synthetic_benchmark/01_synthetic/_m/logs

LIST="01_synthetic_benchmark/01_synthetic/_m/missing_runs.tsv"
if [[ ! -f "${LIST}" ]]; then
    echo "ERROR: ${LIST} missing. Run: python -m isograph_benchmark.benchmark.list_missing_runs"
    exit 1
fi

# +1 skips the header; array is 1-based.
LINE=$(( ${SLURM_ARRAY_TASK_ID:-1} + 1 ))
RUN_ID="$(awk -v n="${LINE}" 'NR==n{print $1}' "${LIST}")"
if [[ -z "${RUN_ID}" ]]; then
    echo "ERROR: no run_id at line ${LINE} of ${LIST}"
    exit 1
fi

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** missing run ${SLURM_ARRAY_TASK_ID:-1}: ${RUN_ID} ****"
python -u -m isograph_benchmark.benchmark.run_one --run-id "${RUN_ID}" "$@"
conda deactivate
log_message "**** Complete (${RUN_ID}) ****"
