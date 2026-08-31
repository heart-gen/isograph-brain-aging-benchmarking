#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=residual-cost
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --array=1-200%50
#SBATCH --time=02:00:00
#SBATCH --output=01_synthetic_benchmark/01_synthetic/_m/logs/%x-%A_%a.log

## Residual-cost ablation: isograph_vae vs isograph_vae_residual on 100 FRESH unconfounded
## datasets, paired within dataset (both arms see byte-identical input; only the
## residualization flag differs). Closes the open benchmark gap "is residualization free
## when there is no confound to remove?".
##
## Uses its own grid, its own dataset root and its own run root, so it cannot touch the
## archived benchmark datasets -- 1,204 of which no longer regenerate from the current
## generator (01_synthetic_benchmark/01_synthetic/_m/SAMPLE_TABLE_REFRESH.md). Model fitting and metrics
## are the production run_one code path, unchanged.
##
## Scope: RIN/neuron_frac/batch are constants here and build_design_matrix drops them, so
## this measures the cost of residualizing LIBRARY DEPTH ONLY.
##
## Prerequisite: python -m isograph_benchmark.benchmark.residual_cost grid
## Follow-up:    python -m isograph_benchmark.benchmark.residual_cost summarize

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

GRID="01_synthetic_benchmark/00_design/_m/residual_cost_grid.parquet"
if [[ ! -f "${GRID}" ]]; then
    echo "ERROR: ${GRID} missing. Run: python -m isograph_benchmark.benchmark.residual_cost grid"
    exit 1
fi

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

IDX=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))
RUN_ID="$(python -c "
import pandas as pd
print(pd.read_parquet('${GRID}').iloc[${IDX}]['run_id'])
")"

log_message "**** residual cost: task ${SLURM_ARRAY_TASK_ID:-1} run_id ${RUN_ID} ****"
python -u -m isograph_benchmark.benchmark.run_one \
    --grid "${GRID}" \
    --run-id "${RUN_ID}" \
    --dataset-root 01_synthetic_benchmark/01_synthetic/_m/datasets_residual_cost \
    --output-root 01_synthetic_benchmark/01_synthetic/_o/runs_residual_cost "$@"
conda deactivate
log_message "**** Complete (${RUN_ID}) ****"
