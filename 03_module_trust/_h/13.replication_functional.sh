#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=repl-functional
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --output=03_module_trust/_m/stability/logs/%x-%A_%a.log

## Reviewer item 7: are the cross-cohort matched modules functionally preserved despite
## their low gene overlap (median gene Jaccard ~0.02)? Compares each matched pair's GO-term
## set, cell-type marker profile and switch-consequence class profile against SIZE-MATCHED
## random pairs (gene-count decile, 1,000 draws, seed 13). Operates on the same 130 pairs as
## the replication permutation test. Array = method x age model: 0 isograph/linear,
## 1 wgcna/linear, 2 isograph/spline, 3 wgcna/spline. Both models run because item 2's
## primary statistic is the spline one while `linear` is what the published table reports.
##
## Not every measure is computable for every method: celltype_composition exists only for
## the isograph tree, and switch_consequence skips regions with no phenotype-significant
## switch genes (3 of the 6 regions here), so structure_r has no data. The report's
## "Data coverage" section names the absent inputs -- n/a is not a null result.
## Usage: sbatch 03_module_trust/_h/13.replication_functional.sh [--n-perm N]
##
## Requires: 09.module_trust_replication.sh (the matched pairs), module_enrichment,
## celltype_composition, and 15.switch_consequence.sh.

set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 03_module_trust/_m/stability/logs

SPECS=(
    "isograph linear"
    "wgcna    linear"
    "isograph spline"
    "wgcna    spline"
)
read -r METHOD MODEL <<< "${SPECS[${SLURM_ARRAY_TASK_ID:-0}]}"

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log "**** functional preservation: ${METHOD} / ${MODEL} ****"
python -u -m isograph_benchmark.real_data.replication_functional \
    --method "${METHOD}" --model "${MODEL}" "$@"
conda deactivate
log "**** Complete (${METHOD} / ${MODEL}) ****"
