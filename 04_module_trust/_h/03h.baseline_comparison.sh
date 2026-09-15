#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=baseline-comparison
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --output=04_module_trust/_m/logs/%x-%j.log

## Three-baseline synthesis: per-module phenotype and GO rates for isograph vs
## wgcna_gene / wgcna_switch_only / wgcna_multiplex, pooled and per region. This is the
## analysis that BOUNDS the paper's claim -- it is what shows IsoGraph is not globally
## superior, and that phenotype signal lives in the switch features rather than the
## inference. Rates, never totals: module counts differ by an order of magnitude between
## methods, so totals would compare granularity rather than quality.
##
## Requires 03_module_enrichment_*.sh (both cohorts) and the matched-feature WGCNA
## baselines from 02_module_discovery/_h/01g-01i to have run first. Light enough for the
## login node, but committed as a wrapper so the published synthesis is reproducible from
## a single command.
## Usage: sbatch 04_module_trust/_h/03h.baseline_comparison.sh

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
mkdir -p 04_module_trust/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** three-baseline comparison starts ****"
python -u -m isograph_benchmark.real_data.baseline_comparison "$@"
conda deactivate
log_message "**** Complete ****"
