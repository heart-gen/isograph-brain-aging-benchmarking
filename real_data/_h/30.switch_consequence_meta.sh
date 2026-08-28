#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=switch-consequence-meta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --output=real_data/_m/logs/%x-%j.log

## Cross-region rollup of the switch coding-consequence enrichment. Pools the per-analysis
## log enrichments by DerSimonian-Laird random effects (with Cochran's Q / I2), separately
## for the aging and disease analyses, and demotes the Fisher-combined p to a secondary
## column. Requires real_data/brainseq/_h/15.switch_consequence.sh to have run first so the
## gene-level block-bootstrap SEs exist. Light enough for the login node, but committed as a
## wrapper so the published rollup is reproducible from a single command.
## Usage: sbatch real_data/_h/30.switch_consequence_meta.sh

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
mkdir -p real_data/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** switch-consequence cross-region meta starts ****"
python -u -m isograph_benchmark.real_data.switch_consequence_meta "$@"
conda deactivate
log_message "**** Complete ****"
