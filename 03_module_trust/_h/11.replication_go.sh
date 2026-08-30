#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=xcohort-replication-go
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=03_module_trust/_m/replication/logs/%x-%j.log

# Cross-cohort GO consistency for replicated, age-associated modules. Requires
# the replication module-match tables (01.replication.sh) and the per-module GO
# enrichment tables (<method>_module_go.parquet) for both cohorts' matched
# regions. Submit with --dependency=afterok on the enrichment jobs.

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
mkdir -p 03_module_trust/_m/replication/logs

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** Cross-cohort GO consistency starts ****"
python -m isograph_benchmark.real_data.replication_go "$@"
conda deactivate
log_message "**** Complete ****"
