#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=xcohort-replication
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=04_module_trust/_m/logs/%x-%j.log

# BrainSEQ vs GTEx cross-cohort replication for the three matched brain regions
# (caudate, hippocampus, dlpfc/BA9): module preservation (cross-cohort best-match
# Jaccard + permutation null) and age-effect concordance, for both methods.
#
# Requires module fits + age_linear for both cohorts and both methods. The WGCNA
# arm needs the corrected GTEx WGCNA re-run (02.wgcna_gene) to have finished;
# submit this with --dependency=afterok:<wgcna_jobid> to chain it.

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

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** Cross-cohort replication starts ****"
python -m isograph_benchmark.real_data.replication "$@"
conda deactivate
log_message "**** Complete ****"
