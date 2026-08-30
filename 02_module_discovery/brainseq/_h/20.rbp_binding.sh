#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=rbp-binding-overlap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --output=02_module_discovery/brainseq/_m/logs/rbp-binding-%j.log

# RBP binding-evidence overlap (ENCODE eCLIP) — Stage 3 of the RBP analysis.
# Overlaps each switch gene's switched exon interval (symmetric-difference exonic space of
# the switch-pair transcripts, GENCODE v47) with downloaded eCLIP peaks, then joins binding
# support into the significant rbp_regulon calls. Single task over all 17 region sets.
#
# Prereq: peaks fetched to inputs/raw/rbp_binding/ (run once, login node w/ network):
#   python -m isograph_benchmark.real_data.rbp_binding fetch
#
# Outputs (real_data/_m/rbp/): rbp_binding_calls.parquet, rbp_binding_regulon.parquet,
#   RBP_BINDING_SUMMARY.md.

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 02_module_discovery/brainseq/_m/logs
log_message "**** RBP binding overlap ****"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

python -m isograph_benchmark.real_data.rbp_binding run "$@"

conda deactivate
log_message "**** Complete ****"
