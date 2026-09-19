#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=gene-universe
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=02_module_discovery/_m/logs/%x-%j.log

# Write the production gene universe (the genes surviving the production transcript
# filter) for every cohort x region, so the gene-level WGCNA baseline can be fit on
# the same universe IsoGraph models instead of the bundle's wider gene list.
# Runs before the fits in tier 01: 01d/01e/01f read what it writes.

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "ERROR: project root does not exist: ${PROJECT_ROOT}"
    exit 1
fi
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the isograph-brain-aging-benchmarking repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 02_module_discovery/_m/logs

log_message "**** production gene universe job starts ****"
echo "Job id: ${SLURM_JOBID:-local} | Node: ${SLURM_NODENAME:-local} | Host: ${HOSTNAME}"

# Guard: some Bridges2 batch nodes start without Lmod initialised.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
module list

log_message "Activating IsoGraph environment"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

python -m isograph_benchmark.real_data.production_gene_universe "$@"

conda deactivate
log_message "**** production gene universe job ends ****"
