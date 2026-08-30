#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=tiers-pilot-caudate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=24   # 24 x 2000M = 48G; caudate VAE fit ~30G, similarity ~5.5G
#SBATCH --time=03:00:00
#SBATCH --output=02_module_discovery/brainseq/_m/logs/%x-%j.log

# One-region pilot for the 4-tier multiplex hierarchy. Step 1 re-fits the caudate
# VAE ONCE (now emitting feature_reconstruction.parquet); step 2 projects that single
# fit into the three IsoGraph tiers (switch_only / switch_primary / full_multiplex)
# without re-fitting. WGCNA (the 4th tier, gene-abundance baseline) is run separately.
# Heavy compute -> SLURM only.
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
mkdir -p 02_module_discovery/brainseq/_m/logs

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

log_message "**** PILOT step 1: re-fit caudate VAE (emits feature_reconstruction) ****"
python -m isograph_benchmark.real_data.run_models brainseq-aging --region caudate

log_message "**** PILOT step 2: project caudate fit into 3 IsoGraph tiers ****"
python -m isograph_benchmark.real_data.project_tiers brainseq-aging --region caudate

conda deactivate
log_message "**** PILOT complete ****"
