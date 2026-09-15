#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=module-coloc-conv
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=05_genetic_anchoring/_m/logs/module-coloc-conv-%j.log
#
# Module-level colocalization convergence for all five traits (AD, PD, LBD, ALS, SCZ).
#
# Generalises the SCZ-only convergence layer in scz_age_projection.py so the four aging
# traits are reported on the same partitions with the same statistics, and adds the two
# comparisons the count-only version lacked: a size-matched permutation null, and the
# CLPP-tested pool as the denominator for the anchored-module enrichment.
#
# Single process, pure joins + a 20k permutation over small vectors; runs in a couple of
# minutes. Deterministic (seed 13).

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
mkdir -p 05_genetic_anchoring/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** Module-level coloc convergence (all traits) ****"
python -m isograph_benchmark.real_data.module_coloc_convergence "$@"
conda deactivate
log_message "**** Complete ****"
