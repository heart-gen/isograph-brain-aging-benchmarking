#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=rbp-motif-families
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=real_data/brainseq/_m/logs/rbp-motif-families-%j.log

## Collapse the ~1,200 redundant ATtRACT human RBP matrices into motif-similarity families
## (reviewer item 6b), so downstream regulon recurrence can be reported per family rather
## than per matrix or per RBP -- N near-identical motifs are not N independent lines of
## evidence. Runs BEFORE 16.rbp_regulon.sh, whose scan tallies hits per family.
## The clustering cut is pinned in configs/rbp_families.yaml; the cut-sensitivity sweep is
## written to real_data/_m/rbp/RBP_MOTIF_FAMILIES.md.
## Usage: sbatch real_data/brainseq/_h/31.rbp_motif_families.sh [--cut 0.25]

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
mkdir -p real_data/brainseq/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log "**** RBP motif families ****"
python -u -m isograph_benchmark.real_data.rbp_motif_families "$@"
conda deactivate
log "**** Complete ****"
