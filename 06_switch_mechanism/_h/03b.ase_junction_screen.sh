#!/usr/bin/env bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=ase-junc-screen
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=01:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/ase-junc-screen-%j.log

## QUEST ONLY (PI item 10a). Pools the 02f shards into junction_allelic_counts.parquet (the
## one table that goes back to Bridges-2, via git-LFS) and applies the pre-registered gate:
## >= 30 donors x >= 5 informative fragments, both isoforms; the arm proceeds only if >= 30
## gate-family pairs pass. Refuses to run on missing shards unless --allow-partial.
##
##   sbatch 06_switch_mechanism/_h/03b.ase_junction_screen.sh [--region dlpfc]
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log "**** Job starts ****"
python -m isograph_benchmark.real_data.ase_junction_switch --stage screen "$@"
log "**** Job ends ****"
