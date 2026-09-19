#!/usr/bin/env bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=ase-junc-targets
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=02:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/ase-junc-targets-%j.log

## QUEST ONLY (PI item 10a). Isoform-specific junctions per switch pair, the switch-QTL
## lead per gene and every donor's phased genotype at it, the fetch regions and the sample
## list the 02f array reads. Inputs are the switch pairs and stage-05 QTL tables in the repo
## (git pull && git lfs pull after wave 5) plus the Quest paths in configs/data_sources.yaml.
##
##   sbatch 06_switch_mechanism/_h/01l.ase_junction_targets.sh [--region dlpfc]
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 06_switch_mechanism/_m/logs

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log "**** Job starts ****"
python -m isograph_benchmark.real_data.ase_junction_switch --stage targets "$@"
log "**** Job ends ****"
