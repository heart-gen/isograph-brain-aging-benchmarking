#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ldsc-annot-prep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=real_data/ldsc/_m/logs/ldsc-annot-prep-%j.log

## S-LDSC step 1 (wrapper) — build the sQTL/eQTL/cis switch-gene SNP annotations
## (hg19 BEDs) for a single analysis or a pooled bundle. Reproducible entry point for
## ldsc_annot_prep.py (was previously run interactively).
##
##   disease : sbatch real_data/ldsc/_h/01.annot_prep.sh --analysis brainseq-sczd
##   aging   : sbatch real_data/ldsc/_h/01.annot_prep.sh --bundle aging --min-recurrence 3
##
## Then 02.make_annot_ldscores.sh <annot> and 03.munge_h2.sh <trait> <annot>.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p real_data/ldsc/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log "**** ldsc_annot_prep $* ****"
python -m isograph_benchmark.real_data.ldsc_annot_prep "$@"
conda deactivate
log "**** ldsc annot prep done ****"
