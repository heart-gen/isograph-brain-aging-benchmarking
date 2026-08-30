#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-clpp
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=32
#SBATCH --time=03:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/coloc-clpp-%j.log
## RM-shared allocates memory per cpu (~2 GB/cpu), so memory is scaled via cpus, not
## --mem: 32 cpus ~ 64 GB, headroom for susie_rss on large per-locus LD matrices (plus
## the MAX_SNPS cap in the R script skips the pathological long-range-LD loci).

## Coloc capstone, step 3 (wrapper) — GWAS SuSiE + eCAVIAR CLPP vs sQTL/eQTL.
## Usage: sbatch 05_genetic_anchoring/_h/10.coloc_clpp.sh brainseq-sczd
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

ANALYSIS="${1:-brainseq-sczd}"
PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
mkdir -p 05_genetic_anchoring/_m/logs

# Some Bridges2 batch nodes start without lmod initialised; source it defensively.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log "**** coloc CLPP: ${ANALYSIS} ****"
Rscript 05_genetic_anchoring/_h/10.coloc_clpp.R "${ANALYSIS}"

log "Summarizing"
conda deactivate
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph
python -m isograph_benchmark.real_data.coloc_summary --analysis "${ANALYSIS}"
conda deactivate
log "**** coloc CLPP + summary done ****"
