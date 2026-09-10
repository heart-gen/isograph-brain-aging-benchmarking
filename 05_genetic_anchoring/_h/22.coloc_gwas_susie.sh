#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-gwas-susie
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=24
#SBATCH --time=16:00:00
#SBATCH --array=0-5
#SBATCH --output=05_genetic_anchoring/_m/logs/coloc-gwas-susie-%A_%a.log
#
# Signal-level coloc, stage A: fit and cache the per-locus GWAS SuSiE, one analysis per
# array task. Stage B (13 tissues x 6 analyses) reads this cache, so the GWAS side is fit
# once instead of 26 times per locus -- and every tissue task colocalizes against the
# identical fit rather than against its own re-convergence.
#
# Memory on PSC is --cpus-per-task x 2000MB; 24 cpus = 48 GB. An N x N float64 LD matrix
# at the MAX_SNPS cap of 12,000 is 1.15 GB and susie_rss plus the eigen decomposition in
# estimate_s_rss hold several copies, so the headroom is real rather than padding.
#
# Time: the eigen decomposition inside estimate_s_rss is ~2 min at p ~ 6k and dominates
# everything else (susie_rss itself is ~8 s). It is therefore run ONLY for loci that
# actually yielded a credible set. aging__scz is the long pole at 274 loci.
#
# Depends on: 08.coloc_prep.sh + 09.locus_ld.sh (per-locus GWAS + LD already on disk).
# Usage: sbatch 05_genetic_anchoring/_h/22.coloc_gwas_susie.sh
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
mkdir -p 05_genetic_anchoring/_m/logs

ANALYSES=(aging__ad aging__als aging__lbd aging__pd aging__scz brainseq-sczd__scz)
IDX="${SLURM_ARRAY_TASK_ID:-${1:-0}}"
ANALYSIS="${ANALYSES[${IDX}]}"

# Guard: some Bridges2 batch nodes start array tasks without Lmod initialised.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log "**** GWAS SuSiE cache: ${ANALYSIS} (task ${IDX}) ****"
Rscript 05_genetic_anchoring/_h/22.coloc_gwas_susie.R "${ANALYSIS}"
conda deactivate
log "**** Complete: ${ANALYSIS} ****"
