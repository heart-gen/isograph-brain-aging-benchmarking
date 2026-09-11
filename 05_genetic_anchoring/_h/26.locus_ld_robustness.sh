#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=locus-ld-robustness
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/locus-ld-robustness-%j.log
#
# Locus-specific LD robustness audit, for a locus that is about to carry a biological
# claim. Five checks: convergence, credible-set purity, boundary perturbation, coloc
# stability, and an LD-mismatch-aware SuSiE-RSS pass (kriging_rss outlier drop +
# estimate_residual_variance = TRUE).
#
# WHEN THIS IS REQUIRED. The primary coloc grid is uniform at MAX_SNPS = 12,000. Loci
# recovered by raising that guard are a scoped SENSITIVITY arm, and they were empirically
# enriched for GWAS-reference-LD inconsistency (aging__ad: all 11 recovered loci at
# s_rss >= 0.310 vs a median of 0.255). So a recovered locus may NOT be promoted into the
# biological narrative on its coloc.susie posterior alone -- it runs through here first,
# and a high-s_rss locus needs substantially more than this before it is a headline.
#
# Memory: 32 cpus = 64 GB. kriging_rss is an O(p^3) eigen on the full locus LD and the
# script holds several fits at once; at p ~ 16k the LD matrix alone is 2 GB and the AD
# stage-A run at this size peaked at 48 GB.
#
# Usage: sbatch 05_genetic_anchoring/_h/26.locus_ld_robustness.sh \
#            <analysis> <LOCUS_ID> <gene_bare> <tissue[,tissue...]>
# Several tissues go in ONE run, comma-separated: the output dir is per locus, so a second
# run for another tissue would overwrite the first, and only check 4 depends on tissue.
# PICALM (the locus this was written for):
#   sbatch 05_genetic_anchoring/_h/26.locus_ld_robustness.sh \
#          aging__ad locus60_chr11 ENSG00000073921 Brain_Cortex
# UNC13A, which colocalizes in both cerebellar tissues:
#   sbatch 05_genetic_anchoring/_h/26.locus_ld_robustness.sh \
#          aging__als locus51_chr19 ENSG00000130477 Brain_Cerebellum,Brain_Cerebellar_Hemisphere
#
# Depends on: the per-locus GWAS, LD and exclusion files under _m/coloc/<analysis>/susie/
# only. The GWAS SuSiE is re-fit here from those inputs with no SNP guard, so this does NOT
# read the stage-22 cache and gives the same answer whether or not a recovery arm exists.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
mkdir -p 05_genetic_anchoring/_m/logs

[[ $# -ge 4 ]] || { echo "ERROR: need <analysis> <LOCUS_ID> <gene_bare> <tissue>"; exit 1; }

# Guard: some Bridges2 batch nodes start array tasks without Lmod initialised.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log "**** locus LD robustness: $1 / $2 / $3 / $4 ****"
Rscript 05_genetic_anchoring/_h/26.locus_ld_robustness.R "$1" "$2" "$3" "$4"
conda deactivate
log "**** Complete: $2 ****"
