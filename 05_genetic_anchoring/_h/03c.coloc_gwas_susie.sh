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
# Depends on: 01g.coloc_prep.sh + 02e.locus_ld.sh (per-locus GWAS + LD already on disk).
#
# MAX_SNPS RECOVERY RUNS. The R script's MAX_SNPS guard (default 12,000) drops 61 of the
# 579 loci before any fit -- including PICALM's AD locus (locus60_chr11, 15,713 SNPs),
# which is why PICALM has no signal-level colocalization in either sQTL arm. Raising it
# recovers them, but memory goes as p^2 and the 16-cpu/24-cpu defaults below are sized
# for 12,000, so a recovery run must override BOTH:
#   sbatch --export=ALL,COLOC_GWAS_MAX_SNPS=30000 --cpus-per-task=32 --array=0 \
#          05_genetic_anchoring/_h/03c.coloc_gwas_susie.sh
# 30,000 covers every skipped locus (max 28,677); --array=0 scopes it to aging__ad.
# A non-default COLOC_GWAS_MAX_SNPS is a scoped SENSITIVITY arm and writes to
# _m/coloc_signal_susie/sensitivity/max_snps_<N>/, never over the primary cache. Stage B
# resolves the same root from the same variable, so submit it with the SAME --export (see
# 06b.coloc_signal_susie.sh), and read it with `--stage meta --max-snps <N>`. Before
# 2026-09-10 there was no such split, and the aging__ad recovery overwrote the primary AD
# cache and shards in place.
# Things to know before trusting what comes back. This script has no resume, so a re-run
# re-fits every locus and rewrites that root's _status.tsv. Loci under the old guard
# reproduce only up to floating point, and that is not always enough: OpenBLAS thread count
# follows --cpus-per-task and is not pinned, so a credible set sitting on the coverage
# boundary can appear in one run and vanish in the next. Measured 2026-09-11: AD
# locus16_chr2 (11,183 SNPs; one 235-variant set at coverage 0.95006 against the 0.95
# request, max PIP 0.087) fit that set at 32 cpus and in the original primary, and none at
# 24 cpus. Its three signal-level cells had PP4 <= 0.005, so no call moved -- but a
# re-run is not bit-identical. Elsewhere coloc.susie PP4 agreed to ~1e-9 across both sqtl
# arms, with 2 lead-variant labels differing at PIP ties / NaN posteriors and 1 eQTL cell
# (PP4 0.003) not re-emitted, most likely the same boundary effect on the QTL side. And loci past the guard
# were empirically enriched for GWAS-reference-LD inconsistency (aging__ad: all 11
# recovered loci at s_rss >= 0.310 vs a median of 0.255; worst 0.739) -- a correlation
# between the compute threshold and statistical difficulty, not a mechanism this run
# establishes. Read s_rss per locus before believing a recovered credible set.
# Stage B (06b.coloc_signal_susie.sh) must then be re-run to pick them up.
#
# SCOPE DECISION, 2026-09-10 (user): MAX_SNPS = 12,000 stays the UNIFORM PRIMARY grid and
# the other five analyses are NOT recovered. The aging__ad recovery is a prespecified,
# scoped SENSITIVITY arm, run to settle PICALM. Any recovered locus entering the
# biological narrative needs locus-specific LD diagnostics and an LD-mismatch-aware
# SuSiE-RSS sensitivity first -- see 05_genetic_anchoring/_h/08d.locus_ld_robustness.sh.
#
# Usage: sbatch 05_genetic_anchoring/_h/03c.coloc_gwas_susie.sh
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
Rscript 05_genetic_anchoring/_h/03c.coloc_gwas_susie.R "${ANALYSIS}"
conda deactivate
log "**** Complete: ${ANALYSIS} ****"
