#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-signal-susie
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --array=0-77
#SBATCH --output=05_genetic_anchoring/_m/logs/coloc-signal-susie-%A_%a.log
#
# Signal-level coloc, stage B: QTL SuSiE + coloc.susie, one (analysis, tissue) task.
# 6 analyses x 13 GTEx brain tissues = 78 tasks, enumerated in the SAME order as the
# work_list.tsv that `--stage prep` writes, so a task id always names the same cell.
#
# Memory on PSC is --cpus-per-task x 2000MB; 16 cpus = 32 GB. One locus LD matrix is held
# at a time (capped at MAX_SNPS = 12,000 -> 1.15 GB as float64) plus one chromosome of
# one tissue's all-pairs for the target genes, which is what the abf layer runs in 16 GB.
#
# Extraction from the 273 GB of GTEx brain all-pairs is by parquet predicate pushdown, so
# the release is never scanned. No per-cell estimate_s_rss: it is an O(p^3) eigen and
# there are ~46k cells (see the R script).
#
# Depends on: 03c.coloc_gwas_susie.sh (the GWAS SuSiE cache) and
#             `coloc_signal_susie --stage prep` (targets + GTEx credible sets).
#
# Gene pool and sQTL phenotype choice are env-selected, per submission:
#   sbatch --export=ALL,COLOC_SIGNAL_ARM=background  05_genetic_anchoring/_h/06b.coloc_signal_susie.sh
#   sbatch --export=ALL,COLOC_SIGNAL_SQTL=all        05_genetic_anchoring/_h/06b.coloc_signal_susie.sh
#
# Each axis writes to its own directory, because both modes name their shards
# <analysis>__<tissue>.parquet: arm -> arms/<arm>/, sqtl=all -> all_introns/. The meta
# stage must be told which set to read, and the flag has to match the env var:
#   python -m isograph_benchmark.real_data.coloc_signal_susie --stage meta --sqtl all
#
# The all-introns arm needs MORE than the directives below and must override them. GTEx
# carries ~12 intron phenotypes per gene, so the sQTL side goes from one SuSiE fit per
# gene to twelve (~6.5x the total work, measured on a chr22 smoke test). Submit it as:
#   sbatch --export=ALL,COLOC_SIGNAL_SQTL=all --time=24:00:00 --cpus-per-task=24 \
#          05_genetic_anchoring/_h/06b.coloc_signal_susie.sh
# The 12h/16cpu defaults are sized for the representative arm, whose slowest task ran
# 66 min; the same task all-introns projects to ~7h, which is too close to that ceiling.
#
# GWAS SNP-guard sensitivity arms. This stage reads COLOC_GWAS_MAX_SNPS exactly as stage A
# does, so a stage-A cache built at a non-default guard is used only by a stage-B run
# exporting the same value, and both land under _m/coloc_signal_susie/sensitivity/
# max_snps_<N>/ instead of on top of the primary shards. aging__ad is tasks 0-12:
#   sbatch --export=ALL,COLOC_GWAS_MAX_SNPS=30000 --array=0-12 --cpus-per-task=32 \
#          05_genetic_anchoring/_h/06b.coloc_signal_susie.sh
#   python -m isograph_benchmark.real_data.coloc_signal_susie --stage meta --max-snps 30000
# Usage: sbatch 05_genetic_anchoring/_h/06b.coloc_signal_susie.sh
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
mkdir -p 05_genetic_anchoring/_m/logs

export COLOC_SIGNAL_ARM="${COLOC_SIGNAL_ARM:-switch}"
export COLOC_SIGNAL_SQTL="${COLOC_SIGNAL_SQTL:-representative}"

BASE=05_genetic_anchoring/_m/coloc_signal_susie
WORK="${BASE}/work_list.tsv"
if [[ "${COLOC_SIGNAL_ARM}" != "switch" ]]; then
    WORK="${BASE}/arms/${COLOC_SIGNAL_ARM}/work_list.tsv"
fi
[[ -f "${WORK}" ]] || { echo "ERROR: ${WORK} not found; run --stage prep first."; exit 1; }

IDX="${SLURM_ARRAY_TASK_ID:-${1:-0}}"
# +2 skips the header and turns the 0-based array id into a 1-based data row.
ROW=$(( IDX + 2 ))
ANALYSIS=$(awk -v r="${ROW}" 'NR==r{print $1}' "${WORK}")
TISSUE=$(awk -v r="${ROW}" 'NR==r{print $2}' "${WORK}")
[[ -n "${ANALYSIS}" && -n "${TISSUE}" ]] || { echo "ERROR: no work row ${ROW}"; exit 1; }

# Guard: some Bridges2 batch nodes start array tasks without Lmod initialised.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log "**** coloc.susie: ${ANALYSIS} / ${TISSUE} (task ${IDX}, arm ${COLOC_SIGNAL_ARM}, sqtl ${COLOC_SIGNAL_SQTL}) ****"
Rscript 05_genetic_anchoring/_h/06b.coloc_signal_susie.R "${ANALYSIS}" "${TISSUE}"
conda deactivate
log "**** Complete: ${ANALYSIS} / ${TISSUE} ****"
