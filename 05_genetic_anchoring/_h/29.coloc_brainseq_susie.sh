#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-brainseq
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --array=0-17
#SBATCH --output=05_genetic_anchoring/_m/logs/coloc-brainseq-%A_%a.log
#
# BrainSEQ signal-level coloc: S_g and A_g QTL SuSiE on in-sample EA LD + coloc.susie, with
# coloc.abf for every cell. One task per row of work_list.tsv -- 6 analyses x 3 regions = 18,
# in the order `coloc_brainseq --stage prep` writes them, so a task id always names one cell.
#
# Same-tissue genetic anchoring, not replication. EA-only by construction: prep refuses any
# other arm, and refuses ea_only unless both BrainSEQ QTL checks passed (28.brainseq_qtl_checks.sh).
#
# Memory on PSC is --cpus-per-task x 2000MB; 16 cpus = 32 GB. One locus LD matrix is held at
# a time (at most the 12,000-SNP stage-A guard -> 1.15 GB as float64, briefly twice while it is
# built), plus one chromosome of the region's nominal pairs for the target genes. plink2 gets
# half the allocation for the LD computation.
#
# Order:
#   python -m isograph_benchmark.real_data.coloc_brainseq --stage prep          # login
#   sbatch --array=0 --export=ALL,COLOC_BRAINSEQ_CHR=22 05_genetic_anchoring/_h/29.coloc_brainseq_susie.sh  # smoke
#   sbatch 05_genetic_anchoring/_h/29.coloc_brainseq_susie.sh
#   python -m isograph_benchmark.real_data.coloc_brainseq --stage meta          # login
#
# A COLOC_BRAINSEQ_CHR run writes under _m/coloc_brainseq/<arm>/smoke/, which meta never reads.
#
# Depends on: 22.coloc_gwas_susie.sh (GWAS SuSiE cache), `coloc_signal_susie --stage prep`
# (target grid), 25.brainseq_switch_qtl.sh --export=ALL,SWQTL_ARM=ea_only (nominal pairs, kept
# on disk though gitignored) and 28.brainseq_qtl_checks.sh for that arm.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
mkdir -p 05_genetic_anchoring/_m/logs

export COLOC_BRAINSEQ_ARM="${COLOC_BRAINSEQ_ARM:-ea_only}"
WORK="05_genetic_anchoring/_m/coloc_brainseq/${COLOC_BRAINSEQ_ARM}/work_list.tsv"
[[ -f "${WORK}" ]] || { echo "ERROR: ${WORK} not found; run --stage prep first."; exit 1; }

IDX="${SLURM_ARRAY_TASK_ID:-${1:-0}}"
# +2 skips the header and turns the 0-based array id into a 1-based data row.
ROW=$(( IDX + 2 ))
ANALYSIS=$(awk -v r="${ROW}" 'NR==r{print $1}' "${WORK}")
REGION=$(awk -v r="${ROW}" 'NR==r{print $2}' "${WORK}")
[[ -n "${ANALYSIS}" && -n "${REGION}" ]] || { echo "ERROR: no work row ${ROW}"; exit 1; }
ONE_CHR="${COLOC_BRAINSEQ_CHR:-}"

# Guard: some Bridges2 batch nodes start array tasks without Lmod initialised.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
export PLINK2=/ocean/projects/bio260021p/shared/opt/envs/eqtl/bin/plink2

log "**** BrainSEQ coloc: ${ANALYSIS} / ${REGION} (task ${IDX}, arm ${COLOC_BRAINSEQ_ARM}${ONE_CHR:+, chr${ONE_CHR} smoke}) ****"
Rscript 05_genetic_anchoring/_h/29.coloc_brainseq_susie.R "${ANALYSIS}" "${REGION}" ${ONE_CHR:+"${ONE_CHR}"}
conda deactivate
log "**** Complete: ${ANALYSIS} / ${REGION} ****"
