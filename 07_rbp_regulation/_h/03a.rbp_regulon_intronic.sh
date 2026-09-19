#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=rbp-regulon-intronic
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8   # 8 x 2000M = 16G; flank sequences are held in memory
#SBATCH --time=08:00:00
#SBATCH --output=07_rbp_regulation/_m/logs/rbp-regulon-intronic-%j.log

## RBP-regulon analysis (intronic splice-site-flank scope), extends 07_rbp_regulation/_h/02a.rbp_regulon.sh from the
## mature transcript to the intronic binding niche of splicing-regulatory RBPs. Three stages:
##   stage 1 (motif env): rbp_scan_intronic.py — derive introns from the GENCODE v47 GTF for every
##            switch-isoform transcript, extract strand-aware intronic splice-site flanks from the
##            genome FASTA (pre-mRNA sense; minus strand reverse-complemented), scan against ATtRACT
##            human PWMs (MOODS, same settings as mature) -> rbp_counts_intronic.parquet.
##   stage 2 (isograph env): rbp_regulon.py --scope intronic  — intronic-only regulon enrichment.
##   stage 3 (isograph env): rbp_regulon.py --scope combined   — mature ∪ intronic presence.
## The canonical mature outputs (07_rbp_regulation/_h/02a.rbp_regulon.sh) are left untouched; this writes *_intronic /
## *_combined suffixed tables + reports. Requires rbp_scan.py's rbp_counts.parquet for --scope
## combined (run 16 first). Genome FASTA + GTF are the staged GENCODE v47 shared resources.
## Usage: sbatch 07_rbp_regulation/_h/03a.rbp_regulon_intronic.sh
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 07_rbp_regulation/_m/logs

MOTIF_ENV=/ocean/projects/bio260021p/shared/opt/envs/motif
ISO_ENV=/ocean/projects/bio260021p/shared/opt/envs/isograph

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"

log "**** stage 1: MOODS intronic-flank motif scan (motif env) ****"
conda activate "${MOTIF_ENV}"
python -u -m isograph_benchmark.real_data.rbp_scan_intronic "$@"
conda deactivate

log "**** stage 2: per-module RBP regulon enrichment, intronic scope (isograph env) ****"
conda activate "${ISO_ENV}"
python -m isograph_benchmark.real_data.rbp_regulon --scope intronic
log "**** stage 3: per-module RBP regulon enrichment, combined scope (isograph env) ****"
python -m isograph_benchmark.real_data.rbp_regulon --scope combined
conda deactivate
log "**** RBP regulon (intronic + combined) done ****"
