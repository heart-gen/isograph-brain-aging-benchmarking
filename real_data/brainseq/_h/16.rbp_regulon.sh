#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=rbp-regulon
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=real_data/brainseq/_m/logs/rbp-regulon-%j.log

## RBP-regulon analysis (light 3'UTR/mature-transcript scope), two stages:
##   stage 1 (motif env): rbp_scan.py — scan switch-isoform sequences (GENCODE v47 transcript
##            FASTA) against ATtRACT human RBP PWMs (MOODS) -> per-(transcript,RBP) hit counts.
##   stage 2 (isograph env): rbp_regulon.py — call RBP site gain/loss between switch-pair
##            isoforms, per-module hypergeometric regulon enrichment vs the switch-gene pool.
## Resources staged under inputs/rbp_motifs (ATtRACT) + inputs/raw/gencode_v47 (transcript FASTA).
## Usage: sbatch real_data/brainseq/_h/16.rbp_regulon.sh
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/brainseq/_m/logs

MOTIF_ENV=/ocean/projects/bio260021p/shared/opt/envs/motif
ISO_ENV=/ocean/projects/bio260021p/shared/opt/envs/isograph

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1

log "**** stage 1: MOODS motif scan (motif env) ****"
conda activate "${MOTIF_ENV}"
python -m isograph_benchmark.real_data.rbp_scan
conda deactivate

log "**** stage 2: per-module RBP regulon enrichment (isograph env) ****"
conda activate "${ISO_ENV}"
python -m isograph_benchmark.real_data.rbp_regulon
conda deactivate
log "**** RBP regulon done ****"
