#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=locus-event-audit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/locus-event-audit-%j.log
#
# Event-level audit of the colocalizing loci: does the colocalizing intron name the splice
# event the literature implicates, or merely the gene?
#
# The audit is deliberately IDENTICAL across nomination sources, so the coloc.abf result
# and the signal-level coloc.susie result are compared on the same rules rather than on
# two separately-tuned analyses:
#   sbatch 05_genetic_anchoring/_h/08b.locus_event_audit.sh abf      # current grid
#   sbatch 05_genetic_anchoring/_h/08b.locus_event_audit.sh susie    # after 23.* + meta
# Arguments after the nomination source are forwarded, which is how the primary (all-introns)
# audit and the SNP-guard sensitivity audit are run:
#   sbatch 05_genetic_anchoring/_h/08b.locus_event_audit.sh susie --sqtl-arm all
#   sbatch 05_genetic_anchoring/_h/08b.locus_event_audit.sh susie --sqtl-arm all --max-snps 30000
#
# Tiers are gated on `configs/known_splice_events.yaml`: only a `reviewed` curated event
# (coordinates verified in-repo against GENCODE v47) can promote a locus to
# `known_mechanism_recovered`. A `proposed` record is reported and never changes a tier.
#
# Memory on PSC is --cpus-per-task x 2000MB; 8 cpus = 16 GB, which holds the GENCODE v47
# exon cache (3.45M rows) plus the nomination and cell tables.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
mkdir -p 05_genetic_anchoring/_m/logs

NOMINATIONS="${1:-abf}"

# Guard: some Bridges2 batch nodes start without Lmod initialised.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log "**** locus event audit: nominations=${NOMINATIONS} ${*:2} ****"
python -m isograph_benchmark.real_data.locus_event_audit --nominations "${NOMINATIONS}" "${@:2}"
conda deactivate
log "**** Complete ****"
