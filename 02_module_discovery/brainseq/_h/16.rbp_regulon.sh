#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=rbp-regulon
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8   # 8 x 2000M = 16G; sequences for ~80k transcripts are
                            # held in memory so each composition bin scans with one
                            # threshold set
#SBATCH --time=08:00:00
#SBATCH --output=02_module_discovery/brainseq/_m/logs/rbp-regulon-%j.log

## RBP-regulon analysis (light 3'UTR/mature-transcript scope), two stages:
##   stage 1 (motif env): rbp_scan.py — scan switch-isoform sequences (GENCODE v47 transcript
##            FASTA) against ATtRACT human RBP PWMs (MOODS), with GC-binned composition
##            backgrounds and a 5'UTR/CDS/3'UTR partition -> per-(transcript,RBP,region)
##            and per-(transcript,family,region) hit counts.
##   stage 2 (isograph env): rbp_regulon.py — call RBP site gain/loss between switch-pair
##            isoforms, per-module hypergeometric regulon enrichment vs the switch-gene pool.
## Resources staged under inputs/rbp_motifs (ATtRACT) + inputs/raw/gencode_v47 (transcript FASTA).
## Requires 02_module_discovery/brainseq/_h/31.rbp_motif_families.sh to have run first (stage 1
## tallies hits per motif family as well as per RBP).
## Usage: sbatch 02_module_discovery/brainseq/_h/16.rbp_regulon.sh [--stage all|scan|regulon] [--no-flat]
##   --stage scan     stage 1 only (the ~8 h MOODS scan)
##   --stage regulon  stage 2 only -- re-tests the frozen Stage-1 count tables without
##                    rescanning.  Use this when only the enrichment/GLM code changed; the
##                    scan output is an expensive frozen input and must not be rewritten
##                    to pick up a stage-2 fix.
## Remaining arguments go to stage 1 when it runs, and to stage 2 otherwise (so
## `--stage regulon --scope intronic --unit family_id` works).
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 02_module_discovery/brainseq/_m/logs

MOTIF_ENV=/ocean/projects/bio260021p/shared/opt/envs/motif
ISO_ENV=/ocean/projects/bio260021p/shared/opt/envs/isograph

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1

STAGE=all
if [[ "${1:-}" == "--stage" ]]; then
    STAGE="${2:?--stage needs a value: all|scan|regulon}"
    shift 2
fi
case "${STAGE}" in
    all|scan|regulon) ;;
    *) echo "ERROR: --stage must be all, scan or regulon (got '${STAGE}')."; exit 1 ;;
esac

if [[ "${STAGE}" == "all" || "${STAGE}" == "scan" ]]; then
    log "**** stage 1: MOODS motif scan (motif env) ****"
    conda activate "${MOTIF_ENV}"
    python -u -m isograph_benchmark.real_data.rbp_scan "$@"
    conda deactivate
else
    log "**** stage 1 skipped (--stage ${STAGE}); consuming the frozen count tables ****"
fi

if [[ "${STAGE}" == "all" || "${STAGE}" == "regulon" ]]; then
    log "**** stage 2: per-module RBP regulon enrichment (isograph env) ****"
    conda activate "${ISO_ENV}"
    # Only stage 1 takes the scan flags; on a stage-2-only run they belong to the regulon.
    if [[ "${STAGE}" == "regulon" ]]; then
        python -u -m isograph_benchmark.real_data.rbp_regulon "$@"
    else
        python -u -m isograph_benchmark.real_data.rbp_regulon
    fi
    conda deactivate
fi
log "**** RBP regulon done (stage ${STAGE}) ****"
