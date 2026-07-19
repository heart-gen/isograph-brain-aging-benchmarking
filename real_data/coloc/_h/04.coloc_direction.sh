#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-direction
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=real_data/coloc/_m/logs/coloc-direction-%j.log

## Coloc step 4 (wrapper) — signed direction + resolved isoform events for every
## colocalized switch gene. Reproducible entry point for the two analysis CLIs that
## were previously run interactively:
##   coloc_direction.py       — align GWAS risk allele to the GTEx effect allele -> signed
##                              risk_qtl_effect; parse the LeafCutter junction for sQTL.
##   coloc_isoform_events.py  — map the junction to the tissue-matched IsoGraph switch pair
##                              (structure_switch_pairs + transcript_polarity, GENCODE v47),
##                              flag GTEx concordance + independent BrainSeq replication.
##
## Deterministic (name-resolved GWAS streaming + annotation joins; no RNG). Runs over all
## analysis dirs under real_data/coloc/_m that have a coloc_colocalized_genes.tsv, then the
## cross-trait rollups. Usage: sbatch real_data/coloc/_h/04.coloc_direction.sh
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/coloc/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log "**** coloc signed direction ****"
python -m isograph_benchmark.real_data.coloc_direction "$@"

log "**** resolve isoform events (+ BrainSeq replication) ****"
python -m isograph_benchmark.real_data.coloc_isoform_events "$@"

conda deactivate
log "**** coloc direction + isoform events done ****"
