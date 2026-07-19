#!/usr/bin/env bash
# Reproducibly build the per-gene deep-dive outputs and the genetic-anchoring figure set.
# Runs, in order:
#   1. gene_deep_dive.py  -> per-gene vignettes, panel, and the supplementary tables
#      (deep_dive_events / deep_dive_rbp / deep_dive_exon_clinical) for every colocalized gene;
#   2. extracts the SNCA switch-pair transcript exons from GENCODE (input to figure panel A);
#   3. renders the four figures (genetic-anchoring main + switch-consequence / RBP-regulon /
#      clinical-consequence supplements).
# Light enough to run interactively or on a login node. Usage (from repo root):
#   bash real_data/_h/build_deep_dive.sh
# Override references via env: GENCODE_GTF=/path/to.gtf
set -euo pipefail
log() { echo "$(date '+%H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${PWD}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: run from the repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"

GENCODE_GTF="${GENCODE_GTF:-/ocean/projects/bio260021p/shared/resources/genomes/human/gencode-v47/gtf/gencode.v47.annotation.gtf}"
RSCRIPT="${RSCRIPT:-/ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript}"
DD_DIR="real_data/_m/deep_dive"

# ---- env init (non-interactive shells need module + conda hook sourced explicitly) ----
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge 2>/dev/null || true
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

# ---- 1. per-gene deep-dive + supplementary tables ----
log "gene_deep_dive.py (vignettes + panel + supplementary tables)"
python -m isograph_benchmark.real_data.gene_deep_dive "$@"

# ---- 2. SNCA switch-pair transcript exons for figure panel A ----
log "extracting SNCA transcript exons from ${GENCODE_GTF##*/}"
mkdir -p "${DD_DIR}"
{
    printf "transcript\tfeature\tstart\tend\tstrand\n"
    # ENST…508895 / …618500 = the alternative-first-exon switch pair; …336904 = canonical
    for tx in ENST00000508895 ENST00000618500 ENST00000336904; do
        grep -F "${tx}" "${GENCODE_GTF}" \
          | awk -F'\t' -v tx="${tx}" '$3=="exon"{print tx"\t"$3"\t"$4"\t"$5"\t"$7}'
    done
} > "${DD_DIR}/snca_transcript_exons.tsv"

# ---- 3. figures ----
log "rendering figures"
"${RSCRIPT}" real_data/_h/genetic_anchoring_figure.R
"${RSCRIPT}" real_data/_h/switch_consequence_figure.R
"${RSCRIPT}" real_data/_h/rbp_regulon_figure.R
"${RSCRIPT}" real_data/_h/clinical_consequence_figure.R

conda deactivate 2>/dev/null || true
log "done -> ${DD_DIR} + real_data/_m/figures"
