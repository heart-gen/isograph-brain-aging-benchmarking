#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=build-deep-dive
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=08_integration/_m/logs/%x-%j.log
# Reproducibly build the per-gene deep-dive panel and the genetic-anchoring figure set.
# Runs, in order:
#   1. gene_deep_dive.py --part panel -> per-gene vignettes, panel, and the supplementary tables
#      (deep_dive_rbp / deep_dive_exon_clinical / deep_dive_literature) for every colocalized
#      gene. The per-event table (deep_dive_events) is written in stage 05 by
#      05_genetic_anchoring/_h/05c.deep_dive_events.sh, because stage 06 reads it;
#   2. extracts the SNCA switch-pair transcript exons from GENCODE (input to figure panel A);
#   3. renders the four figures (genetic-anchoring main + switch-consequence / RBP-regulon /
#      clinical-consequence supplements).
# Needs stage 05 (coloc events + direction), stage 06 (clinical consequence constraint/exons) and
# stage 07 (rbp_switch_calls, rbp_regulon). Usage (from repo root):
#   sbatch 08_integration/_h/01a.build_deep_dive.sh
# Light enough to run with `bash` on a login node too. Override references via env:
# GENCODE_GTF=/path/to.gtf
set -euo pipefail
log() { echo "$(date '+%H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: run from the repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 08_integration/_m/logs

GENCODE_GTF="${GENCODE_GTF:-/ocean/projects/bio260021p/shared/resources/genomes/human/gencode-v47/gtf/gencode.v47.annotation.gtf}"
RSCRIPT="${RSCRIPT:-/ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript}"
DD_DIR="08_integration/_m/deep_dive"

# ---- env init (non-interactive shells need module + conda hook sourced explicitly) ----
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge 2>/dev/null || true
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

# ---- 1. per-gene deep-dive panel + supplementary tables ----
log "gene_deep_dive.py --part panel (vignettes + panel + supplementary tables)"
python -m isograph_benchmark.real_data.gene_deep_dive --part panel "$@"

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
"${RSCRIPT}" manuscript/_h/genetic_anchoring_figure.R
"${RSCRIPT}" manuscript/_h/switch_consequence_figure.R
"${RSCRIPT}" manuscript/_h/rbp_regulon_figure.R
"${RSCRIPT}" manuscript/_h/clinical_consequence_figure.R

conda deactivate 2>/dev/null || true
log "done -> ${DD_DIR} + manuscript/_m/figures"
