#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=dl-clinical
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/dl-clinical-%j.log

## Download the light-scope clinical-consequence reference data into inputs/raw/clinical/
## (gitignored; regenerable from this script):
##   - ClinVar GRCh38 VCF (Pathogenic/Likely_pathogenic density over switched exons)
##   - gnomAD v4.1 constraint metrics TSV (LOEUF gene-level anchor)
## Small enough for a login/transfer node too; run as `bash 06_switch_mechanism/_h/08.download_clinical.sh`
## or submit with sbatch. Safe to re-run (skips files already present).
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
DEST="inputs/raw/clinical"
mkdir -p "${DEST}" 02_module_discovery/brainseq/_m/logs

CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
GNOMAD_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"

fetch() {  # url dest
    local url="$1" out="$2"
    if [[ -s "${out}" ]]; then log "exists, skip: ${out}"; return; fi
    log "downloading ${url}"
    if command -v wget >/dev/null 2>&1; then wget -q -O "${out}.part" "${url}";
    else curl -fsSL -o "${out}.part" "${url}"; fi
    mv "${out}.part" "${out}"
    log "-> ${out} ($(du -h "${out}" | cut -f1))"
}

fetch "${CLINVAR_URL}"  "${DEST}/clinvar.vcf.gz"
fetch "${CLINVAR_URL}.tbi" "${DEST}/clinvar.vcf.gz.tbi" || log "note: .tbi optional (not required by the parser)"
fetch "${GNOMAD_URL}"   "${DEST}/gnomad.v4.1.constraint_metrics.tsv"
log "**** clinical reference data ready in ${DEST} ****"
