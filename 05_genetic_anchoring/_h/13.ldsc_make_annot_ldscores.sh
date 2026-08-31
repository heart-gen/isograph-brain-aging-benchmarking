#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ldsc-annot
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=03:00:00
#SBATCH --array=1-22
#SBATCH --output=05_genetic_anchoring/_m/logs/ldsc-annot-%A_%a.log

## S-LDSC step 2 — make_annot + LD scores for the IsoGraph sQTL/eQTL/cis switch
## annotations, one SLURM task per chromosome. Adapts ancestry-aging 14_gwas_ldsc
## steps 3-4 to PSC paths + the genomics env (ldsc is py3 there). Reads the BED
## files from ldsc_annot_prep.py (step 1).
##
## Usage: sbatch 05_genetic_anchoring/_h/13.ldsc_make_annot_ldscores.sh brainseq-sczd
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

ANALYSIS="${1:-brainseq-sczd}"
CHR="${SLURM_ARRAY_TASK_ID:-1}"
PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 05_genetic_anchoring/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
# make_annot.py -> pybedtools needs the bedtools `sortBed` binary on PATH.
module load bedtools/2.30.0
conda activate /ocean/projects/bio250020p/shared/opt/env/genomics

LDSC_DIR=/ocean/projects/bio250020p/shared/opt/ldsc
RES=/ocean/projects/bio250020p/shared/resources/ldsc
BIM_DIR="${RES}/1000G_EUR_Phase3_plink"
HAPMAP3_SNPS="${RES}/hm3_no_MHC.list.txt"
WRAP="05_genetic_anchoring/_h/ldsc_wrapper.py"

MDIR="${PROJECT_ROOT}/05_genetic_anchoring/_m/ldsc/${ANALYSIS}"
BED_DIR="${MDIR}/beds"
BIM="${BIM_DIR}/1000G.EUR.QC.${CHR}"
[[ -f "${BIM}.bim" ]] || { echo "ERROR: bim not found ${BIM}.bim"; exit 1; }

log "**** make_annot + ldscores: ${ANALYSIS} chr${CHR} ****"
for ANNOT in sqtl_switch eqtl_switch cis_switch; do
    BED="${BED_DIR}/${ANNOT}_hg19.bed"
    [[ -f "${BED}" ]] || { log "  ${ANNOT}: BED missing, skip"; continue; }
    OUTD="${MDIR}/ldscores/${ANNOT}"; mkdir -p "${OUTD}"
    PREFIX="${OUTD}/${ANNOT}.${CHR}"

    if [[ ! -f "${PREFIX}.annot.gz" ]]; then
        python "${WRAP}" "${LDSC_DIR}" make_annot.py \
            --bed-file "${BED}" --bimfile "${BIM}.bim" \
            --annot-file "${PREFIX}.annot.gz"
    fi
    if [[ ! -f "${PREFIX}.l2.ldscore.gz" ]]; then
        log "  ldscores: ${ANNOT} chr${CHR}"
        python "${WRAP}" "${LDSC_DIR}" ldsc.py \
            --l2 --bfile "${BIM}" --ld-wind-cm 1 \
            --annot "${PREFIX}.annot.gz" --thin-annot \
            --print-snps "${HAPMAP3_SNPS}" --out "${PREFIX}"
    fi
done
conda deactivate
log "**** chr${CHR} done ****"
