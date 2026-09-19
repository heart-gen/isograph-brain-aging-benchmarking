#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-ld
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/coloc-ld-%j.log

## Coloc capstone, step 2 — per-locus LD matrices for GWAS SuSiE + a variant->rsID
## bridge for the QTL credible sets. Adapts the organoid pipeline (colocalization
## _h/04.locus_ld.sh) to the IsoGraph coloc prep (coloc_prep.py, step 1).
##
## For each locus in <coloc_m>/<analysis>/susie/loci_testable.tsv, compute the
## signed LD correlation matrix among the locus SNPs from the 1000G EUR Phase3 panel
## (hg19; same build as PGC3 + NCBI37.3 gene coords). plink2 '--r-unphased square
## ref-based bin4' writes a float32 matrix anchored to the REF allele; step 3
## re-signs the GWAS z to that REF allele before susie_rss.
##
## Usage (submit from repo root):
##   sbatch 05_genetic_anchoring/_h/02e.locus_ld.sh brainseq-sczd
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

ANALYSIS="${1:-brainseq-sczd}"
PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 05_genetic_anchoring/_m/logs

PANEL_DIR=/ocean/projects/bio250020p/shared/resources/ldsc/1000G_EUR_Phase3_plink
V8_LOOKUP=/ocean/projects/bio250020p/shared/resources/public-data/gtex-v8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz
MDIR="${PROJECT_ROOT}/05_genetic_anchoring/_m/coloc/${ANALYSIS}"
SUSIE_DIR="${MDIR}/susie"
LOCI="${SUSIE_DIR}/loci_testable.tsv"
[[ -f "${LOCI}" ]] || { echo "ERROR: ${LOCI} not found; run coloc_prep.py first."; exit 1; }

# Guard: some Bridges2 batch nodes start without Lmod initialised.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/genomics
command -v plink2 >/dev/null || { echo "ERROR: plink2 not on PATH."; exit 1; }

## --- QTL credible-set variant_id (b38) -> rsID bridge -----------------------
## GTEx v11 SuSiE_summary carries no rsID; the GWAS/LD side is rsID-keyed. Filter
## the GTEx v8 WGS lookup (b38 variant_id -> dbSNP rsID) to the CS variants of BOTH
## sQTL and eQTL. Positions/alleles are identical across v8/v11 (both GRCh38).
MAP="${MDIR}/variant_rsid_map.tsv"
log "Building variant_id -> rsID map (sQTL + eQTL CS variants)"
# variant_id is column 4 of qtl_credible_sets.tsv (gene gene_name phenotype_id variant_id ...)
tail -n +2 "${MDIR}/qtl_credible_sets.tsv" | cut -f4 | sort -u > "${MDIR}/_cs_variant_ids.txt"
zcat "${V8_LOOKUP}" | awk -F'\t' 'NR==FNR{want[$1]=1; next} FNR==1{next} ($1 in want) && $7!="." && $7!="" {print $1"\t"$7}' \
    "${MDIR}/_cs_variant_ids.txt" - > "${MAP}"
log "  mapped $(wc -l < "${MAP}") / $(wc -l < "${MDIR}/_cs_variant_ids.txt") CS variants to rsID"

## --- per-locus LD matrices ---------------------------------------------------
# loci_testable.tsv columns: LOCUS_ID chr start stop genes n_gwas_snp ...
# A matrix on disk is reused only if it belongs to THIS locus: every variant in its .vars is
# in the current SNP list and the file holds exactly n^2 float32s. Locus ids are positional
# (locusNN_chrC), so a re-prep with a different gene set re-uses the same ids for different
# regions; a bare existence test silently fed coloc stale LD (2026-09-17: 229 of 1,146
# aging loci). The size test also catches a matrix truncated by an OOM kill.
ld_valid() {
    local bin="$1.unphased.vcor1.bin" vars="$1.unphased.vcor1.bin.vars" n
    [[ -s "${bin}" && -s "${vars}" ]] || return 1
    n=$(wc -l < "${vars}")
    [[ $(stat -c %s "${bin}") -eq $(( n * n * 4 )) ]] || return 1
    awk 'NR==FNR{s[$1]; next} !($1 in s){exit 1}' "$2" "${vars}"
}
# plink2 must know the cgroup budget, not the node's RAM (a 54k-SNP square is 11.7 GB).
PLINK_MEM=$(( ${SLURM_CPUS_PER_TASK:-4} * 2000 * 9 / 10 ))
tail -n +2 "${LOCI}" | while IFS=$'\t' read -r LID CHR START STOP REST; do
    SNPS="${SUSIE_DIR}/${LID}.snps.txt"
    OUT="${SUSIE_DIR}/${LID}"
    [[ -s "${SNPS}" ]] || { log "skip ${LID} (no snp list)"; continue; }
    if ld_valid "${OUT}" "${SNPS}"; then continue; fi
    rm -f "${OUT}.unphased.vcor1.bin" "${OUT}.unphased.vcor1.bin.vars"
    log "${LID} (chr${CHR}, $(wc -l < "${SNPS}") SNPs)"
    if ! plink2 --bfile "${PANEL_DIR}/1000G.EUR.QC.${CHR}" \
                --extract "${SNPS}" \
                --r-unphased square ref-based bin4 \
                --threads "${SLURM_CPUS_PER_TASK:-4}" --memory "${PLINK_MEM}" \
                --out "${OUT}" > "${OUT}.plink.out" 2>&1; then
        tail -5 "${OUT}.plink.out"
        rm -f "${OUT}.unphased.vcor1.bin" "${OUT}.unphased.vcor1.bin.vars"
        log "ERROR: plink2 failed on ${LID}"; exit 1
    fi
    grep -E "variants remaining|Matrix written" "${OUT}.plink.out" || true
    ld_valid "${OUT}" "${SNPS}" || { log "ERROR: ${LID} matrix fails validation"; exit 1; }
done

n_ld=$(ls -1 "${SUSIE_DIR}"/*.unphased.vcor1.bin 2>/dev/null | wc -l)
log "LD matrices written: ${n_ld}"
conda deactivate
log "**** coloc LD step done ****"
