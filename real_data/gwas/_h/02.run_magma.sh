#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=magma-isograph
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=real_data/gwas/_m/logs/magma-isograph-%j.log

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "ERROR: project root does not exist: ${PROJECT_ROOT}"
    exit 1
fi
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the isograph-brain-aging-benchmarking repo root or set ISOGRAPH_BENCHMARK_ROOT."
    echo "Current project root candidate: ${PROJECT_ROOT}"
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

MAGMA_HPC=/ocean/projects/bio250020p/shared/opt/magma-v1.10/magma
# The g1000_eur LD reference panel is build hg19/b37 (verified: bim chr1 max pos
# 249,240,539 -> hg19, and the rsID-keyed positions match the b37 coords embedded
# in the CLOZUK SCZ SNP ids, not the hg38 BP column). The gene location file MUST
# match that build, so use NCBI37.3 -- pairing it with NCBI38 mis-maps SNPs to genes.
GENE_LOC_HPC=/ocean/projects/bio250020p/shared/opt/magma-v1.10/NCBI37.3.gene.loc
LD_REF_HPC=/ocean/projects/bio250020p/shared/opt/magma-v1.10/g1000_eur

MAGMA="${MAGMA:-${MAGMA_HPC}}"
GENE_LOC="${GENE_LOC:-${GENE_LOC_HPC}}"
LD_REF="${LD_REF:-${LD_REF_HPC}}"
CONFIG="${CONFIG:-configs/gwas_magma.yaml}"

TRAITS=(scz mdd bp ad pd stroke lbd als)
# IsoGraph backend dir is selectable (default canonical isograph_vae). Set
# MAGMA_ISOGRAPH_BACKEND=isograph_vae_res5 to run the resolution-5.0 partition;
# the gene analysis (.genes.out/.genes.raw) is module-independent and reused,
# only the gene-set analysis (GSA) differs. Output files embed ${backend} so a
# non-canonical backend never clobbers the canonical results.
ISOGRAPH_BACKEND="${MAGMA_ISOGRAPH_BACKEND:-isograph_vae}"
BACKENDS=("${ISOGRAPH_BACKEND}" wgcna_gene)

SETS_DIR="${PROJECT_ROOT}/real_data/gwas/_m/gene_sets"
PVAL_DIR="${PROJECT_ROOT}/real_data/gwas/_m/tmp/magma_pvals"
REF_TMP_DIR="${PROJECT_ROOT}/real_data/gwas/_m/tmp/magma_ref"
GENE_RESULTS_DIR="${PROJECT_ROOT}/real_data/gwas/_m/gene_analysis"
OUT_DIR="${PROJECT_ROOT}/real_data/gwas/_m/results"
LOG_DIR="${PROJECT_ROOT}/real_data/gwas/_m/logs"
ANNOT_PREFIX="${REF_TMP_DIR}/g1000_eur_NCBI37"
SNP_LOC="${REF_TMP_DIR}/g1000_eur.snps.loc"

mkdir -p "${PVAL_DIR}" "${REF_TMP_DIR}" "${GENE_RESULTS_DIR}" "${OUT_DIR}" "${LOG_DIR}"

log_message "**** MAGMA job starts ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID:-local}"
echo "Job name: ${SLURM_JOB_NAME:-magma-isograph}"
echo "Node name: ${SLURM_NODENAME:-local}"
echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "Activating IsoGraph environment"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "Checking Python analysis dependencies"
python -c "import isograph_benchmark, numpy, pandas, scipy, patsy"

log_message "Preparing MAGMA SNP p-value inputs"
python -m isograph_benchmark.gwas.prepare_magma_inputs --config "${CONFIG}"

log_message "Preparing MAGMA SNP location file"
if [[ ! -f "${SNP_LOC}" ]]; then
    awk 'BEGIN {OFS="\t"} {print $2, $1, $4}' "${LD_REF}.bim" > "${SNP_LOC}"
fi

log_message "Building MAGMA gene annotation"
if [[ ! -f "${ANNOT_PREFIX}.genes.annot" ]]; then
    "${MAGMA}" \
        --annotate \
        --snp-loc "${SNP_LOC}" \
        --gene-loc "${GENE_LOC}" \
        --out "${ANNOT_PREFIX}" \
        2>&1 | tee "${LOG_DIR}/annotate.log"
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        log_message "ERROR: MAGMA annotation failed"
        exit 1
    fi
fi

log_message "Running MAGMA gene analysis"
for trait in "${TRAITS[@]}"; do
    PVAL_FILE="${PVAL_DIR}/${trait}_noMHC.tsv"
    GENE_PREFIX="${GENE_RESULTS_DIR}/${trait}_noMHC"

    if [[ ! -f "${PVAL_FILE}" ]]; then
        log_message "ERROR: missing p-value file ${PVAL_FILE}"
        exit 1
    fi

    if [[ -f "${GENE_PREFIX}.genes.raw" ]]; then
        log_message "  ${trait}: existing ${GENE_PREFIX}.genes.raw"
        continue
    fi

    log_message "  ${trait}: MAGMA gene analysis"
    # NB: do NOT pass --genes-only here. It suppresses the .genes.raw file (the
    # gene-gene correlation matrix), which the downstream gene-set analysis below
    # consumes via --gene-results. We need both .genes.out and .genes.raw.
    "${MAGMA}" \
        --bfile "${LD_REF}" \
        --pval "${PVAL_FILE}" ncol=N \
        --gene-annot "${ANNOT_PREFIX}.genes.annot" \
        --out "${GENE_PREFIX}" \
        2>&1 | tee "${LOG_DIR}/gene_${trait}.log"

    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        log_message "ERROR: MAGMA gene analysis failed for ${trait}"
        exit 1
    fi
done

log_message "Running MAGMA gene-set enrichment"
for backend in "${BACKENDS[@]}"; do
    SET_FILE="${SETS_DIR}/${backend}_gene_sets.txt"
    if [[ ! -f "${SET_FILE}" ]]; then
        log_message "ERROR: gene set file not found: ${SET_FILE}. Run 01.prep_module_gene_sets.R first."
        exit 1
    fi

    for trait in "${TRAITS[@]}"; do
        GENE_FILE="${GENE_RESULTS_DIR}/${trait}_noMHC.genes.raw"
        OUTFILE="${OUT_DIR}/${trait}_noMHC_${backend}"

        if [[ ! -f "${GENE_FILE}" ]]; then
            log_message "ERROR: gene results not found: ${GENE_FILE}"
            exit 1
        fi

        log_message "  ${trait} x ${backend}"
        "${MAGMA}" \
            --gene-results "${GENE_FILE}" \
            --set-annot "${SET_FILE}" \
            --out "${OUTFILE}" \
            2>&1 | tee "${LOG_DIR}/gsa_${trait}_${backend}.log"

        if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
            log_message "ERROR: MAGMA GSA failed for ${trait} x ${backend}"
            exit 1
        fi
    done
done

conda deactivate
log_message "MAGMA results in ${OUT_DIR}"
log_message "**** MAGMA job ends ****"
