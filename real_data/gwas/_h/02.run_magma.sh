#!/bin/bash
# Run MAGMA gene-set enrichment for IsoGraph and WGCNA brain modules.
#
# Prerequisites (run 01.prep_module_gene_sets.R first):
#   real_data/gwas/_m/gene_sets/{isograph_vae,wgcna_gene}_gene_sets.txt
#
# This script can run locally if MAGMA is installed, or via SLURM on Bridges-2.
# To submit on Bridges-2: sbatch 02.run_magma.sh
#
# Steps:
#   1. SNP annotation (skip if already done in organoid project — reuse annot files)
#   2. Gene analysis (skip if .genes.raw already computed)
#   3. Gene-set enrichment for SCZ, MDD, BP
#
#SBATCH --partition=RM-shared
#SBATCH --job-name=magma_isograph
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=4
#SBATCH --time=02:00:00
#SBATCH --output=logs/magma_isograph_%j.log
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

# ── MAGMA binary and reference paths ──────────────────────────────────────────
# On Bridges-2 HPC:
MAGMA_HPC=/ocean/projects/bio250020p/shared/opt/magma-v1.10/magma
GENE_LOC_HPC=/ocean/projects/bio250020p/shared/opt/magma-v1.10/NCBI38.gene.loc
LD_REF_HPC=/ocean/projects/bio250020p/shared/opt/magma-v1.10/g1000_eur

# Reuse pre-computed gene-analysis results from organoid project (SCZ, MDD, BP)
# These were run on Bridges-2 with the same gene boundaries and LD reference.
ORGANOID_GENE_RESULTS=/ocean/projects/bio250020p/kbenjamin/cerebral_organoid_schizophrenia/\
differential_expression/gwas_overlap/_m/gene_analysis

MAGMA="${MAGMA:-${MAGMA_HPC}}"
GENE_LOC="${GENE_LOC:-${GENE_LOC_HPC}}"
LD_REF="${LD_REF:-${LD_REF_HPC}}"
GENE_RESULTS_DIR="${GENE_RESULTS_DIR:-${ORGANOID_GENE_RESULTS}}"

SETS_DIR="${PROJECT_ROOT}/real_data/gwas/_m/gene_sets"
OUT_DIR="${PROJECT_ROOT}/real_data/gwas/_m/results"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${OUT_DIR}" "${LOG_DIR}"

TRAITS=(scz mdd bp)
BACKENDS=(isograph_vae wgcna_gene)

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ── Step 3: Gene-set enrichment ────────────────────────────────────────────────
log "=== Gene-set enrichment ==="
for backend in "${BACKENDS[@]}"; do
    SET_FILE="${SETS_DIR}/${backend}_gene_sets.txt"
    if [[ ! -f "${SET_FILE}" ]]; then
        log "WARNING: gene set file not found: ${SET_FILE} — run 01.prep_module_gene_sets.R first"
        continue
    fi

    for trait in "${TRAITS[@]}"; do
        # Use noMHC results to avoid MHC inflation
        GENE_FILE="${GENE_RESULTS_DIR}/${trait}_noMHC.genes.raw"
        if [[ ! -f "${GENE_FILE}" ]]; then
            log "WARNING: gene results not found: ${GENE_FILE}"
            continue
        fi

        OUTFILE="${OUT_DIR}/${trait}_noMHC_${backend}"
        log "  ${trait} × ${backend} → ${OUTFILE}.gsa.out"

        "${MAGMA}" \
            --gene-results "${GENE_FILE}" \
            --set-annot    "${SET_FILE}" \
            --out          "${OUTFILE}" \
            2>&1 | tee "${LOG_DIR}/gsa_${trait}_${backend}.log"

        if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
            log "ERROR: MAGMA failed for ${trait} × ${backend}"; exit 1
        fi
    done
done

log "Gene-set enrichment complete."
log "Results in: ${OUT_DIR}"
