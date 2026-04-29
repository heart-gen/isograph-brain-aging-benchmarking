#!/usr/bin/env bash
# Compute SNP principal components for BrainSEQ from updated TOPMed-imputed
# genotypes (PLINK2 pgen/psam/pvar format, per chromosome).
#
# Steps:
#   1. LD prune each autosomal chromosome (MAF ≥ 0.05, geno ≤ 0.05, HWE p > 1e-6)
#   2. Merge pruned variants across chr1-chr22
#   3. Run PCA (10 PCs, approximate algorithm)
#
# Output: inputs/raw/brainseq/genotypes/brainseq_snp_pcs.eigenvec
#         inputs/raw/brainseq/genotypes/brainseq_snp_pcs.eigenval
#
# Usage:
#   cd <repo-root>
#   bash inputs/_h/compute_snp_pcs.sh [n_threads]
#
# Requirements: plink2 in PATH (~/.local/bin/plink2)
# Runtime: ~20-30 min on 16 cores

set -euo pipefail

GENO_DIR="/projects/b1213/resources/processed-data/genotypes/qtl/all_samples"
OUT_DIR="./"
TMP_DIR="${OUT_DIR}/tmp_pca"
THREADS="${1:-16}"
N_PCS=10

PLINK2="${HOME}/.local/bin/plink2"

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

echo "=== BrainSEQ SNP PCA ==="
echo "Genotype dir : ${GENO_DIR}"
echo "Output dir   : ${OUT_DIR}"
echo "Threads      : ${THREADS}"
echo ""

# ---------------------------------------------------------------------------
# Step 1: LD prune per autosomal chromosome
# ---------------------------------------------------------------------------
echo "[1/3] LD pruning chr1-chr22 (MAF>=0.05, geno<=0.05, HWE p>1e-6, r2<0.2)..."
PRUNE_IN="${TMP_DIR}/all.prune.in"
> "${PRUNE_IN}"

for chr in $(seq 1 22); do
    "${PLINK2}" \
        --pfile "${GENO_DIR}/chr${chr}" \
        --maf 0.05 \
        --geno 0.05 \
        --hwe 1e-6 midp \
        --indep-pairwise 200 100 0.2 \
        --threads "${THREADS}" \
        --out "${TMP_DIR}/chr${chr}_prune" \
        --silent
    cat "${TMP_DIR}/chr${chr}_prune.prune.in" >> "${PRUNE_IN}"
done

N_PRUNED=$(wc -l < "${PRUNE_IN}")
echo "  Pruned variant set: ${N_PRUNED} variants"

# ---------------------------------------------------------------------------
# Step 2: Merge all autosomes and extract pruned variants
# ---------------------------------------------------------------------------
echo "[2/3] Merging autosomes and extracting pruned variants..."
MERGE_LIST="${TMP_DIR}/pmerge.list"
> "${MERGE_LIST}"
for chr in $(seq 1 22); do
    echo "${GENO_DIR}/chr${chr}" >> "${MERGE_LIST}"
done

"${PLINK2}" \
    --pmerge-list "${MERGE_LIST}" pfile \
    --extract "${PRUNE_IN}" \
    --make-pgen \
    --threads "${THREADS}" \
    --out "${TMP_DIR}/merged_pruned"

# ---------------------------------------------------------------------------
# Step 3: PCA
# ---------------------------------------------------------------------------
echo "[3/3] Running PCA (${N_PCS} PCs, approx algorithm)..."
"${PLINK2}" \
    --pfile "${TMP_DIR}/merged_pruned" \
    --pca ${N_PCS} approx \
    --threads "${THREADS}" \
    --out "${OUT_DIR}/brainseq_snp_pcs"

echo ""
echo "Done."
echo "Eigenvectors : ${OUT_DIR}/brainseq_snp_pcs.eigenvec"
echo "Eigenvalues  : ${OUT_DIR}/brainseq_snp_pcs.eigenval"
echo ""
echo "Samples with PCs:"
tail -n +2 "${OUT_DIR}/brainseq_snp_pcs.eigenvec" | wc -l

echo "Cleaning up temporary files..."
rm -rf "${TMP_DIR}"
