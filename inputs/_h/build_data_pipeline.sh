#!/usr/bin/env bash
# =============================================================================
# IsoGraph Brain Aging — Data Preparation Pipeline
# =============================================================================
#
# This script documents and reproduces the full path from raw data sources to
# the IsoGraph dataset bundles used in all analyses.  Run it from the
# repository root:
#
#   cd <repo-root>
#   bash inputs/_h/build_data_pipeline.sh [--skip-raw] [--skip-pcs]
#
# PIPELINE OVERVIEW
# -----------------
#   Step 0  Copy raw data          (requires separate download; skip with --skip-raw)
#   Step 1  Convert to Parquet     (local; requires inputs/raw/)
#   Step 2  Compute SNP PCs        (requires genotypes (controlled access); skip with --skip-pcs)
#   Step 3  Build IsoGraph bundles (local; requires inputs/processed/)
#
# OUTPUTS (at completion)
# -----------------------
#   inputs/raw/          Raw TSV/GCT copies (gitignored; Zenodo-archived)
#   inputs/processed/    Per-region Parquet matrices and metadata
#   inputs/bundles/      IsoGraph dataset bundles:
#     brainseq_v1/
#       caudate/         BrainSEQ Phase 3 caudate — Control only (n≈279)
#       hippocampus/     BrainSEQ Phase 2 hippocampus — Control only
#       dlpfc/           BrainSEQ Phase 2 DLPFC — Control only
#     brainseq_sczd/
#       caudate/         BrainSEQ Phase 3 caudate — Control + SCZD (n≈439)
#     gtex_v11_brain/
#       amygdala/ … substantia_nigra/   GTEx v11 — 13 brain regions
#
# REQUIREMENTS
# ------------
#   Python  >= 3.11  (with isograph-brain-aging-benchmarking env active)
#   R       >= 4.3   (WGCNA, arrow, dplyr — for downstream analyses)
#   plink2  >= 2.00  (for Step 2; path: ~/.local/bin/plink2)
#
# REPRODUCIBILITY NOTE
# --------------------
# Raw data files are intentionally excluded from git (.gitignore).
# They will be archived on Zenodo with DOI <pending>. SNP PCs will be included.
# Steps 1–3 are fully reproducible from the Zenodo deposit.
# =============================================================================

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}"

# ── Argument parsing ──────────────────────────────────────────────────────────
SKIP_RAW=0
SKIP_PCS=0
for arg in "$@"; do
    case "${arg}" in
        --skip-raw) SKIP_RAW=1 ;;
        --skip-pcs) SKIP_PCS=1 ;;
        --help|-h)
            awk '/^set -euo pipefail/{exit} NR>1 && /^#/{sub(/^# ?/,""); print}' "${BASH_SOURCE[0]}"
            exit 0 ;;
        *) echo "Unknown argument: ${arg}"; exit 1 ;;
    esac
done

log() { printf '\n\033[1;34m[%s] %s\033[0m\n' "$(date '+%H:%M:%S')" "$*"; }
ok()  { printf '  \033[1;32m✓\033[0m %s\n' "$*"; }
warn(){ printf '  \033[1;33m⚠\033[0m %s\n' "$*"; }
die() { printf '  \033[1;31m✗\033[0m %s\n' "$*"; exit 1; }

check_python() { python3 -c "import isograph_benchmark" 2>/dev/null || die "isograph_benchmark not importable. Activate the project env first."; }
check_dir()    { [[ -d "$1" ]] && ok "$1" || warn "Not found: $1"; }

# =============================================================================
# STEP 0 — Copy raw data from HPC storage
# =============================================================================
# Source paths (defined in configs/data_sources.yaml):
#
# Key raw files produced:
#   inputs/raw/brainseq/counts/{caudate,hippocampus,dlpfc}/tx-counts.tsv
#   inputs/raw/brainseq/counts/{caudate,hippocampus,dlpfc}/gene-counts.tsv
#   inputs/raw/brainseq/counts/{caudate,hippocampus,dlpfc}/psi-events.tsv.gz
#   inputs/raw/brainseq/metadata/libd_rnaseq_metadata.tab
#   inputs/raw/brainseq/annotations/transcript-annotation.tsv
#   inputs/raw/gtex_v11/counts/GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_tpm.txt.gz
#   inputs/raw/gtex_v11/counts/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_reads.gct.gz
#   inputs/raw/gtex_v11/counts/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz
#   inputs/raw/gtex_v11/metadata/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt
#   inputs/raw/gtex_v11/metadata/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt
#   inputs/raw/gtex_v11/metadata/GTEx_Analysis_v11_Sample_Tissue_Changes_From_v8.txt
#   inputs/raw/gtex_v11/metadata_v8/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz
#   reports/raw_copy_manifest.parquet  (SHA-256 checksums for every copied file)
# =============================================================================
log "STEP 0 — Copy raw data"
if [[ "${SKIP_RAW}" -eq 1 ]]; then
    warn "Skipping (--skip-raw). Assuming inputs/raw/ is already populated."
    check_dir "inputs/raw/brainseq/counts/caudate"
    check_dir "inputs/raw/gtex_v11/counts"
else
    check_python
    echo "  Copying raw files from HPC storage (see configs/data_sources.yaml for paths)..."
    python3 -m isograph_benchmark.inputs.copy_raw
    ok "Raw copy complete. Manifest: reports/raw_copy_manifest.parquet"
fi

# =============================================================================
# STEP 1 — Convert raw files to Parquet
# =============================================================================
# Reads from inputs/raw/ and writes compressed Parquet to inputs/processed/.
#
# BrainSEQ regions processed: caudate, hippocampus, dlpfc
#   tx_counts.parquet      Transcript-level Salmon counts (wide: tx × samples)
#   gene_counts.parquet    Gene-level counts (wide: gene × samples)
#   psi_events.parquet     Percent-spliced-in events
#
# BrainSEQ metadata:
#   libd_rnaseq_metadata.parquet    Master sample manifest (all regions, all Dx)
#   {region}_rnaseq_metrics.parquet Per-sample QC metrics (RIN, mapping rate, etc.)
#
# GTEx v11 — per region (13 brain regions):
#   transcript_tpm.parquet    RSEM transcript TPM (wide: transcript × samples)
#   gene_reads.parquet        RNASeQC raw gene read counts
#   gene_tpm.parquet          RNASeQC gene TPM
#   sample_attributes.parquet Sample QC + subject phenotype metadata
#                             NOTE: AGE is resolved from v8 exact ages when
#                             available and falls back to v11 age-band midpoints.
# =============================================================================
log "STEP 1 — Convert raw files to Parquet"
check_python
echo "  Converting BrainSEQ TSVs and GTEx GCTs to Parquet..."
python3 -m isograph_benchmark.inputs.build_parquet
ok "Processed Parquet files in inputs/processed/"

# =============================================================================
# STEP 2 — Compute BrainSEQ SNP principal components
# =============================================================================
# Genotype source: TOPMed-imputed per-chromosome pgen files (controlled access).
#
# Steps inside compute_snp_pcs.sh:
#   1. LD prune each autosome (MAF>=0.05, geno<=0.05, HWE p>1e-6, r²<0.2)
#   2. Merge pruned variants across chr1-chr22
#   3. Compute 10 PCs with plink2 approximate algorithm
#
# Output (written from inside the _h/ directory):
#   inputs/processed/brainseq/genetic_similarity/_m/TOPMed_LIBD.eigenvec
#   inputs/processed/brainseq/genetic_similarity/_m/TOPMed_LIBD.eigenval
#
# The .eigenvec file is merged into the sample table at bundle-build time
# (build_bundles.py: _load_snp_pcs).  Samples without genotype data receive
# NaN PCs and are retained.
#
# Runtime: ~20-30 min on 16 cores.
# =============================================================================
log "STEP 2 — Compute BrainSEQ SNP PCs"
PC_OUT="inputs/processed/brainseq/genetic_similarity/_m/TOPMed_LIBD.eigenvec"
if [[ "${SKIP_PCS}" -eq 1 ]]; then
    warn "Skipping (--skip-pcs). Samples without PCs will have NaN covariate columns."
    check_dir "inputs/processed/brainseq/genetic_similarity/_m" 2>/dev/null || true
elif [[ -f "${PC_OUT}" ]]; then
    ok "SNP PCs already computed: ${PC_OUT}"
else
    echo "  Running plink2 PCA (~20-30 min, 16 threads)..."
    pushd "inputs/processed/brainseq/genetic_similarity/_h" > /dev/null
    bash compute_snp_pcs.sh 16
    popd > /dev/null
    ok "SNP PCs: ${PC_OUT}"
fi

# =============================================================================
# STEP 3 — Build IsoGraph dataset bundles
# =============================================================================
# Reads from inputs/processed/ (and SNP PCs from Step 2) to produce
# self-contained IsoGraph bundles in inputs/bundles/.
#
# Each bundle directory contains:
#   manifest.json          Schema, provenance, and column descriptions
#   samples.parquet        Sample metadata table (sample_id + all covariates)
#   genes.parquet          Gene feature table (gene_id, genomic coords)
#   transcripts.parquet    Transcript feature table (transcript_id → gene_id)
#   gene_counts.npz        Gene count matrix [genes × samples]
#   transcript_counts.npz  Transcript count matrix [transcripts × samples]
#
# Bundles produced:
#
#   brainseq_v1/caudate      BrainSEQ Phase 3 caudate
#                            Dx=Control, dropped=f, Age>=18
#                            Covariates: Sex, MoD, RIN, mapping_rate, mito_rate,
#                                        SNP_PC1–PC5
#
#   brainseq_v1/hippocampus  BrainSEQ Phase 2 hippocampus (same filters)
#
#   brainseq_v1/dlpfc        BrainSEQ Phase 2 DLPFC (BSP2 RiboZeroGold only)
#
#   brainseq_sczd/caudate    BrainSEQ Phase 3 caudate — Control + SCZD
#                            Used for DRD2 isoform case study and disease
#                            trait associations (Dx as trait column).
#                            Dx=Control+SCZD, dropped=f, Age>=18
#
#   gtex_v11_brain/<region>  GTEx v11 brain regions (13 total)
#                            Covariates: SEX, SMRIN, SMTSISCH, SMMAPRT
#                            Trait: AGE (exact age; v8 phenotypes preferred)
# =============================================================================
log "STEP 3 — Build IsoGraph bundles"
check_python
echo "  Building BrainSEQ and GTEx bundles..."
python3 -m isograph_benchmark.inputs.build_bundles

echo ""
log "Pipeline complete."
echo ""
echo "  Bundle inventory:"
for bundle_dir in inputs/bundles/*/*; do
    [[ -f "${bundle_dir}/manifest.json" ]] || continue
    n=$(python3 -c "
import json, numpy as np
m = json.load(open('${bundle_dir}/manifest.json'))
import pandas as pd
s = pd.read_parquet('${bundle_dir}/samples.parquet', columns=['sample_id'])
print(len(s))
" 2>/dev/null || echo "?")
    printf "    %-50s  %s samples\n" "${bundle_dir}" "${n}"
done
echo ""
echo "  Next steps:"
echo "    BrainSEQ aging (controls):  python3 -m isograph_benchmark.real_data.run_models"
echo "    GTEx aging:                 bash real_data/gtex/_h/01.run_isograph.sh"
echo "    GTEx WGCNA comparison:      Rscript real_data/gtex/_h/02.wgcna_gene.R"
echo "    BrainSEQ SCZD WGCNA:        Rscript real_data/brainseq/_h/02.run_wgcna_sczd.R"
echo "    DRD2 case study:            Rscript real_data/brainseq/_h/03.drd2_case_study.R"
echo "    MAGMA gene sets (HPC):      Rscript real_data/gwas/_h/01.prep_module_gene_sets.R"
echo "    MAGMA enrichment (HPC):     sbatch real_data/gwas/_h/02.run_magma.sh"
