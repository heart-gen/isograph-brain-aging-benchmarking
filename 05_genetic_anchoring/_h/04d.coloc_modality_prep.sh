#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-mod-prep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/coloc-mod-prep-%j.log
#
# Per-gene sQTL-vs-eQTL colocalization contrast, step 1: variant bridge + targets.
#
#   bridge — GTEx v8 WGS lookup (46.6M rows) -> per-chromosome variant_id/rsID parquet
#            under inputs/raw/gtex_v11/variant_bridge/. One-off; cached and gitignored.
#            GTEx v11 ships no lookup of its own; v8 positions/alleles are identical
#            (both GRCh38) and cover 98.6% of v11 variants (measured on chr22).
#   prep   — (analysis, locus, gene) coloc targets from the CLPP layer's testable loci,
#            plus GTEx's grouped-permutation representative sQTL intron per gene.
#
# Streams rather than loads: memory stays flat at one chromosome. ~10 min.
#
# Gene-pool arm (switch | background | wgcna_switch | wgcna_multiplex) for the prep stage;
# the bridge is arm-independent and cached, so later arms only pay for prep:
#   sbatch --export=ALL,COLOC_MODALITY_ARM=background 05_genetic_anchoring/_h/04d.coloc_modality_prep.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 05_genetic_anchoring/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

ARM="${COLOC_MODALITY_ARM:-switch}"
log_message "**** coloc modality contrast: bridge + prep (arm ${ARM}) ****"
python -m isograph_benchmark.real_data.coloc_modality_contrast --stage bridge "$@"
python -m isograph_benchmark.real_data.coloc_modality_contrast --stage prep --arm "${ARM}"
conda deactivate
log_message "**** Complete ****"
