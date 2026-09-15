#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-mod-abf
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --array=0-12
#SBATCH --output=05_genetic_anchoring/_m/logs/coloc-mod-abf-%A_%a.log
#
# Per-gene sQTL-vs-eQTL colocalization contrast, step 2: coloc.abf, one brain tissue
# per array task. Extraction from GTEx v11 cis all-pairs is by parquet predicate
# pushdown (a chromosome's candidate genes read in ~1 s eQTL / ~4 s sQTL), so the
# 273 GB of brain all-pairs is never scanned.
#
# Memory on PSC is --cpus-per-task x 2000MB; 8 cpus = 16 GB, which holds one
# chromosome of both modalities plus the chromosome's variant bridge.
#
# Requires 04d.coloc_modality_prep.sh to have run (targets + variant bridge).

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
mkdir -p 05_genetic_anchoring/_m/logs

TISSUES=(Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia \
         Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 \
         Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia \
         Brain_Putamen_basal_ganglia Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra)
IDX="${SLURM_ARRAY_TASK_ID:-${1:-0}}"
TISSUE="${TISSUES[${IDX}]}"

# Guard: older wrappers lose array tasks to "module: command not found" when the Lmod
# init is absent in the array-task shell.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

# Gene pool. `switch` (default) is the IsoGraph arm and writes to the top-level dir;
# background / wgcna_switch / wgcna_multiplex write under arms/<arm>/. Loci and GWAS are
# identical across arms. Override per submission:
#   sbatch --export=ALL,COLOC_MODALITY_ARM=background 05_genetic_anchoring/_h/05a.coloc_modality_abf.sh
export COLOC_MODALITY_ARM="${COLOC_MODALITY_ARM:-switch}"

log_message "**** coloc.abf modality contrast: ${TISSUE} (task ${IDX}, arm ${COLOC_MODALITY_ARM}) ****"
Rscript 05_genetic_anchoring/_h/05a.coloc_modality_abf.R "${TISSUE}"
conda deactivate
log_message "**** Complete: ${TISSUE} ****"
