#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=gtex-celltype-composition
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-8
#SBATCH --time=03:00:00
#SBATCH --output=real_data/gtex/_m/logs/celltype-composition-%A_%a.log

# Cell-type composition confound test — GTEx aging replication arm. Unlike BrainSEQ (which
# reuses committed MuSiC fractions), GTEx has no prior deconvolution, so each task runs the
# full pipeline for one region with a defensibly matched Tran reference:
#   1. export-gtex   (isograph env) -> genes×samples bulk parquet for R.
#   2. 09.gtex_music_deconv.R (rnaseq env) -> MuSiC proportions + marker stats
#      (real_data/gtex/_m/composition/music-proportions-gtex-<region>.tsv).
#   3. celltype_composition fractions gtex-aging (isograph env) -> celltype_fractions.parquet
#      (sample_id-keyed) + marker-depletion cut.
#   4. incremental_association gtex-aging --composition -> repeats the de-confounded switch-
#      vs-abundance test with fractions as inference covariates, writing to
#      isograph_vae/incremental_association_composition/. Compare vs the canonical
#      incremental_association/ for the with-vs-without contrast.
#
# Only the 8 regions with a matched reference are deconvolved (striatum->NAc, cortex->DLPFC,
# AMY/sACC/HPC direct); cerebellum/hypothalamus/spinal cord/substantia nigra are skipped.
#
# Requires the canonical fit (modules.parquet + feature_scores.parquet) and the canonical
# incremental_association/ (05.incremental_association.sh) already on disk.

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/gtex/_m/logs

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

# Regions with a matched Tran reference (see GTEX_REF in celltype_composition.py).
REGIONS=(
    amygdala anterior_cingulate_cortex_ba24 cortex frontal_cortex_ba9
    hippocampus caudate_basal_ganglia putamen_basal_ganglia
    nucleus_accumbens_basal_ganglia
)
REGION="${REGIONS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}"
log_message "**** GTEx composition: ${REGION} ****"

ISO_ENV=/ocean/projects/bio260021p/shared/opt/envs/isograph
RNA_ENV=/ocean/projects/bio260021p/shared/opt/envs/rnaseq

module purge
module load anaconda3/2024.10-1

# 1. bulk export (isograph env)
conda activate "${ISO_ENV}"
python -m isograph_benchmark.real_data.celltype_composition export-gtex --region "${REGION}"
conda deactivate

# 2. MuSiC deconvolution (rnaseq env: MuSiC + DeconvoBuddies + arrow)
conda activate "${RNA_ENV}"
log_message "MuSiC deconvolution ${REGION}"
Rscript real_data/gtex/_h/09.gtex_music_deconv.R "${REGION}"
conda deactivate

# 3-4. fractions join + composition-adjusted incremental association (isograph env)
conda activate "${ISO_ENV}"
python -m isograph_benchmark.real_data.celltype_composition fractions gtex-aging --region "${REGION}"
python -m isograph_benchmark.real_data.incremental_association gtex-aging --region "${REGION}" --composition
conda deactivate

# After all array tasks finish, roll up BrainSEQ + GTEx into the shared summary:
#   python -m isograph_benchmark.real_data.celltype_composition meta
log_message "**** Complete ****"
