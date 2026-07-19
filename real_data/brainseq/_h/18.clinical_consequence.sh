#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=clinical-consequence
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --array=0-16
#SBATCH --output=real_data/brainseq/_m/logs/clinical-consequence-%A_%a.log

## Clinical-consequence of switched exons (clinical_consequence.py) across all 17 switch-layer
## regions: do the exons IsoGraph switches carry higher ClinVar P/LP density than the gene's
## constitutive exons (within-gene permutation), and are the switch genes gnomAD-constrained?
## Requires inputs/raw/clinical/ (run 17.download_clinical.sh first). Deterministic (--seed 13).
## One region per array task. Usage: sbatch real_data/brainseq/_h/18.clinical_consequence.sh
## After the array finishes, roll up with:
##   python -m isograph_benchmark.real_data.clinical_consequence_meta
set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/brainseq/_m/logs

# array index -> (tree, region): same 17-region switch layer as 15.switch_consequence.sh.
SPECS=(
    "brainseq caudate_sczd"
    "brainseq caudate"
    "brainseq hippocampus"
    "brainseq dlpfc"
    "gtex amygdala"
    "gtex anterior_cingulate_cortex_ba24"
    "gtex caudate_basal_ganglia"
    "gtex cerebellar_hemisphere"
    "gtex cerebellum"
    "gtex cortex"
    "gtex frontal_cortex_ba9"
    "gtex hippocampus"
    "gtex hypothalamus"
    "gtex nucleus_accumbens_basal_ganglia"
    "gtex putamen_basal_ganglia"
    "gtex spinal_cord_cervical_c_1"
    "gtex substantia_nigra"
)
read -r TREE REGION <<< "${SPECS[${SLURM_ARRAY_TASK_ID:-0}]}"
ARTIFACT="real_data/${TREE}/${REGION}/_m/isograph_vae"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** clinical consequence: ${TREE}/${REGION} ****"
python -m isograph_benchmark.real_data.clinical_consequence \
    --region "${TREE}_${REGION}" --artifact-dir "${ARTIFACT}" "$@"
conda deactivate
log_message "**** Complete ****"
