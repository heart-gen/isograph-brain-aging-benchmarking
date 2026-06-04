#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=isograph-scale-realistic
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=64
#SBATCH --array=1-90%10
#SBATCH --time=48:00:00
#SBATCH --output=benchmark/01_synthetic/_m/logs/%x-%A_%a.log

# Runs scale_realistic scenario: 16,000 genes / 300 samples at 15 seeds.
# Array upper bound matches the number of batches in batches_scale_realistic.tsv.
# Generate the batch file first:
#   python -m isograph_benchmark.benchmark.make_batches \
#     --scenario-filter scale_realistic \
#     --out benchmark/01_synthetic/_m/batches_scale_realistic.tsv

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID} array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: ${SLURM_NODENAME:-local}"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

BATCH_FILE="${BATCH_FILE:-benchmark/01_synthetic/_m/batches_scale_realistic.tsv}"
BATCH_INDEX="${SLURM_ARRAY_TASK_ID:-1}"

log_message "**** Running scale_realistic batch ${BATCH_INDEX} ****"
python -m isograph_benchmark.benchmark.run_batch \
  --batch-file "${BATCH_FILE}" \
  --batch-index "${BATCH_INDEX}"

if [ $? -ne 0 ]; then
    log_message "Error: Python execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
