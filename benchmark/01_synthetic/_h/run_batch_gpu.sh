#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=GPU-shared
#SBATCH --job-name=isograph-synth-gpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:1
#SBATCH --array=1-197%50
#SBATCH --time=01:00:00
#SBATCH --output=benchmark/01_synthetic/_m/logs/%x-%A_%a.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Bridges-2 info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module load cuda
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

BATCH_FILE="${BATCH_FILE:-benchmark/01_synthetic/_m/batches_gpu.tsv}"
BATCH_INDEX="${SLURM_ARRAY_TASK_ID:-1}"

log_message "**** Subsetting benchmark run ****"
python -m isograph_benchmark.benchmark.run_batch \
  --batch-file "${BATCH_FILE}" \
  --batch-index "${BATCH_INDEX}"

if [ $? -ne 0 ]; then
    log_message "Error: Python execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
