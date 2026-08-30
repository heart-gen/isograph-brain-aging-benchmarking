#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=isograph-synth-multiplex
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --array=1-300%50
#SBATCH --time=02:00:00
#SBATCH --output=01_synthetic_benchmark/01_synthetic/_m/logs/%x-%A_%a.log
# Runs isograph_vae_multiplex and abundance_switch_mixed scenario jobs.
# Update --array upper bound to match the row count in batches_multiplex_vae.tsv
# before submitting (run make_multiplex_batches.sh first).

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
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

BATCH_FILE="${BATCH_FILE:-01_synthetic_benchmark/01_synthetic/_m/multiplex_batches_vae.tsv}"
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
