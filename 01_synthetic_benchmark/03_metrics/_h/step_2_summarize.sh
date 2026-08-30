#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=isograph-summarize
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=01_synthetic_benchmark/03_metrics/_m/logs/summarize-%j.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Step 2: Summarize benchmark metrics ****"
echo "User: ${USER}"
echo "Host: ${HOSTNAME}"
mkdir -p 01_synthetic_benchmark/03_metrics/_m/logs

module purge
module load anaconda3/2024.10-1

log_message "Activating isograph environment"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "Running bootstrap CI summary"
python -m isograph_benchmark.stats.summarize

if [ $? -ne 0 ]; then
    log_message "ERROR: summarize failed"
    exit 1
fi

conda deactivate
log_message "**** Step 2 complete ****"
