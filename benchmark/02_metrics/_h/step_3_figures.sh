#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=isograph-figures
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=benchmark/02_metrics/_m/logs/figures-%j.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Step 3: Generate figures and tables (R/ggplot2) ****"
echo "User: ${USER}"
echo "Host: ${HOSTNAME}"
mkdir -p benchmark/02_metrics/_m/logs

module purge
module load anaconda3/2024.10-1

log_message "Activating R_env (ggplot2 + patchwork + arrow)"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "Running R figure script"
Rscript isograph_benchmark/figures/synthetic_benchmark.R

if [ $? -ne 0 ]; then
    log_message "ERROR: R figure generation failed"
    exit 1
fi

conda deactivate
log_message "**** Step 3 complete ****"
