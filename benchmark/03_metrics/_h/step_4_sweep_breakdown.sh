#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=isograph-sweep-breakdown
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=benchmark/03_metrics/_m/logs/sweep-breakdown-%j.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Step 4: Per-sweep-point metric breakdown ****"
echo "User: ${USER}"
echo "Host: ${HOSTNAME}"
mkdir -p benchmark/03_metrics/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1

log_message "Activating isograph environment"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "Re-grouping the long metric table by within-scenario sweep axes"
python -m isograph_benchmark.stats.sweep_breakdown

if [ $? -ne 0 ]; then
    log_message "ERROR: sweep_breakdown failed"
    exit 1
fi

conda deactivate
log_message "**** Step 4 complete ****"
