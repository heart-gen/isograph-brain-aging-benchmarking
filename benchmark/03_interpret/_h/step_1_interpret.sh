#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=isograph-interpret
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=benchmark/03_interpret/_m/logs/interpret-%j.log

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Synthetic module interpretation benchmark ****"
mkdir -p benchmark/03_interpret/_m/logs

module purge
module load anaconda3/2024.10-1

log_message "Activating isograph environment"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

python -m isograph_benchmark.benchmark.interpret_modules "$@"

conda deactivate
log_message "**** Complete ****"
