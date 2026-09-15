#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-interpret
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --output=03_module_characterization/_m/logs/interpret-%j.log

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** BrainSeq module interpretation ****"
mkdir -p 03_module_characterization/_m/logs

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

python -m isograph_benchmark.real_data.interpret_modules --scope brainseq "$@"

conda deactivate
log_message "**** Complete ****"
