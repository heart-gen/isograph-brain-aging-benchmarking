#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=isograph-interpret-collect
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=benchmark/02_interpret/_m/logs/interpret-collect-%j.log

# Merge the per-shard partials written by step_1.sh (array) into the
# canonical synthetic_interpret_{results,module_metrics}.parquet, then run the
# bootstrap summary once over the full set. Submit with:
#   sbatch --dependency=afterok:<step_1_array_jobid> step_2.sh

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Synthetic interpretation: collect shards + summarize ****"
mkdir -p benchmark/02_interpret/_m/logs

module purge
module load anaconda3/2024.10-1

log_message "Activating isograph environment"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

python -m isograph_benchmark.benchmark.interpret_modules --collect "$@"

conda deactivate
log_message "**** Complete ****"
