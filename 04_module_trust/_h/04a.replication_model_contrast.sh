#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=replication-model-contrast
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --output=04_module_trust/_m/logs/replication-model-contrast-%j.log

## Linear-vs-spline decomposition of the Q3 concordance count, both methods
## (module_trust `replication-model-contrast` subcommand): splits the count into sign
## agreement and both-cohort significance, so a reader can see which conjunct the age model
## moves.
##
## Had no committed launcher. It reads the per-pair tables of BOTH age models and re-fits
## nothing, so run it after the linear and spline arrays of _h/03c:
##   l=$(sbatch --parsable 04_module_trust/_h/03c.module_trust_replication.sh)
##   s=$(sbatch --parsable 04_module_trust/_h/03c.module_trust_replication.sh --model spline)
##   sbatch --dependency=afterok:${l}:${s} 04_module_trust/_h/04a.replication_model_contrast.sh
##
## Writes 04_module_trust/_m/stability/module_trust/replication_model_contrast__<method>.parquet
## Calls the env interpreter directly rather than `module load` + `conda activate`.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 04_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

for METHOD in isograph wgcna; do
    log "replication-model-contrast --method ${METHOD}"
    "${PY}" -u -m isograph_benchmark.real_data.module_trust replication-model-contrast \
        --method "${METHOD}"
done
log "**** Complete ****"
