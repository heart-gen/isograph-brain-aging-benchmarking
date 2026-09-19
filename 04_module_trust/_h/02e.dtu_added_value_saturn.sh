#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=dtu-added-value-saturn
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-30
#SBATCH --time=02:00:00
#SBATCH --output=04_module_trust/_m/logs/%x-%A_%a.log

## Module context and held-out DTU evidence (PI review item 12a), step 1 of 3: satuRn age-DTU on each
## stage-04 split half, plus per-gene technical covariates and the abundance co-expression
## neighbours the comparator needs. The halves are the ones 01a fitted IsoGraph on
## (stability._split_indices), so this waits on 01a only for the partitions the next step reads.
##
## Array: 6 split-half regions x 5 seeds; each task runs both halves (~1-2 min of satuRn each).
## Skips halves already written; pass --force to redo. Bridges memory is --cpus-per-task x
## 2000MB (16G); do NOT pass --mem. Calls the env interpreter directly.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p 04_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

REGIONS=(brainseq:caudate brainseq:hippocampus brainseq:dlpfc
         gtex:caudate_basal_ganglia gtex:hippocampus gtex:frontal_cortex_ba9)
IDX=$((${SLURM_ARRAY_TASK_ID:-1} - 1))
CR="${REGIONS[$((IDX / 5))]}"
SEED=$((IDX % 5))
log "module_dtu_added_value saturn --cohort ${CR%%:*} --region ${CR##*:} --seed ${SEED}"
"${PY}" -u -m isograph_benchmark.real_data.module_dtu_added_value saturn \
    --cohort "${CR%%:*}" --region "${CR##*:}" --seed "${SEED}" \
    --cores "${SLURM_CPUS_PER_TASK:-8}" "$@"
log "**** complete ****"
