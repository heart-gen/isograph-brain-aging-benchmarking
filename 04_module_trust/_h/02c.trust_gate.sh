#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=trust-gate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-12
#SBATCH --time=02:00:00
#SBATCH --output=04_module_trust/_m/logs/trust-gate-%A_%a.log

## Q1 trust gate (module_trust `stability` subcommand): for each production module, the mean
## split-half co-assignment density of its genes against a size-matched permutation null;
## modules with BH FDR < 0.05 are "trusted".
##
## Every downstream trust arm reads the trusted sets this writes -- cross-cohort replication
## (_h/03c, 03d, 03e), eigengene projection (_h/03f) and complementarity (_h/03g) -- yet it had no
## committed launcher and was run from an uncommitted funnel script. Run it after the
## split-half partitions (_h/01a-01b) and the production fits it validates:
##
##   sbatch 04_module_trust/_h/02c.trust_gate.sh
##
## Array = the six trust-funnel regions x {isograph, wgcna}. No model fits.
## Writes 04_module_trust/_m/stability/module_trust/module_stability__<cohort>__<region>__<method>.parquet
##
## Bridges memory is --cpus-per-task x 2000MB, so 8 cpus = 16G; do NOT pass --mem.
## Calls the env interpreter directly rather than `module load` + `conda activate`, so an
## array task cannot die on "module: command not found".
## Extra arguments are forwarded, e.g. `--n-perm 1000 --fdr 0.05`.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 04_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
REGIONS=(
    brainseq:caudate brainseq:dlpfc brainseq:hippocampus
    gtex:caudate_basal_ganglia gtex:frontal_cortex_ba9 gtex:hippocampus
)
METHODS=(isograph wgcna)
IDX=$((${SLURM_ARRAY_TASK_ID:-1} - 1))
PAIR="${REGIONS[$((IDX / 2))]}"
METHOD="${METHODS[$((IDX % 2))]}"
COHORT="${PAIR%%:*}"
REGION="${PAIR##*:}"

log "**** Q1 trust gate: ${COHORT}/${REGION}/${METHOD} ****"
"${PY}" -u -m isograph_benchmark.real_data.module_trust stability \
    --cohort "${COHORT}" --region "${REGION}" --method "${METHOD}" "$@"
log "**** Complete ****"
