#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=within-cohort
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --array=1-6
#SBATCH --time=01:00:00
#SBATCH --output=03_module_trust/_m/logs/within-cohort-%A_%a.log

## Within-cohort split-half agreement: Q3 aging-sign concordance and Q2 driver
## reproducibility (module_trust `within` subcommand).
##
## Like the pooled arm before _h/15, this subcommand had no committed launcher, so nothing
## re-ran it when its inputs moved: the IsoGraph modules_meta it reads was built on the
## resolution-2.0 split halves and never rebuilt after the 2026-09-09 re-fit, leaving the
## Q2 driver statistics joined to a different partition. Run it after `meta` (_h/04):
##
##   j=$(sbatch --parsable --export=ALL,COHORT=brainseq,REGION=caudate 03_module_trust/_h/04.module_meta.sh)
##   sbatch --dependency=afterok:${j} 03_module_trust/_h/17.within_cohort.sh
##
## Array covers the six trust-funnel regions. METHOD env selects isograph (default) or wgcna.
## Writes 03_module_trust/_m/stability/module_trust/within_cohort__<cohort>__<region>__<method>.parquet
##
## Calls the env interpreter directly rather than `module load` + `conda activate`, so an
## array task cannot die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 03_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
METHOD="${METHOD:-isograph}"
REGIONS=(
    brainseq:caudate brainseq:dlpfc brainseq:hippocampus
    gtex:caudate_basal_ganglia gtex:frontal_cortex_ba9 gtex:hippocampus
)
PAIR="${REGIONS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}"
COHORT="${PAIR%%:*}"
REGION="${PAIR##*:}"

log "**** within-cohort: ${COHORT}/${REGION}/${METHOD} ****"
"${PY}" -u -m isograph_benchmark.real_data.module_trust within \
    --cohort "${COHORT}" --region "${REGION}" --method "${METHOD}" "$@"
log "**** Complete ****"
