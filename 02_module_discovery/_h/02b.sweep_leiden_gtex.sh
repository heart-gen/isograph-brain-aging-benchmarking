#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=gtex-leiden-sweep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16
#SBATCH --array=1-13
#SBATCH --time=06:00:00
#SBATCH --output=02_module_discovery/_m/logs/gtex-leiden-sweep-%A_%a.log

## Leiden resolution sweep on the 13 GTEx v11 brain fits, extending _h/02a (BrainSEQ only).
##
## Re-clusters the saved IsoGraph edges at each resolution WITHOUT refitting the VAE, and
## reports n_modules / giant fraction / NMI to the previous resolution / GO-enriched module
## count / AGE-associated module counts (linear and df=3 spline on the GTEx QC covariates).
## Writes 02_module_discovery/gtex/<region>/_m/isograph_vae/leiden_sweep_results.parquet.
##
## This is a DISCLOSED SENSITIVITY, not a re-selection: resolution 5.0 stays canonical on the
## >= 900-gene giant-module criterion (PI decision 2026-09-12), and the CLI refuses
## --write-best for gtex-aging so the production partitions cannot be overwritten.
## Grid defaults to the aging grid plus the production 5.0; override with --resolutions.
##
## Bridges memory is --cpus-per-task x 2000MB (32G here); do NOT pass --mem. Calls the env
## interpreter directly rather than `module load` + `conda activate`, so an array task cannot
## die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 02_module_discovery/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16

REGIONS=(
    amygdala anterior_cingulate_cortex_ba24 caudate_basal_ganglia
    cerebellar_hemisphere cerebellum cortex frontal_cortex_ba9
    hippocampus hypothalamus nucleus_accumbens_basal_ganglia
    putamen_basal_ganglia spinal_cord_cervical_c_1 substantia_nigra
)
REGION="${REGIONS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}"

log "**** GTEx Leiden resolution sweep: ${REGION} ****"
"${PY}" -u -m isograph_benchmark.real_data.sweep_leiden gtex-aging --region "${REGION}" "$@"
log "**** Complete ****"
