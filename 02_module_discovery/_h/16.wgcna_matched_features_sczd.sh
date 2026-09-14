#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=sczd-wgcna-feat
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=02_module_discovery/_m/logs/%x-%j.log

## Matched-feature WGCNA baselines (`wgcna_switch_only`, `wgcna_multiplex`) for the SCZD
## caudate store, the one store _h/10–11 never covered. Same recipe as the aging baselines:
## the IsoGraph feature matrix rebuilt from the bundle (unfiltered, exactly as
## run_models.run_brainseq_caudate_sczd feeds the fit), residualized on the BrainSEQ QC and
## genotype-PC covariates (Dx is never regressed), then WGCNA; diagnosis association downstream
## on Age + those covariates. Writes 02_module_discovery/brainseq/caudate_sczd/_m/wgcna_{switch_only,multiplex}/.
##
## Afterwards: enrichment (04/_h/03 task 1, `--method all`) and the baseline comparison (04/_h/14).
##
## Calls the env directly rather than `module load` + `conda activate` (a non-login shell has no
## `module`); the env's bin on PATH supplies the Rscript that _check_rscript() needs.
set -euo pipefail

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from the repo root."; exit 1; }
mkdir -p 02_module_discovery/_m/logs
ENV_ROOT=/ocean/projects/bio260021p/shared/opt/envs/isograph
export PATH="${ENV_ROOT}/bin:${PATH}"
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
echo "$(date '+%Y-%m-%d %H:%M:%S') - matched-feature WGCNA: brainseq-sczd"
"${ENV_ROOT}/bin/python" -m isograph_benchmark.real_data.run_matched_wgcna brainseq-sczd "$@"
echo "$(date '+%Y-%m-%d %H:%M:%S') - **** complete ****"
