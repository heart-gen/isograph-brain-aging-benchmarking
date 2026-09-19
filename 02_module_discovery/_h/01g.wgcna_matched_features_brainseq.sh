#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-wgcna-feat
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=02_module_discovery/_m/logs/%x-%A_%a.log

set -euo pipefail

REGIONS=(caudate hippocampus dlpfc)
IDX=${SLURM_ARRAY_TASK_ID:-0}
REGION=${REGIONS[$IDX]}

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
mkdir -p 02_module_discovery/_m/logs
# Call the env directly instead of `module load` + `conda activate`: in a non-login shell
# `module` is undefined, and with `set -euo pipefail` an ARRAY TASK then dies before it
# reaches the analysis -- silently losing regions from the run. Putting the env's bin on
# PATH also supplies the Rscript that _check_rscript() looks for.
ENV_ROOT=/ocean/projects/bio260021p/shared/opt/envs/isograph
export PATH="${ENV_ROOT}/bin:${PATH}"
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
"${ENV_ROOT}/bin/python" -m isograph_benchmark.real_data.run_matched_wgcna brainseq-aging --region "${REGION}" "$@"
