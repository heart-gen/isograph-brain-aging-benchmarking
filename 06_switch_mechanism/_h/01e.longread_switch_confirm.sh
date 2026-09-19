#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=longread-confirm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=02:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/longread-confirm-%j.log

# A2: independent long-read confirmation of IsoGraph aging switch transcript pairs.
# Aguzzoli-Heberle 2024 (Nat Biotechnol) deep ONT DLPFC BA9/46, n=12, Bambu quant
# (Zenodo 10.5281/zenodo.8180677). Region-matched GTEx cortical aging switches
# (frontal_cortex_ba9 + cortex; BrainSEQ DLPFC aging is vacuous). Tier-1 = processed
# matrix only (no ~936 GB raw pull): isoform-level existence + switch-like usage.
#
# Submit from the repository root:
#   sbatch 06_switch_mechanism/_h/01e.longread_switch_confirm.sh

set -euo pipefail
umask 0027

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repository root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi

PYTHON="/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python"
if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: missing the isograph Python interpreter."
    exit 1
fi

export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p 06_switch_mechanism/_m/logs

log_message "**** A2 long-read switch confirmation (Bambu ONT DLPFC BA9/46) ****"

# Fetch the processed Bambu matrix (idempotent; ~92 MB) then confirm.
"${PYTHON}" -m isograph_benchmark.real_data.longread_switch_confirm fetch
"${PYTHON}" -m isograph_benchmark.real_data.longread_switch_confirm confirm --trait age "$@"

log_message "outputs: 06_switch_mechanism/_m/longread_switch_confirm/"
log_message "**** complete ****"
