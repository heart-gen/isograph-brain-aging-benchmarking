#!/usr/bin/env bash
## Long-read confirmation of the module *driver* switches (local; no scheduler).
##
## Scores the top-k driver transcripts of the age-associated GTEx cortical co-switch
## modules on the ONT DLPFC matrix, reusing the pairs already scored by 01e
## (`longread_switch_confirm`) and the abundance-matched null of 02b. Produces the citable
## "driver switches detected / switch-like in long read" numbers with their null.
##
## Run from the repository root:
##   bash 06_switch_mechanism/_h/02f.longread_driver_confirm.sh            # k = 5 (primary)
##   bash 06_switch_mechanism/_h/02f.longread_driver_confirm.sh --k 10     # one k only
##   SENSITIVITY=0 bash 06_switch_mechanism/_h/02f.longread_driver_confirm.sh  # skip k = 1, 10
## Extra arguments are forwarded to the CLI (e.g. --alpha, --n-draws, --regions).
##
## Python: $PYTHON if set, else ~/.venvs/isograph, else `uv run` with
## env/requirements-analysis.txt (the CLI needs only pandas/pyarrow/numpy/scipy/pyhere).
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${PWD}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: run from the repository root."; exit 1; }
mkdir -p 06_switch_mechanism/_m/logs
LOG="06_switch_mechanism/_m/logs/longread-driver-$(date '+%Y%m%d-%H%M%S').log"
exec > >(tee -a "${LOG}") 2>&1

export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

PY=()
if [[ -n "${PYTHON:-}" ]]; then
    PY=("${PYTHON}")
elif [[ -x "${HOME}/.venvs/isograph/bin/python" ]]; then
    PY=("${HOME}/.venvs/isograph/bin/python")
elif command -v uv >/dev/null; then
    PY=(uv run --quiet --no-project --with-requirements env/requirements-analysis.txt python)
else
    echo "ERROR: no Python found; set PYTHON or install uv."; exit 1
fi
"${PY[@]}" -c "import pandas, pyarrow, scipy" || { echo "ERROR: ${PY[*]} lacks pandas/pyarrow/scipy."; exit 1; }

PAIRS=06_switch_mechanism/_m/longread_switch_confirm/pair_confirmation.parquet
[[ -s "${PAIRS}" ]] || { echo "ERROR: missing ${PAIRS}; run 01e first."; exit 1; }

MOD=isograph_benchmark.real_data.longread_driver_confirm
log "**** long-read driver confirmation ($*) ****"
"${PY[@]}" -m "${MOD}" "$@"

## k sensitivity unless a --k was given or SENSITIVITY=0
if [[ "${SENSITIVITY:-1}" == "1" && " $* " != *" --k "* ]]; then
    for K in 1 10; do
        log "sensitivity: k = ${K}"
        "${PY[@]}" -m "${MOD}" --k "${K}" "$@"
    done
fi
log "outputs: 06_switch_mechanism/_m/longread_driver_confirm/"
log "**** complete ****"
