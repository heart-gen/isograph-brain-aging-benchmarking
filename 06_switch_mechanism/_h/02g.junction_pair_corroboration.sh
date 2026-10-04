#!/usr/bin/env bash
# Local, table-only Figure 6c analysis. No cluster paths or scheduler required.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
cd "${PROJECT_ROOT}"
if [[ -z "${PYTHON_BIN:-}" ]]; then
    if [[ -x "${HOME}/.venvs/isograph/bin/python" ]]; then
        PYTHON_BIN="${HOME}/.venvs/isograph/bin/python"
    else
        PYTHON_BIN=python3
    fi
fi
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
"${PYTHON_BIN}" -B -u -m isograph_benchmark.real_data.junction_pair_corroboration --root "${PROJECT_ROOT}" "$@"
# The plotter uses defaults; custom --output-dir analyses should invoke it explicitly.
if [[ " $* " != *" --output-dir"* ]]; then
    "${RSCRIPT_BIN}" manuscript/_h/junction_pair_corroboration_figure.R "${PROJECT_ROOT}"
fi
