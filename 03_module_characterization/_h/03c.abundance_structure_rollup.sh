#!/usr/bin/env bash
## Pool the per-store abundance-vs-switch orthogonality tables (from the _h/03b array) into
## 03_module_characterization/_m/axis_orthogonality_{all.parquet,summary.parquet,summary.csv}
## + AXIS_ORTHOGONALITY.md, then rebuild figSeparation (S-real-4), whose panel A is now
## faceted over all 17 analyses.
##
## A LOCAL script on purpose (no #SBATCH header): the rollup is a concat of 17 small tables
## and runs in seconds. Run it after _h/03a (example gene + incremental summary, which the
## figure's other panels read) and after every task of _h/03b has finished.
##
## Interpreters come from the environment so the same script serves the cluster and a local
## checkout: set ISOGRAPH_PY / ISOGRAPH_RSCRIPT to override the defaults below. Extra
## arguments are forwarded to the rollup, e.g. `--variant with-abundance`.
set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${PWD}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: run from the repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

PY="${ISOGRAPH_PY:-python}"
RSCRIPT="${ISOGRAPH_RSCRIPT:-Rscript}"
command -v "${PY}" >/dev/null || { echo "ERROR: python not found (set ISOGRAPH_PY)."; exit 1; }
command -v "${RSCRIPT}" >/dev/null || { echo "ERROR: Rscript not found (set ISOGRAPH_RSCRIPT)."; exit 1; }

log_message "**** axis orthogonality rollup ****"
"${PY}" -u -m isograph_benchmark.real_data.abundance_structure_separation --rollup "$@"

log_message "**** figure ****"
"${RSCRIPT}" manuscript/_h/abundance_structure_figure.R

log_message "**** Complete ****"
