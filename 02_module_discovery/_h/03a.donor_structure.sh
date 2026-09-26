#!/usr/bin/env bash
## Donor sharing across the 17 analyses: the evidence behind "17 applications, not 17
## independent cohorts".
##
## A LOCAL script on purpose (no #SBATCH header): this reads only the bundle sample tables
## already in inputs/bundles, so it runs in seconds on a laptop and needs no allocation. It
## depends on nothing downstream -- the bundles are enough -- so it can be run at any point
## after 00a.
##
## Writes 02_module_discovery/_m/donor_structure/{donor_overlap.parquet,donor_overlap.csv,
## analysis_donor_counts.csv,donor_analysis_counts.csv,donor_incidence.csv,
## donor_structure_summary.json,DONOR_STRUCTURE.md} and the figure
## manuscript/_m/figures/figDonorStructure.{pdf,png}.
##
## Interpreters come from the environment so the same script serves the cluster and a local
## checkout: set ISOGRAPH_PY / ISOGRAPH_RSCRIPT to override the defaults below.
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

log_message "**** donor structure ****"
"${PY}" -u -m isograph_benchmark.real_data.donor_structure "$@"

log_message "**** figure ****"
"${RSCRIPT}" manuscript/_h/donor_structure_figure.R

log_message "**** Complete ****"
