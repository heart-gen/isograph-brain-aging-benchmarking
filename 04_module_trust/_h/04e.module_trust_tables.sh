#!/usr/bin/env bash
## Display tables and figures for the module-reproducibility subsection.
##
## A LOCAL script on purpose (no #SBATCH header): it re-reads committed ledgers that are a
## few hundred rows each, so it runs in seconds on a laptop and needs no allocation. It
## fits nothing -- every number it writes is copied or counted from a parquet already in
## the tree.
##
## Waits on: _h/02c (trust gate), _h/03b (within-cohort), _h/03d (permutation nulls),
## _h/04c (projection aggregate) and 08_integration/_h/01c (functional preservation).
## Re-run it after any of those, or the figure will carry stale counts.
##
## Writes 04_module_trust/_m/stability/module_trust_tables/{funnel_claims,
## split_half_modules,split_half_pairs,projection_modules,projection_summary,
## functional_preservation,crosscohort_permutation,region_funnel,projection_sign_scale,
## resolution_sensitivity,driver_structure}.csv and
## MODULE_TRUST_TABLES.md, then the figures manuscript/_m/figures/figTrustFunnel.{pdf,png}
## and figDriverStructure.{pdf,png}.
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

log_message "**** module trust display tables ****"
"${PY}" -u -m isograph_benchmark.real_data.module_trust_tables "$@"

log_message "**** figure: module reproducibility ****"
"${RSCRIPT}" manuscript/_h/trust_funnel_figure.R

log_message "**** figure: driver structural classes ****"
"${RSCRIPT}" manuscript/_h/driver_structure_figure.R

log_message "**** Complete ****"
