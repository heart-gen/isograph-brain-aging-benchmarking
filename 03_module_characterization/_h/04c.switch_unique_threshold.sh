#!/usr/bin/env bash
## Threshold sensitivity of the switch-unique classification: does the regional pattern
## survive moving the FDR line off 0.10?
##
## A LOCAL script on purpose (no #SBATCH header): this is a recount over the gene-level
## tables `incremental_association` already wrote. The BH-adjusted values do not depend on
## alpha, so nothing is refit and it runs in seconds on a laptop. Run it after the
## incremental-association arrays (`_h/01e`, `_h/01f`) and, for the adjusted arm, after the
## composition arrays (`_h/01g`, `_h/01h`).
##
## Writes 03_module_characterization/_m/switch_unique_threshold{.parquet,.csv,_rank.csv,
## _calls.csv,_summary.json}, SWITCH_UNIQUE_THRESHOLD.md, and the figure
## manuscript/_m/figures/figSwitchUniqueThreshold.{pdf,png}.
##
## Interpreters come from the environment so the same script serves the cluster and a local
## checkout: set ISOGRAPH_PY / ISOGRAPH_RSCRIPT to override the defaults below. Extra
## arguments are forwarded, e.g. `--alpha 0.01`.
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

log_message "**** switch-unique threshold sensitivity ****"
"${PY}" -u -m isograph_benchmark.real_data.switch_unique_threshold "$@"

log_message "**** figure ****"
"${RSCRIPT}" manuscript/_h/switch_unique_threshold_figure.R

log_message "**** Complete ****"
