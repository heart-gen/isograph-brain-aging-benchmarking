#!/usr/bin/env bash
## Composition-vs-age coupling: is the composition the adjustment conditions on itself
## age-coupled, and does that coupling explain which regions lose their switch-unique genes?
##
## A LOCAL script on purpose (no #SBATCH header): this stage is a pure read over artifacts
## already on disk -- the committed MuSiC fractions and the bundle sample tables -- so it
## runs in seconds on a laptop and needs no allocation. Run it after the composition arrays
## and the rollup:
##   bash 03_module_characterization/_h/01g.celltype_composition_brainseq.sh   # or on SLURM
##   bash 03_module_characterization/_h/01h.celltype_composition_gtex.sh
##   bash 03_module_characterization/_h/02c.composition_meta.sh
##   bash 03_module_characterization/_h/04b.composition_age_coupling.sh
##
## Writes 03_module_characterization/_m/composition_age_{coupling.parquet,coupling.csv,
## samples.csv,persistence.csv}, COMPOSITION_AGE_COUPLING.md, and the figure
## manuscript/_m/figures/figCompositionAgeCoupling.{pdf,png}.
##
## Interpreters are taken from the environment so the same script serves the cluster and a
## local checkout: set ISOGRAPH_PY / ISOGRAPH_RSCRIPT to override the defaults below.
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

log_message "**** composition-vs-age coupling ****"
"${PY}" -u -m isograph_benchmark.real_data.composition_age_coupling "$@"

log_message "**** figure ****"
"${RSCRIPT}" manuscript/_h/composition_age_coupling_figure.R

log_message "**** Complete ****"
