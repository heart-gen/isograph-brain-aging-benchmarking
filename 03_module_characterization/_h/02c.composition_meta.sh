#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=composition-meta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=03_module_characterization/_m/logs/%x-%j.log

## Cell-type composition rollup (celltype_composition `meta`): the with-vs-without
## composition adjustment contrast across BrainSEQ and GTEx.
##
## Previously a comment at the bottom of _h/01g and _h/01h and run by hand. Run after BOTH
## composition arrays finish:
##   b=$(sbatch --parsable 03_module_characterization/_h/01g.celltype_composition_brainseq.sh)
##   g=$(sbatch --parsable 03_module_characterization/_h/01h.celltype_composition_gtex.sh)
##   sbatch --dependency=afterok:${b}:${g} 03_module_characterization/_h/02c.composition_meta.sh
##
## Writes 03_module_characterization/_m/{composition_adjustment.parquet,COMPOSITION_ADJUSTMENT_SUMMARY.md}
## and the GTEx rollup composition_adjustment_gtex.parquet in the GTEx MuSiC directory.
## Calls the env interpreter directly rather than `module load` + `conda activate`.
## Extra arguments are forwarded, e.g. `--variant standard`.
set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 03_module_characterization/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log_message "**** composition adjustment rollup ****"
"${PY}" -u -m isograph_benchmark.real_data.celltype_composition meta "$@"
log_message "**** Complete ****"
