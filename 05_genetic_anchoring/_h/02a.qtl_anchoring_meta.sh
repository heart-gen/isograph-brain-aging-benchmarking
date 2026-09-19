#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=qtl-anchoring-meta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%j.log
## Random-effects pooling of the per-analysis xQTL anchoring (01, 02, 17) and the paired
## splicing-specificity contrast -> _m/qtl_anchoring_meta/. The arm is named by flags, exactly
## as in the per-analysis runs; each sensitivity arm lands under
## sensitivity/<outcome>_<covariate_set>/ and never over the primary:
##
##   primary      sbatch --dependency=afterok:<01>:<02> 05_genetic_anchoring/_h/02a.qtl_anchoring_meta.sh
##   constraint   sbatch ... 02a.qtl_anchoring_meta.sh --outcome binary --covariate-set constraint
##   continuous   sbatch ... 02a.qtl_anchoring_meta.sh --outcome continuous
##   dose         sbatch ... 02a.qtl_anchoring_meta.sh --outcome dose
##
## The matched WGCNA baselines have no dose arm; the meta skips a method whose per-analysis
## file is absent. Pure joins over small tables. Calls the env interpreter directly, so the
## job cannot die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 05_genetic_anchoring/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "**** qtl_anchoring_meta ${*:-(primary)} ****"
"${PY}" -u -m isograph_benchmark.real_data.qtl_anchoring_meta "$@"
log "**** complete ****"
