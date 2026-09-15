#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=xcohort-pooled
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=04_module_trust/_m/logs/xcohort-pooled-%j.log

## Q3 POOLED cross-cohort aging replication, both methods.
##
## The per-pair arm (09/10) is structurally underpowered: ~3 reproducible matched modules
## per region pair puts the binomial sign floor at 0.5^3 = 0.125, so perfect concordance
## cannot reach p < 0.05. Pooling every region pair into one test breaks that floor via a
## Stouffer directional meta-Z plus two permutation arms.
##
## This existed as a CLI with no wrapper, so it was only ever run by hand -- which is how
## both output parquets came to sit at 0 rows without anyone noticing. It is in the
## numbered sequence now so it runs, and fails visibly, with everything else.
##
## Runs after 09/10 (needs the Q1 trusted sets, the production module tables and the
## age_linear effects for both cohorts and both methods). No model fits.
##
## Bridges memory is --cpus-per-task x 2000MB, so 4 cpus = 8G; do NOT pass --mem.
## Calls the env interpreter directly rather than `module load` + `conda activate`, so it
## cannot die on "module: command not found" in a non-login shell.
##
## Extra arguments are forwarded to BOTH methods, e.g.
##   sbatch 04_module_trust/_h/03e.replication_pooled.sh --min-jaccard 0.25
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 04_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

for METHOD in isograph wgcna; do
    log "replication-pooled --method ${METHOD} $*"
    "${PY}" -m isograph_benchmark.real_data.module_trust replication-pooled \
        --method "${METHOD}" "$@"
done
log "**** complete ****"
