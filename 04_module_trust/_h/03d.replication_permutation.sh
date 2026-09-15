#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=rep-perm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --array=0-11
#SBATCH --output=04_module_trust/_m/logs/rep-perm-%A_%a.log

## Empirical null for the cross-cohort aging-replication count (reviewer item 2).
##
## Array = {isograph, wgcna} x {pearson, partial_linear, spline_f} x {age, matching}:
##   statistic  pearson        -- reproduces the PUBLISHED covariate-free Pearson statistic
##              partial_linear -- covariate-adjusted signed age term (primary)
##              spline_f       -- df=3 spline block F-test (sensitivity)
##   null       age            -- Freedman-Lane residual permutation under the age-null,
##                                preserving GTEx covariate structure; matching held fixed
##              matching       -- age statistics held fixed, BrainSEQ->GTEx assignment
##                                reshuffled within each pair's GTEx module pool
##
## B = 10,000, seed 13.  Module construction, trusted sets and the gene-Jaccard matching are
## never permuted.  Each task asserts the reconstructed eigengenes reproduce the published
## age_linear.parquet before permuting.
##
## Production runs one array per covariate mode, after the Q1 trust gate (02c) and module
## interpretation (03_module_characterization/_h/01a-01b):
##   for cov in full complement none; do
##       sbatch 04_module_trust/_h/03d.replication_permutation.sh --covariates ${cov}
##   done
## (add --min-jaccard 0.05 for the matching-threshold sensitivity). The markdown report is its own
## step, 04_module_trust/_h/04b.replication_permutation_report.sh, after all three arrays.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 04_module_trust/_m/logs

METHODS=(isograph wgcna)
STATS=(pearson partial_linear spline_f)
NULLS=(age matching)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
METHOD="${METHODS[$((IDX / 6))]}"
STAT="${STATS[$(((IDX % 6) / 2))]}"
NULL="${NULLS[$((IDX % 2))]}"

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log "**** replication permutation: method=${METHOD} statistic=${STAT} null=${NULL} ****"
python -m isograph_benchmark.real_data.replication_permutation \
  --method "${METHOD}" --statistic "${STAT}" --null "${NULL}" \
  --n-perm 10000 --seed 13 "$@"

conda deactivate
log "**** done (${METHOD}/${STAT}/${NULL}) ****"
