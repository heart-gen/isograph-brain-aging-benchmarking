#!/usr/bin/env bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=ase-junc-allelic
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=04:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/ase-junc-allelic-%j.log

## QUEST ONLY (PI item 10a, step 4). Runs only where 03b's gate passed (exits cleanly as the
## pre-specified negative otherwise): orients each lead-heterozygous donor to the lead's ALT
## haplotype, fits the beta-binomial GLMM (donor random intercept, LRT on beta) per pair, with
## the homozygous-at-lead null, a model-free score test and the between-donor check.
## allelic_test.parquet + ASE_JUNCTION_ALLELIC.md go back to Bridges-2 via git.
##
##   sbatch 06_switch_mechanism/_h/04a.ase_junction_allelic.sh [--region dlpfc]
##          [--risk-alleles <rsid,risk_allele TSV>]
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
# one BLAS thread per worker process
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

log "**** Job starts ****"
python -m isograph_benchmark.real_data.ase_junction_allelic --stage test \
    --threads "${SLURM_CPUS_PER_TASK:-1}" "$@"
log "**** Job ends ****"
