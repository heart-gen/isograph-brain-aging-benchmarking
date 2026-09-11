#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-qtl-checks
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/brainseq-qtl-checks-%j.log
#
# Validation checks on a finished BrainSEQ switch-QTL arm (_h/25), required before its results
# are read:
#   signpin  -- orient the recomputed S_g like the discovery fit; record n_genes_sign_flipped and
#               write pinned slopes beside the permutation results (never over them).
#   control  -- positive control for the A_g eQTL arm against tissue-matched GTEx v11 brain
#               eGenes: gene-level pi1 and direction concordance at GTEx's lead variant. The
#               pass rule is fixed in isograph_benchmark/real_data/brainseq_qtl_checks.py.
#
# Memory on PSC is --cpus-per-task x 2000MB; 8 cpus = 16 GB. The control reads each region's
# nominal cis shards (~6 GB per region) by predicate pushdown on the GTEx lead rsIDs and scans
# the per-chromosome .pvar once for their alleles.
#
# Usage (from the repo root):
#   sbatch 05_genetic_anchoring/_h/28.brainseq_qtl_checks.sh
#   sbatch --export=ALL,SWQTL_ARM=ea_only 05_genetic_anchoring/_h/28.brainseq_qtl_checks.sh
# Smoke test on one chromosome (login-node safe):
#   python -m isograph_benchmark.real_data.brainseq_qtl_checks --stage control --chrom 22
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 05_genetic_anchoring/_m/logs

ARM="${SWQTL_ARM:-all_samples}"
STAGE="${CHECK_STAGE:-all}"

# Guard: some Bridges2 batch nodes start without Lmod initialised.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log "**** BrainSEQ QTL checks: arm=${ARM} stage=${STAGE} ****"
python -m isograph_benchmark.real_data.brainseq_qtl_checks --stage "${STAGE}" --arm "${ARM}"
conda deactivate
log "**** Complete: arm=${ARM} ****"
