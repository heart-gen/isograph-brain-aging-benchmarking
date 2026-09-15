#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=smr-meta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=05_genetic_anchoring/_m/logs/smr-meta-%j.log
#
# SMR + HEIDI, assembly step: turn the per-task .smr/.msmr tables that 10a.smr_heidi.sh wrote
# into smr_results.parquet, smr_runs.tsv, snp_attrition.parquet and SMR_HEIDI.md.
#
# Split out of 10a so that every meta run is a committed, re-runnable command rather than
# something typed on the login node. The stage is a pure function of what is already on disk:
# it reads the .smr/.msmr tables, targets.parquet and the coloc hierarchy, and re-running it
# overwrites its own outputs in place. No RNG, no compute -- joins over a small table, a
# couple of minutes -- so it is also safe to run directly with `bash` on a login node.
#
# The arm is named by flags, exactly as in the --stage smr run that produced the inputs; the
# flags pick the output directory, so a mismatch reads an empty tree and exits rather than
# silently summarizing the wrong arm. Pass them through "$@":
#
#   Primary (5e-8, single-SNP):
#     sbatch 05_genetic_anchoring/_h/11a.smr_heidi_meta.sh
#     sbatch 05_genetic_anchoring/_h/11a.smr_heidi_meta.sh --qtl-source brainseq --arm ea_only
#
#   Sensitivity, relaxed instrument threshold (-> sensitivity/peqtl_smr_1e-06/):
#     sbatch 05_genetic_anchoring/_h/11a.smr_heidi_meta.sh --peqtl-smr 1e-6
#     sbatch 05_genetic_anchoring/_h/11a.smr_heidi_meta.sh --peqtl-smr 1e-6 --qtl-source brainseq --arm ea_only
#
#   Sensitivity, multi-SNP SMR (-> sensitivity/smr_multi/; reads .msmr, adds p_SMR_multi):
#     sbatch 05_genetic_anchoring/_h/11a.smr_heidi_meta.sh --smr-multi
#     sbatch 05_genetic_anchoring/_h/11a.smr_heidi_meta.sh --smr-multi --qtl-source brainseq --arm ea_only
#
# Depends on: 10a.smr_heidi.sh for the same arm, and -- because the agreement column joins SMR
# onto the colocalization it is corroborating -- coloc_signal_susie --stage meta --sqtl all
# (GTEx) or coloc_brainseq --stage meta (BrainSEQ).

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
mkdir -p 05_genetic_anchoring/_m/logs

# Guard: some Bridges2 batch nodes start array tasks without Lmod initialised.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** SMR meta: ${*:-(gtex, primary)} ****"
python -m isograph_benchmark.real_data.smr_heidi --stage meta "$@"
conda deactivate
log_message "**** Complete ****"
