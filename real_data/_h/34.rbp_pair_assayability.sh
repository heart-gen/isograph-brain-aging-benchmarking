#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=rbp-pair-assayability
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=real_data/_m/logs/%x-%j.log

## Score the frozen IsoGraph switch pairs of each perturbation panel for bench assayability:
## whether an isoform-ratio assay can resolve the two transcripts (unique splice junction or
## unique exonic segment, from the GENCODE v47 structures), whether both isoforms are actually
## expressed in GTEx brain at a ratio with headroom, and whether that ratio moves with age
## (logit carrier fraction ~ age + SEX + RIN + ischemic time, BH within RBP).
##
## Depends on real_data/_h/build_rbp_target_panel.sh having been run for each --rbp, so
## real_data/_m/rbp_target_panel/<RBP>/rbp_target_switch_pairs.parquet exists. Reads the 13
## GTEx RSEM TPM tables (~85 MB each), which is why this is a batch job and not a login-node
## script.
##
## Usage: sbatch real_data/_h/34.rbp_pair_assayability.sh
##        sbatch real_data/_h/34.rbp_pair_assayability.sh --rbp QKI --min-iqr 0.08

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

# The three panel arms of WETLAB_PERTURBATION_DESIGN.md; override wholesale by passing --rbp.
ARGS=("$@")
if [[ ${#ARGS[@]} -eq 0 ]]; then
    ARGS=(--rbp NONO --rbp ELAVL1 --rbp KHDRBS1)
fi

log_message "**** RBP switch-pair assayability starts ****"
python -u -m isograph_benchmark.real_data.rbp_pair_assayability "${ARGS[@]}"
conda deactivate
log_message "**** Complete -> real_data/_m/rbp_target_panel/<RBP>/rbp_pair_assayability.* ****"
