#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-mod-meta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=05_genetic_anchoring/_m/logs/coloc-mod-meta-%j.log
#
# Per-gene sQTL-vs-eQTL colocalization contrast, step 3: the paired statistics.
#
# Primary binary test is an exact McNemar on the DISCORDANT genes (splicing-only vs
# expression-only); primary continuous test is a paired Wilcoxon on the conditional
# posterior PP4/(PP3+PP4), which divides out the eQTL-vs-sQTL discovery-power
# difference. Sensitivity arms vary the p12 prior, the PP4 call and the shared-SNP
# floor. Pure joins over a small table; a couple of minutes.

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
mkdir -p 05_genetic_anchoring/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

# Arm defaults to `switch`; pass --arm background (etc.) through "$@", or via
#   sbatch 05_genetic_anchoring/_h/21.coloc_modality_meta.sh --arm background
log_message "**** coloc modality contrast: meta ${*:-(switch)} ****"
python -m isograph_benchmark.real_data.coloc_modality_contrast --stage meta "$@"
conda deactivate
log_message "**** Complete ****"
