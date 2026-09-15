#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=rbp-target-panel
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=08_integration/_m/logs/%x-%j.log
# Reproducibly build the RBP perturbation target panel (default KHDRBS1) from the Stage-2
# RBP-regulon outputs. Joins rbp_switch_calls (gene-level motif gain/loss) to rbp_regulon
# (module-level enrichment), ranks the eligible genes, tiers them, and writes the
# collaborator-facing expression-check list.
#
# Depends on: stage 07 (rbp_scan -> rbp_regulon -> rbp_binding), the per-region isograph_vae
# artifacts, and 08_integration/_m/deep_dive/deep_dive_panel.parquet (01a.build_deep_dive.sh).
# A few minutes; the motif count tables are large (8 cpus = 16 GB; do NOT pass --mem).
#
# Usage (from repo root):
#   sbatch 08_integration/_h/02a.build_rbp_target_panel.sh
#   sbatch 08_integration/_h/02a.build_rbp_target_panel.sh --rbp QKI --n-core 10
set -euo pipefail
log() { echo "$(date '+%H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: run from the repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 08_integration/_m/logs

# ---- env init (non-interactive shells need module + conda hook sourced explicitly) ----
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge 2>/dev/null || true
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log "rbp_target_panel.py ${*:-(defaults: --rbp KHDRBS1)}"
python -m isograph_benchmark.real_data.rbp_target_panel "$@"

conda deactivate 2>/dev/null || true
log "done -> 08_integration/_m/rbp_target_panel"
