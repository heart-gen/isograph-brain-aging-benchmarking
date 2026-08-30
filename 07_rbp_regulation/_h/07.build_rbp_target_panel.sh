#!/usr/bin/env bash
# Reproducibly build the RBP perturbation target panel (default KHDRBS1) from the Stage-2
# RBP-regulon outputs. Joins rbp_switch_calls (gene-level motif gain/loss) to rbp_regulon
# (module-level enrichment), ranks the eligible genes, tiers them, and writes the
# collaborator-facing expression-check list.
#
# Depends on: rbp_scan.py -> rbp_regulon.py -> rbp_binding.py (all already run), plus the
# per-region isograph_vae artifacts and 05_genetic_anchoring/_m/deep_dive/deep_dive_panel.parquet.
# Light enough for a login node (a few minutes; the motif count tables are large).
#
# Usage (from repo root):
#   bash 07_rbp_regulation/_h/07.build_rbp_target_panel.sh
#   bash 07_rbp_regulation/_h/07.build_rbp_target_panel.sh --rbp QKI --n-core 10
set -euo pipefail
log() { echo "$(date '+%H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${PWD}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: run from the repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"

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
log "done -> 07_rbp_regulation/_m/rbp_target_panel"
