#!/usr/bin/env bash
# Launch validation gate B.2 — the single-LR/optimizer config test.
#
# Submits one FULL-data IsoGraph fit per region (6 trust-funnel + the diverging GTEx
# nucleus_accumbens) at a single fixed LR with gradient clipping, then a chained aggregate
# that emits the PASS/FAIL gate. Orthogonal to the giant-cap A/B (05.gcap_ab.sh): that
# tests module *detection*, this tests VAE *optimizer* stability — they can run in parallel.
#
# Gate (SOFTWARE_ROBUSTNESS_PLAN.md S6.B): one documented LR trains all 7 regions to
# BrainSEQ-range RMSE (band 1.00-1.12) with no divergence and no OOM. PASS means the merged
# grad_clip_norm + divergence guard remove the need for per-dataset LR tuning; FAIL means
# S2/S3/S4 (auto-LR-backoff / free-bits / density cap) are still needed.
#
# Usage (LOGIN NODE — this only submits jobs, it is not itself a SLURM job):
#   bash real_data/stability/_h/08.lr_validation_launch.sh
#   STABILITY_LR=5e-4 bash real_data/stability/_h/08.lr_validation_launch.sh   # try another LR
set -euo pipefail

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${PWD}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: run from the repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
mkdir -p real_data/stability/_m/logs real_data/stability/_m/lr_validation

# Forward an optional single-LR / grad-clip override to the array driver.
EXPORTS="ALL"
[[ -n "${STABILITY_LR:-}" ]]        && EXPORTS="${EXPORTS},STABILITY_LR=${STABILITY_LR}"
[[ -n "${STABILITY_GRAD_CLIP:-}" ]] && EXPORTS="${EXPORTS},STABILITY_GRAD_CLIP=${STABILITY_GRAD_CLIP}"

JID=$(sbatch --export="${EXPORTS}" \
    real_data/stability/_h/06.lr_validation.sh | awk '{print $NF}')
echo "single-LR fits (array 1-7): ${JID}"

JID_AGG=$(sbatch --dependency=afterok:${JID} \
    real_data/stability/_h/07.lr_aggregate.sh | awk '{print $NF}')
echo "aggregate (afterok ${JID}): ${JID_AGG}"
echo "lr_validation_summary.{parquet,json} -> real_data/stability/_m/lr_validation/"
