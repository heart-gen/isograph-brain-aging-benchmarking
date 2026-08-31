#!/usr/bin/env bash
# Stage the collapse-fix C (resolution-sweep giant cap) split-half A/B.
#
# Runs a SELF-CONTAINED A/B in an isolated sandbox partitions dir so the committed
# production partitions (03_module_trust/_m/stability/partitions) are untouched, and the
# baseline is regenerated on CURRENT main (weighted/seeded Leiden) so the only
# difference vs the candidate is the giant cap. The Jun-12 production baseline
# predates the determinism fix (d603938, Jun-21) and would otherwise confound the
# comparison (two changes, not one).
#
#   Arm A (baseline):  isograph         leiden_resolution=2.0, no cap
#   Arm B (candidate): isograph_gcap15  leiden_max_giant_frac=0.15 (sweep 2,4,8,16,32,64)
# across the 6 trust-funnel regions x 5 seeds x 2 halves.
#
# Acceptance (see SOFTWARE_ROBUSTNESS_PLAN.md S1.3 / S6.B): cap arm giant-fraction
# <= 0.15 on all 6 regions AND within-cohort ARI/NMI NOT worse than baseline. (The
# mechanically-similar post-hoc max_module_frac REGRESSED ARI -- this resolution
# sweep is the re-test before any production default flip.)
#
# Usage (LOGIN NODE -- this only submits jobs, it is not itself a SLURM job):
#   bash 03_module_trust/_h/08.gcap_ab.sh
set -euo pipefail

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${PWD}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: run from the repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"

SANDBOX="${PROJECT_ROOT}/03_module_trust/_m/stability_gcap_ab/partitions"
mkdir -p "${SANDBOX}" 03_module_trust/_m/stability/logs
echo "A/B sandbox partitions: ${SANDBOX}"

DRIVER=03_module_trust/_h/01.stability_isograph.sh
# 48G (24 x 2000M) so dlpfc (array task 3) does not OOM; 6h covers 10 serial fits on
# the largest region. CLI overrides the driver's 32G/3h directives; --array=1-6 stays.
COMMON=(--account=bio260021p --cpus-per-task=24 --time=06:00:00)

JID_A=$(sbatch "${COMMON[@]}" --job-name=gcapAB-baseline \
    --export=ALL,STABILITY_PARTITIONS_DIR=${SANDBOX} \
    "${DRIVER}" | awk '{print $NF}')
echo "baseline arm  (isograph):        ${JID_A}"

JID_B=$(sbatch "${COMMON[@]}" --job-name=gcapAB-cap15 \
    --export=ALL,STABILITY_PARTITIONS_DIR=${SANDBOX},STABILITY_LEIDEN_GIANT_FRAC=0.15 \
    "${DRIVER}" | awk '{print $NF}')
echo "candidate arm (isograph_gcap15): ${JID_B}"

JID_AGG=$(sbatch --account=bio260021p \
    --dependency=afterok:${JID_A}:${JID_B} \
    --export=ALL,STABILITY_PARTITIONS_DIR=${SANDBOX} \
    03_module_trust/_h/03.stability_aggregate.sh | awk '{print $NF}')
echo "aggregate (afterok ${JID_A}:${JID_B}): ${JID_AGG}"
echo "summary.json will land under: ${SANDBOX%/partitions}"
