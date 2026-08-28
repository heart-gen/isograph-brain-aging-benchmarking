#!/usr/bin/env bash
# ENCODE eCLIP peak acquisition for the RBP binding-evidence stage.
#
# NOT a SLURM job: PSC compute nodes have no outbound network, so this must run on the login
# node. It is committed anyway because the analysis it feeds is paper-facing and the repo
# rule is that every stage is reproducible from a command rather than from a note in a
# comment (20.rbp_binding.sh previously documented this step in prose only).
#
# Downloads narrowPeak BEDs for every RBP nominated by the motif regulon analysis into
# inputs/raw/rbp_binding/<RBP>.bed.gz, skipping any already present. Re-run it whenever the
# regulon nominations change; then run real_data/brainseq/_h/20.rbp_binding.sh.
#
# Usage (from the repo root, on a login node):
#   bash real_data/brainseq/_h/32.rbp_binding_fetch.sh

set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${PWD}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: run from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo "ERROR: compute nodes have no outbound network; run this on a login node."
    exit 1
fi

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "**** ENCODE eCLIP fetch ****"
"${PY}" -u -m isograph_benchmark.real_data.rbp_binding fetch "$@"
log "**** Complete — now run real_data/brainseq/_h/20.rbp_binding.sh ****"
