#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=rbp-binding-fetch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --output=07_rbp_regulation/_m/logs/rbp-binding-fetch-%j.log

## ENCODE eCLIP peak acquisition for the RBP binding-evidence stage.
##
## Downloads narrowPeak BEDs for every RBP nominated by the motif regulon analysis into
## inputs/raw/rbp_binding/<RBP>.bed.gz, skipping any already present. Re-run it whenever the
## regulon nominations change; 07_rbp_regulation/_h/04a.rbp_binding.sh consumes the result.
##
## Compute nodes reach the network, so this is an ordinary batch step in the stage DAG and
## needs no manual login-node hop. (It refused to run under SLURM until 2026-09-15 on the
## belief that compute nodes had no outbound route; verified false against ENCODE, NCBI and
## Google Storage from RM-shared.)
##
## Usage: sbatch 07_rbp_regulation/_h/03b.rbp_binding_fetch.sh
##    or: bash 07_rbp_regulation/_h/03b.rbp_binding_fetch.sh   (login node, still fine)

set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: run from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "**** ENCODE eCLIP fetch ****"
"${PY}" -u -m isograph_benchmark.real_data.rbp_binding fetch "$@"
log "**** Complete — now run 07_rbp_regulation/_h/04a.rbp_binding.sh ****"
