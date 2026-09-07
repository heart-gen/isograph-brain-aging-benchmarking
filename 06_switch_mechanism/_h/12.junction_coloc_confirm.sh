#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=junction-coloc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/junction-coloc-%j.log

## Short-read BrainSEQ junction confirmation of the genetically anchored switch pairs
## (SNCA, CTSH) -- the test the ONT long-read data could not run: the junction the sQTL
## actually tags, in the region the coloc was found in, at ~40x the long-read n.
##
## Memory on Bridges is --cpus-per-task x 2000MB, so 16 cpus = 32G. The three PSI event
## tables are 540-620 MB of parquet each and are loaded whole (one region at a time), so
## this needs the headroom; do NOT pass --mem.
##
## Extra arguments are forwarded, e.g.
##   sbatch 06_switch_mechanism/_h/12.junction_coloc_confirm.sh --regions hippocampus
##   sbatch 06_switch_mechanism/_h/12.junction_coloc_confirm.sh --min-usage 0.02
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 06_switch_mechanism/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

## The PSI tables are git-LFS tracked; a fresh clone leaves pointers behind.
for REG in caudate hippocampus dlpfc; do
    F="inputs/processed/brainseq/${REG}/psi_events.parquet"
    if [[ -s "${F}" ]] && head -c 40 "${F}" | grep -q 'git-lfs'; then
        log "restoring ${F} from git-lfs"
        git lfs pull --include="${F}"
    fi
done

log "junction_coloc_confirm $*"
"${PY}" -m isograph_benchmark.real_data.junction_coloc_confirm "$@"
log "**** complete ****"
