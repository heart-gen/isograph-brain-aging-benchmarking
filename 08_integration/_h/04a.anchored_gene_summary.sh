#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=anchored-summary
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=08_integration/_m/logs/%x-%j.log

## One table for the 30 genes whose colocalizing junction resolves to an IsoGraph switch
## pair: all of the genetics (coloc.abf as PRIMARY, SMR/HEIDI, eCAVIAR CLPP as a trailing
## secondary column) beside all of the orthogonal validation (ONT long-read, the BrainSEQ
## allele-aware junction recount, the LIBD PSI arm), sorted on coloc.
##
## This lives in 08_integration because it reads stage 05 (coloc, coloc_signal_susie, SMR)
## and stage 06 (longread_switch_confirm, ase_junction_switch, junction_coloc_confirm)
## together. Re-run it whenever any of those is re-run: it is the only place their numbers
## are compared, so a stale copy here is worse than none.
##
## Every number is generated. `provenance.json` records the git commit, whether the tree was
## dirty, package versions, and each input's size, mtime and digest, so a reader can tell
## whether the table still matches its inputs. Nothing in ANCHORED_GENE_SUMMARY.md is to be
## edited by hand.
##
## ~8 GB peak: the junction recount count tables are ~16.5M rows per BrainSEQ region and are
## filtered to the relevant switch pairs on read.
## Usage: sbatch 08_integration/_h/04a.anchored_gene_summary.sh
##    or: bash  08_integration/_h/04a.anchored_gene_summary.sh   (interactive; it is light)
##
## Requires: 05_genetic_anchoring coloc + smr_heidi, and 06_switch_mechanism
## longread_switch_confirm + ase_junction_switch + junction_coloc_confirm.

set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 08_integration/_m/logs

# Non-interactive shells do not have `module` on the path; older wrappers silently lost
# array tasks to "module: command not found".
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log "**** anchored gene summary: genetics + orthogonal validation ****"
python -u -m isograph_benchmark.real_data.anchored_gene_summary "$@"

conda deactivate 2>/dev/null || true
log "done -> 08_integration/_m/anchored_gene_summary"
