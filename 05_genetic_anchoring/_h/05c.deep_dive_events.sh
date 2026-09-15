#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=deep-dive-events
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%j.log
## Per-event table of the colocalized isoform events (gene_deep_dive --part events): one row per
## event, anchor -> switch -> consequence, merged with the signed risk-allele direction. Written to
## 05_genetic_anchoring/_m/deep_dive/deep_dive_events.{parquet,tsv}.
##
## Split out of the per-gene deep dive because it needs only the coloc layer, and stage 06 reads
## it (06_switch_mechanism/_h/01b-01c validate_switch_splicing, 01k junction_coloc_confirm). The
## panel, vignettes and RBP/clinical tables need stages 06 and 07 and are
## 08_integration/_h/01a.build_deep_dive.sh.
##
## Run after 04b.coloc_direction.sh:
##   sbatch 05_genetic_anchoring/_h/05c.deep_dive_events.sh
## Calls the env interpreter directly, so the job cannot die on "module: command not found".
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 05_genetic_anchoring/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

log "**** gene_deep_dive --part events ****"
"${PY}" -u -m isograph_benchmark.real_data.gene_deep_dive --part events "$@"
log "**** complete ****"
