#!/usr/bin/env bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=ase-junc-count
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=04:00:00
#SBATCH --array=1-498%60
#SBATCH --output=06_switch_mechanism/_m/logs/ase-junc-count-%A_%a.log

## QUEST ONLY (PI item 10a). One WASP BAM per task: fragment-level isoform x haplotype
## counts per switch pair, read in place from star-wasp/aligned_bams (never copied). The
## task id is the 1-based row of _m/ase_junction_switch/<region>/samples.tsv written by 01l.
## The header's --array fits DLPFC (498); ase_junction_quest.sh sizes it per region from
## the manifest, and tasks past the end of samples.tsv exit cleanly.
##
##   sbatch --array=1-486%60 06_switch_mechanism/_h/02f.ase_junction_count.sh --region caudate
##   sbatch --array=17,203 06_switch_mechanism/_h/02f.ase_junction_count.sh --region caudate  # reruns
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log "**** Job starts: task ${SLURM_ARRAY_TASK_ID} on ${HOSTNAME} ****"
python -m isograph_benchmark.real_data.ase_junction_switch --stage count \
    --task "${SLURM_ARRAY_TASK_ID}" "$@"
log "**** Job ends ****"
