#!/usr/bin/env bash
#SBATCH --job-name=ase-gwas-lead
#SBATCH --partition=RM-shared
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=06_switch_mechanism/_m/logs/ase-gwas-lead-%j.log

## BRIDGES-2 (PI item 10a, second arm). Step 5 fits at the switch-QTL lead and transfers the
## sign to the disease allele through signed LD, which fails for most rows (336 of 556 below
## the |r| >= 0.8 gate or with no LD), so its negative reads as "we could not ask". This arm
## anchors the same within-donor contrast on the GWAS lead ITSELF: the risk haplotype comes
## straight from the donor's phased genotype, orientation is exact by construction, and the
## cost is power rather than validity.
##
## No recount: it relabels the haplotype counts already in allelic_donor_counts.parquet. The
## anchor genotypes are read from phASER's own per-sample VCFs (shared/resources/
## processed-data/ase-files/<region>/phASER/), which is the frame the counts are phased in --
## the TOPMed QTL panel is phased too, but in a DIFFERENT frame (agreement with phASER 0.51),
## so it must not be used for this. The run self-checks that before fitting anything.
##
##   sbatch 06_switch_mechanism/_h/05b.ase_gwas_lead_arm.sh [--region caudate]
##
## Writes gwas_lead_arm.parquet + GWAS_LEAD_ARM.md beside the allelic outputs. Reads ~1,400
## small indexed VCFs, so it is IO-bound rather than CPU-bound.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 06_switch_mechanism/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge 2>/dev/null || true
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

log "**** Job starts ****"
python -m isograph_benchmark.real_data.ase_gwas_lead_arm "$@"
log "**** Complete -> 06_switch_mechanism/_m/ase_junction_switch/<region>/ ****"
