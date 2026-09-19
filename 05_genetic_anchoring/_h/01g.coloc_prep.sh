#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc-prep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
## RM-shared allocates ~2 GB/cpu, so memory is scaled via cpus, not --mem: 12 cpus ~ 24 GB.
#SBATCH --output=05_genetic_anchoring/_m/logs/coloc-prep-%j.log

## Coloc step 1 (wrapper) — select switch genes under a trait's GWAS peaks + their GTEx
## brain QTL credible sets, and write per-locus GWAS z for SuSiE. Reproducible entry
## point for coloc_prep.py (was previously run interactively). Build-agnostic: hg38
## neurodegeneration sumstats are joined to the hg19 LD panel by rsID.
##
##   disease : sbatch 05_genetic_anchoring/_h/01g.coloc_prep.sh --gene-source brainseq-sczd --trait scz
##   aging   : sbatch 05_genetic_anchoring/_h/01g.coloc_prep.sh --gene-source aging --trait ad --min-recurrence 1
##
## Output dir <gene_source>__<trait>; then 05_genetic_anchoring/_h/02e.locus_ld.sh <dir> and 05_genetic_anchoring/_h/03b.coloc_clpp.sh <dir>.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 05_genetic_anchoring/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log "**** coloc_prep $* ****"
python -m isograph_benchmark.real_data.coloc_prep "$@"
conda deactivate
log "**** coloc prep done ****"
