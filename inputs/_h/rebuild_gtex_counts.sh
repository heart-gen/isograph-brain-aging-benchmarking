#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=gtex-rebuild-counts
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --output=real_data/gtex/_m/logs/rebuild-counts-%j.log

# Rebuild the GTEx v11 brain bundles from RSEM transcript expected_count (count
# scale) instead of TPM, with a CPM>=1 gene filter for BrainSEQ parity.
#   Step 1: convert the 7GB expected_count GCT -> per-region transcript_reads.parquet
#   Step 2: rebuild all 13 bundles (transcript_counts + gene_counts)
# Fixes the VAE divergence (TPM) and OOM (loose TPM>=0.1 filter, ~2x features).

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/gtex/_m/logs

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "**** Step 1: convert transcript expected_count GCT -> transcript_reads.parquet ****"
python -c "from isograph_benchmark.inputs.build_parquet import convert_gtex_transcript_reads; convert_gtex_transcript_reads()"

log_message "**** Step 2: rebuild 13 GTEx bundles from counts ****"
python -c "
from isograph_benchmark.inputs.build_bundles import build_gtex_bundle
from isograph_benchmark.paths import rel
root = rel('inputs', 'processed', 'gtex_v11')
regions = sorted(p.name for p in root.iterdir() if p.is_dir())
for r in regions:
    print(f'  building {r} ...', flush=True)
    build_gtex_bundle(r)
print(f'rebuilt {len(regions)} GTEx bundles')
"

conda deactivate
log_message "**** Complete ****"
