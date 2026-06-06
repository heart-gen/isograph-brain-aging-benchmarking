#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-leiden-sweep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=real_data/brainseq/_m/logs/leiden-sweep-%j.log

# Re-clusters saved IsoGraph edges at a grid of Leiden resolutions WITHOUT
# refitting the VAE (reuses edges.parquet + feature_scores.parquet). Reports
# n_modules / giant fraction / trait-association counts (and DRD2 for SCZD
# caudate) at each resolution so we can pick a biologically interpretable
# module count (WGCNA reference: ~23 modules).
#
# Requires the base runs (01.run_isograph_aging.sh / 02.run_isograph_sczd.sh)
# to have already written artifacts to real_data/brainseq/<region>/_m/.
#
# Flags forwarded to the sweep (pass after the script name):
#   --dry-run      compute + print only, write nothing
#   --write-best   commit the best-resolution artifacts (selected by GO density)
#   --no-go        skip GO enrichment (use trait-association fallback for --write-best)
#   --resolutions  override the default per-analysis resolution grid

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

mkdir -p real_data/brainseq/_m/logs
log_message "**** BrainSeq Leiden resolution sweep ****"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "Sweeping SCZD caudate"
python -m isograph_benchmark.real_data.sweep_leiden brainseq-sczd "$@"

log_message "Sweeping aging regions (caudate, hippocampus, dlpfc)"
python -m isograph_benchmark.real_data.sweep_leiden brainseq-aging "$@"

conda deactivate
log_message "**** Complete ****"
