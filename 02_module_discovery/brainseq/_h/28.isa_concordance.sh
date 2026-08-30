#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=isa-concordance
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=08:00:00
#SBATCH --output=02_module_discovery/brainseq/_m/logs/isa-concordance-%j.log

# Orthogonal switch-detection concordance: IsoGraph co-switch modules vs satuRn
# DTU (the switch test IsoformSwitchAnalyzeR v2 uses internally). Runs the three
# BrainSEQ aging regions + the caudate SCZD dx contrast by default. Set
# ISA_RUN_GTEX=1 to also run the 13 GTEx aging regions.
#
# Submit from the repository root:
#   sbatch 02_module_discovery/brainseq/_h/28.isa_concordance.sh
#   ISA_RUN_GTEX=1 sbatch 02_module_discovery/brainseq/_h/28.isa_concordance.sh

set -euo pipefail
umask 0027

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repository root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi

PYTHON="/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python"
RSCRIPT="/ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript"
if [[ ! -x "${PYTHON}" || ! -x "${RSCRIPT}" ]]; then
    echo "ERROR: missing the isograph Python or rnaseq Rscript interpreter."
    exit 1
fi

export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
CORES="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p 02_module_discovery/brainseq/_m/logs

log_message "**** ISA(satuRn) switch-detection concordance ****"

# BrainSEQ aging (3 regions) + caudate SCZD dx.
"${PYTHON}" -m isograph_benchmark.real_data.isa_concordance \
    --cohort brainseq --trait age --cores "${CORES}" "$@"
"${PYTHON}" -m isograph_benchmark.real_data.isa_concordance \
    --cohort brainseq --trait dx --cores "${CORES}" "$@"

if [[ "${ISA_RUN_GTEX:-0}" == "1" ]]; then
    log_message "**** GTEx aging (13 regions) ****"
    "${PYTHON}" -m isograph_benchmark.real_data.isa_concordance \
        --cohort gtex --trait age --cores "${CORES}" "$@"
fi

log_message "outputs: real_data/_m/isa_concordance/<cohort>_<region>_<trait>/"
log_message "**** complete ****"
