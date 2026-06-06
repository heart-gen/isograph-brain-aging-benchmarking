#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=GPU-shared
#SBATCH --job-name=isograph-scale-test
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=5
#SBATCH --gpus=v100-32:1
#SBATCH --time=04:00:00
#SBATCH --output=develop/_m/logs/%x-%j.log
# IsoGraph scale test: troubleshoot parameter choices at intermediate and real
# data scale before committing to publication-quality benchmarks.
# Runs both allow_abundance_abundance=False and =True variants at 5k and 17k genes.
# Results written to develop/_m/scale_test_results.parquet.
#
# Override options (passed directly to scale_test.py):
#   sbatch scale_test.sh --scale intermediate   # skip large scale
#   sbatch scale_test.sh --leiden-resolution 3.0
#   sbatch scale_test.sh --replicates 1         # single quick pass
set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "ERROR: project root does not exist: ${PROJECT_ROOT}"
    exit 1
fi
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the isograph-brain-aging-benchmarking repo root or set ISOGRAPH_BENCHMARK_ROOT."
    echo "Current project root candidate: ${PROJECT_ROOT}"
    exit 1
fi
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p develop/_m/logs

log_message "**** IsoGraph scale test starts ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID:-local}"
echo "Node: ${SLURM_NODENAME:-local}"

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "CUDA diagnostics"
nvidia-smi || echo "nvidia-smi not found or no GPU"
python - <<'PYEOF'
import torch
print(f"PyTorch version  : {torch.__version__}")
print(f"CUDA available   : {torch.cuda.is_available()}")
if torch.cuda.is_available():
    for i in range(torch.cuda.device_count()):
        props = torch.cuda.get_device_properties(i)
        print(f"  GPU {i}: {torch.cuda.get_device_name(i)}, {props.total_memory/1e9:.1f} GB")
PYEOF

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

python develop/scale_test.py "$@"

conda deactivate
log_message "**** IsoGraph scale test ends ****"
