#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=brainseq-iso-aging
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=32
#SBATCH --array=1-3
#SBATCH --time=02:00:00
#SBATCH --output=02_module_discovery/_m/logs/%x-%A_%a.log
# Run IsoGraph VAE on BrainSEQ adult control aging bundles:
# caudate, hippocampus, and DLPFC.
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
mkdir -p 02_module_discovery/_m/logs

log_message "**** BrainSEQ aging IsoGraph job starts ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID:-local}"
echo "Job name: ${SLURM_JOB_NAME:-brainseq-iso-aging}"
echo "Node name: ${SLURM_NODENAME:-local}"
echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "Activating IsoGraph environment"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

log_message "Checking Python analysis dependencies"
python -c "import isograph_benchmark, numpy, pandas, scipy, patsy"

log_message "CUDA diagnostics"
nvidia-smi || echo "nvidia-smi not found or no GPU"
python - <<'PYEOF'
import torch
print(f"PyTorch version  : {torch.__version__}")
print(f"PyTorch CUDA ver : {torch.version.cuda}")
print(f"CUDA available   : {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"GPU count        : {torch.cuda.device_count()}")
    for i in range(torch.cuda.device_count()):
        props = torch.cuda.get_device_properties(i)
        print(f"  GPU {i}: {torch.cuda.get_device_name(i)}, {props.total_memory/1e9:.1f} GB")
    try:
        torch.zeros(1).cuda()
        print("CUDA smoke test  : PASSED")
    except Exception as e:
        print(f"CUDA smoke test  : FAILED — {e}")
else:
    print("WARNING: CUDA not available — check driver/toolkit compatibility")
PYEOF

REGIONS=(caudate hippocampus dlpfc)
RUN_ARGS=("$@")
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    REGION="${REGIONS[$((SLURM_ARRAY_TASK_ID - 1))]}"
    RUN_ARGS=(--region "${REGION}" "$@")
    log_message "Running region ${REGION}"
fi

# Limit BLAS/OpenMP threads to avoid per-thread workspace accumulation
# over the 17k+ SVD calls in gene_switch_coordinates (main CPU OOM cause).
# 8 threads balances memory footprint with training throughput on RM-shared.
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

python -m isograph_benchmark.real_data.run_models brainseq-aging "${RUN_ARGS[@]}"

conda deactivate
log_message "**** BrainSEQ aging IsoGraph job ends ****"
