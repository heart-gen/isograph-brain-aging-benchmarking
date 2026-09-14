#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=GPU-shared
#SBATCH --job-name=gpu-repro-probe
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:1
#SBATCH --array=1-24
#SBATCH --time=02:00:00
#SBATCH --output=01_synthetic_benchmark/01_synthetic/_m/logs/%x-%A_%a.log

## GPU VAE reproducibility probe: re-run one stored isograph_vae_gpu run twice on GPU and its
## isograph_vae CPU twin once, from the current grid config, into ISOLATED output roots under
## 01_synthetic/_o/gpu_repro_probe/{gpu_rep1,gpu_rep2,cpu_rep1}/. The stored runs are never
## touched. Rows come from `gpu_reproducibility_probe select`, which only admits runs whose
## cached dataset exists -- run_one would otherwise regenerate it, and the synthetic datasets
## are archived, not regenerable. Then: `gpu_reproducibility_probe compare`.
##
## Calls the env interpreter directly (no `module load`), so a task cannot die on
## "module: command not found"; the torch cu126 wheel ships its own CUDA runtime.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 01_synthetic_benchmark/01_synthetic/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
PROBE=01_synthetic_benchmark/03_metrics/_m/gpu_repro_probe/probe_runs.tsv
ROOT=01_synthetic_benchmark/01_synthetic/_o/gpu_repro_probe
IDX="${SLURM_ARRAY_TASK_ID:-1}"

row=$(awk -F'\t' -v i="${IDX}" 'NR>1 && $1==i' "${PROBE}")
[[ -n "${row}" ]] || { echo "ERROR: no probe row ${IDX} in ${PROBE}"; exit 1; }
RUN_GPU=$(echo "${row}" | cut -f6)
RUN_CPU=$(echo "${row}" | cut -f7)
log "probe ${IDX}: $(echo "${row}" | cut -f2) gpu=${RUN_GPU} cpu=${RUN_CPU}"

for rep in gpu_rep1 gpu_rep2; do
    "${PY}" -m isograph_benchmark.benchmark.run_one --run-id "${RUN_GPU}" \
        --output-root "${ROOT}/${rep}" --force
done
"${PY}" -m isograph_benchmark.benchmark.run_one --run-id "${RUN_CPU}" \
    --output-root "${ROOT}/cpu_rep1" --force
log "**** complete ****"
