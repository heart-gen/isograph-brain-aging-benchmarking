#!/usr/bin/env bash

set -euo pipefail

source ~/.venvs/isograph/bin/activate

python -m isograph_benchmark.benchmark.run_synthetic
python -m isograph_benchmark.benchmark.make_batches

batch_dir="benchmark/01_synthetic/_m"
script="benchmark/01_synthetic/_h/run_batch.sh"

for batch_file in "${batch_dir}"/batches_*.tsv; do
  resource_class="$(basename "${batch_file}" .tsv)"
  resource_class="${resource_class#batches_}"
  n_tasks="$(($(wc -l < "${batch_file}") - 1))"
  if [[ "${n_tasks}" -le 0 ]]; then
    continue
  fi

  case "${resource_class}" in
    cpu_short)
      sbatch_opts=(--time=02:00:00 --cpus-per-task=16 --array="1-${n_tasks}%100")
      ;;
    vae)
      sbatch_opts=(--time=02:00:00 --cpus-per-task=32 --array="1-${n_tasks}%50")
      ;;
    wgcna_cpu)
      sbatch_opts=(--time=02:00:00 --cpus-per-task=50 --array="1-${n_tasks}%80")
      ;;
    scale)
      sbatch_opts=(--time=04:00:00 --cpus-per-task=64 --array="1-${n_tasks}%20")
      ;;
    *)
      sbatch_opts=(--time=02:00:00 --cpus-per-task=32 --array="1-${n_tasks}%50")
      ;;
  esac

  echo "BATCH_FILE=${batch_file} sbatch ${sbatch_opts[*]} ${script}"
done
