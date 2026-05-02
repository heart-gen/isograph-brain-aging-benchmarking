#!/usr/bin/env bash
# Run IsoGraph VAE on all 13 GTEx v11 brain regions.
# Each region is independent; set PARALLEL=1 to launch background jobs.
# Outputs land in: real_data/gtex/<region>/_m/isograph_vae/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

PARALLEL="${PARALLEL:-0}"
PIDS=()

REGIONS=(
    amygdala
    anterior_cingulate_cortex_ba24
    caudate_basal_ganglia
    cerebellar_hemisphere
    cerebellum
    cortex
    frontal_cortex_ba9
    hippocampus
    hypothalamus
    nucleus_accumbens_basal_ganglia
    putamen_basal_ganglia
    spinal_cord_cervical_c_1
    substantia_nigra
)

run_region() {
    local region="$1"
    echo "[$(date +%H:%M:%S)] Starting $region"
    python3 - <<EOF
import sys
sys.path.insert(0, "${PROJECT_ROOT}")
from isograph_benchmark.real_data.run_models import run_gtex_region
run_gtex_region("${region}")
EOF
    echo "[$(date +%H:%M:%S)] Done    $region"
}

for region in "${REGIONS[@]}"; do
    if [[ "${PARALLEL}" == "1" ]]; then
        run_region "$region" &
        PIDS+=($!)
    else
        run_region "$region"
    fi
done

if [[ "${PARALLEL}" == "1" ]]; then
    for pid in "${PIDS[@]}"; do
        wait "$pid"
    done
fi

echo "All GTEx regions complete."
