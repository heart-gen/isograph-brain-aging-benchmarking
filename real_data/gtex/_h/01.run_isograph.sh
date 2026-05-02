#!/usr/bin/env bash
# Run IsoGraph VAE on all 13 GTEx v11 brain aging regions.
# Outputs land in: real_data/gtex/<region>/_m/isograph_vae/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"

python3 -m isograph_benchmark.real_data.run_models gtex-aging "$@"
