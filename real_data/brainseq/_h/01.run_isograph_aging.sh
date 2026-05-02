#!/usr/bin/env bash
# Run IsoGraph VAE on BrainSEQ adult control aging bundles:
# caudate, hippocampus, and DLPFC.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"

python3 -m isograph_benchmark.real_data.run_models brainseq-aging "$@"
