#!/usr/bin/env bash
# Stage the heavy regenerable real-data artifacts for Zenodo upload.
#
# These artifacts are git-ignored (too large for the GitHub/LFS repo) and are
# distributed via Zenodo instead. This writes zenodo/MANIFEST.tsv (relative path,
# bytes [, sha256]) listing exactly what to upload, and optionally a single tar.
#
# Usage:
#   00_scripts/stage_zenodo_bundle.sh                 # manifest with sizes (fast)
#   00_scripts/stage_zenodo_bundle.sh --checksums     # add sha256 (slower)
#   00_scripts/stage_zenodo_bundle.sh --tar           # also build the upload tarball
set -euo pipefail

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${PWD}}"
cd "${PROJECT_ROOT}"
[[ -f .here ]] || { echo "ERROR: run from the repo root."; exit 1; }

CHECKSUMS=0; MAKE_TAR=0
for a in "$@"; do
    case "$a" in
        --checksums) CHECKSUMS=1 ;;
        --tar) MAKE_TAR=1 ;;
        *) echo "unknown arg: $a"; exit 1 ;;
    esac
done

# Heavy artifacts -> Zenodo (keep in sync with .gitignore).
NAMES=(feature_scores.parquet feature_reconstruction.parquet \
       high_vs_low_table.parquet edges.parquet)

mkdir -p zenodo
MANIFEST=zenodo/MANIFEST.tsv
if [[ "${CHECKSUMS}" == "1" ]]; then
    printf 'path\tbytes\tsha256\n' > "${MANIFEST}"
else
    printf 'path\tbytes\n' > "${MANIFEST}"
fi

mapfile -t FILES < <(
    for n in "${NAMES[@]}"; do
        # real-data stages 02-07 (the synthetic benchmark in 01 is archived separately)
        find 02_module_discovery 04_module_trust 03_module_characterization \
             05_genetic_anchoring 06_switch_mechanism 07_rbp_regulation \
             -type f -name "${n}" -not -path '*/_m/tmp/*' -not -path '*/logs/*'
    done | sort
)

total=0
for f in "${FILES[@]}"; do
    bytes=$(stat -c '%s' "${f}")
    total=$((total + bytes))
    if [[ "${CHECKSUMS}" == "1" ]]; then
        sha=$(sha256sum "${f}" | cut -d' ' -f1)
        printf '%s\t%s\t%s\n' "${f}" "${bytes}" "${sha}" >> "${MANIFEST}"
    else
        printf '%s\t%s\n' "${f}" "${bytes}" >> "${MANIFEST}"
    fi
done

printf 'Staged %d files, %.2f GB -> %s\n' "${#FILES[@]}" \
    "$(echo "${total}/1073741824" | bc -l)" "${MANIFEST}"

if [[ "${MAKE_TAR}" == "1" ]]; then
    TAR=zenodo/isograph_realdata_heavy_artifacts.tar
    tar -cf "${TAR}" "${FILES[@]}"
    printf 'Wrote %s (%.2f GB)\n' "${TAR}" "$(echo "$(stat -c '%s' "${TAR}")/1073741824" | bc -l)"
fi
