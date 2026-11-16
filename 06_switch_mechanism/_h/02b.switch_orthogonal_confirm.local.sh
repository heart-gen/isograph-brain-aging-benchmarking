#!/usr/bin/env bash
## Long-read switch-like rates against their matched nulls (local; no scheduler).
##
## Reruns the two 02b arms the manuscript reports, with the core packages pinned to
## env/environment.lock.txt so the draws are reproducible off the cluster:
##   --mode anchored --events signal   colocalization-anchored pairs (signal_coloc/; Data S21,
##                                     Fig. 6a anchored rows)
##   --mode global-null                every IsoGraph switch pair (Fig. 6a genome-wide row)
## Each arm reports the abundance-decile null (`matched_null`) and the stricter joint null
## (`matched_null_joint`: abundance quintile x detection-rate, isoform-count and usage-SD
## tertiles; `longread_null.JOINT_KEYS`/`JOINT_BINS`), 2,000 draws, seed 13.
##
## Run from the repository root (about 1 min for the anchored arm; the global null draws
## ~16,000 pairs per draw and takes longer):
##   bash 06_switch_mechanism/_h/02b.switch_orthogonal_confirm.local.sh            # both arms
##   ARMS=anchored bash 06_switch_mechanism/_h/02b.switch_orthogonal_confirm.local.sh
## Needs uv and the Bambu matrix inputs/_m/longread_aged_dlpfc/counts_transcript.txt.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${PWD}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: run from the repository root."; exit 1; }
mkdir -p 06_switch_mechanism/_m/logs
LOG="06_switch_mechanism/_m/logs/switch-orthogonal-local-$(date '+%Y%m%d-%H%M%S').log"
exec > >(tee -a "${LOG}") 2>&1

export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

## Versions from env/environment.lock.txt. The draws depend on numpy's generator and on
## pandas' grouping order, so these are pinned rather than taken from requirements ranges.
command -v uv >/dev/null || { echo "ERROR: uv is required."; exit 1; }
PY=(uv run --quiet --no-project --python 3.12
    --with isograph==0.1.5 --with numpy==2.4.4 --with pandas==2.3.3
    --with pyarrow==24.0.0 --with scipy==1.17.1 --with pyhere --with statsmodels==0.14.6
    python)

COUNTS=inputs/_m/longread_aged_dlpfc/counts_transcript.txt
if [[ ! -s "${COUNTS}" ]] || head -c 40 "${COUNTS}" | grep -q 'git-lfs'; then
    log "restoring ${COUNTS}"
    git lfs pull --include="${COUNTS}" || \
        "${PY[@]}" -m isograph_benchmark.real_data.longread_switch_confirm fetch
fi

MOD=isograph_benchmark.real_data.switch_orthogonal_confirm
OC=06_switch_mechanism/_m/switch_orthogonal_confirm
for ARM in ${ARMS:-anchored global-null}; do
    case "${ARM}" in
        anchored)    log "anchored arm (signal layer)"
                     "${PY[@]}" -m "${MOD}" --mode anchored --events signal > /dev/null ;;
        global-null) log "global null"
                     "${PY[@]}" -m "${MOD}" --mode global-null > /dev/null ;;
        *)           echo "ERROR: unknown arm ${ARM}"; exit 1 ;;
    esac
done

## Headline numbers, read back from the summaries (distinct pairs are the primary unit).
"${PY[@]}" - "${OC}" <<'EOF'
import json, sys
from pathlib import Path
oc = Path(sys.argv[1])
rows = []
s = oc / "signal_coloc" / "anchored_summary.json"
if s.exists():
    d = json.loads(s.read_text())["distinct_pairs"]
    for key, label in (("switch_like_rate", "anchored, all detected"),
                       ("switch_like_rate_usable_only", "anchored, well expressed")):
        rows.append((label, d["matched_null"][key], d["matched_null_joint"][key]))
g = oc / "global_null_summary.json"
if g.exists():
    d = json.loads(g.read_text())
    if "matched_null_joint" in d:
        rows.append(("genome-wide switch pairs", d["matched_null"]["switch_like_rate"],
                     d["matched_null_joint"]["switch_like_rate"]))
for label, ab, jt in rows:
    print(f"{label:26s} n={jt['n_focal']:<6d} observed {jt['observed']:.3f} | "
          f"abundance null {ab['null_mean']:.3f} (P = {ab['p_empirical_two_sided']:.2g}) | "
          f"joint null {jt['null_mean']:.3f} [{jt['null_q025']:.3f}-{jt['null_q975']:.3f}] "
          f"(P = {jt['p_empirical_two_sided']:.2g})")
EOF
log "outputs: ${OC}/{signal_coloc/,global_null_summary.json}"
log "**** complete ****"
