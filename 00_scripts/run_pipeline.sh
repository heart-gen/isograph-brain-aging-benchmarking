#!/usr/bin/env bash
## Submit the real-data analysis, stage by stage. Each stage's root steps wait on every job the
## previous stage submitted, and each stage's internal order comes from <stage>/_h/run_stage.sh.
##
##   bash run_pipeline.sh --dry-run                        # print every stage's plan
##   bash run_pipeline.sh                                  # stages 02-08
##   bash run_pipeline.sh --stages 05,06,07,08 --after <job ids of work already queued>
##   bash run_pipeline.sh --login-done 06:01f,07:03b       # login-node steps already run
##
## A stage with held steps (a login-node step not yet done) ends the chain there: later stages
## read those outputs. Run the login step, re-run that stage with --login-done, then continue
## with --stages <remaining>.
##
## Not covered: inputs/ (inputs/_h/build_data_pipeline.sh) and the archived synthetic benchmark in
## 01_synthetic_benchmark. Memory on Bridges-2 is --cpus-per-task x 2000MB throughout.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
[[ -f .here ]] || { echo "ERROR: run from the repo"; exit 1; }

STAGES=(02_module_discovery 03_module_characterization 04_module_trust 05_genetic_anchoring
        06_switch_mechanism 07_rbp_regulation 08_integration)
DRY=(); AFTER=""; WANT=""; LOGIN=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY=(--dry-run) ;;
        --after) AFTER=$2; shift ;;
        --stages) WANT=$2; shift ;;
        --login-done) LOGIN=$2; shift ;;
        -h|--help) sed -n '2,17p' "$0"; exit 0 ;;
        *) echo "unknown option $1"; exit 1 ;;
    esac
    shift
done

for s in "${STAGES[@]}"; do
    n=${s:0:2}
    if [[ -n "${WANT}" && ",${WANT}," != *",${n},"* ]]; then continue; fi
    args=(${DRY[@]+"${DRY[@]}"})
    [[ -n "${AFTER}" ]] && args+=(--after "${AFTER}")
    done_ids=$(tr ',' '\n' <<< "${LOGIN}" | sed -n "s/^${n}://p" | paste -sd, -)
    [[ -n "${done_ids}" ]] && args+=(--login-done "${done_ids}")

    out=$(bash "${s}/_h/run_stage.sh" ${args[@]+"${args[@]}"})
    jobs=$(sed -n 's/^STAGE_JOBS=//p' <<< "${out}" | tail -1)
    held=$(sed -n 's/^STAGE_HELD=//p' <<< "${out}" | tail -1)
    [[ -n "${jobs}" ]] && AFTER="${jobs}"
    if [[ "${held:-0}" != "0" ]]; then
        echo "== ${s} has ${held} held step(s); stopping the chain here (see the LOGIN lines above)." >&2
        exit 0
    fi
done
echo "== pipeline submitted${DRY:+ (dry run)}" >&2
