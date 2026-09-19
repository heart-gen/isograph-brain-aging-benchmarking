#!/usr/bin/env bash
# Minimal SLURM dependency-graph submitter, sourced by every <stage>/_h/run_stage.sh.
#
# A stage runner declares its steps in run order; each step names the step ids it waits on, and
# the library turns that into `sbatch --dependency=afterok:...`. Step ids are the wrapper's tier
# prefix (`03d`), optionally with a variant suffix (`03d.complement`); a dependency on `03d`
# means every submitted variant of it.
#
#   source scripts/slurm_dag.sh
#   dag_init 04_module_trust "$@"
#   step 01a ""        04_module_trust/_h/01a.stability_isograph.sh
#   step 02c "01a 01b" 04_module_trust/_h/02c.trust_gate.sh
#   step 03d.full "02c" --export=ALL,FOO=1 04_module_trust/_h/03d.replication_permutation.sh --covariates full
#   login_step 03b "02a" 07_rbp_regulation/_h/03b.rbp_binding_fetch.sh
#   manual_step 08d "per-locus; see the stage README"
#   dag_finish
#
# sbatch options go between the dependency list and the script and must use the `--opt=value`
# form (the first argument not starting with `-` is the script).
#
# Runner options (all optional):
#   --dry-run              print the sbatch commands, submit nothing
#   --after JOBIDS         colon-separated job ids every root step waits on (stage chaining)
#   --from TIER --to TIER  submit only tiers in [TIER, TIER] (e.g. --from 03); a dependency on a
#                          step outside the range is treated as already satisfied
#   --only IDS / --skip IDS  comma-separated step ids or id prefixes
#   --login-done IDS       login-node steps already run; steps waiting on an unfinished one are
#                          held (listed, not submitted) until the runner is re-run with --from
#   --ledger FILE          where to append "<id>\t<label>\t<jobid>" (default: <stage>/_m/logs/)
# The last stdout lines are `STAGE_HELD=<n>` and `STAGE_JOBS=<id:id:...>` (the stage's leaf jobs,
# i.e. the ones nothing else in the stage waits on) for run_pipeline.sh.

set -euo pipefail

DAG_STAGE=""; DAG_DRY=0; DAG_AFTER=""; DAG_FROM=0; DAG_TO=99
DAG_ONLY=""; DAG_SKIP=""; DAG_LOGIN_DONE=""; DAG_LEDGER=""
declare -A DAG_JOB=()      # step id -> job id (submitted this run)
declare -A DAG_HELD=()     # step id -> reason it was not submitted
declare -A DAG_USED=()     # step id -> 1 once a later step waits on it (non-leaves)
DAG_ALL=()

_dag_die() { echo "slurm_dag: $*" >&2; exit 1; }

dag_init() {
    DAG_STAGE=$1; shift
    local here; here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    cd "${ISOGRAPH_BENCHMARK_ROOT:-${here}}"
    [[ -f .here && -d "${DAG_STAGE}/_h" ]] || _dag_die "run from the repo (no ${DAG_STAGE}/_h here)"
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --dry-run) DAG_DRY=1 ;;
            --after) DAG_AFTER=$2; shift ;;
            --from) DAG_FROM=$((10#$2)); shift ;;
            --to) DAG_TO=$((10#$2)); shift ;;
            --only) DAG_ONLY=$2; shift ;;
            --skip) DAG_SKIP=$2; shift ;;
            --login-done) DAG_LOGIN_DONE=$2; shift ;;
            --ledger) DAG_LEDGER=$2; shift ;;
            -h|--help) sed -n '2,30p' "${BASH_SOURCE[0]}"; exit 0 ;;
            *) _dag_die "unknown option $1" ;;
        esac
        shift
    done
    if [[ -z "${DAG_LEDGER}" ]]; then
        mkdir -p "${DAG_STAGE}/_m/logs"
        DAG_LEDGER="${DAG_STAGE}/_m/logs/run_stage-$(date +%Y%m%d-%H%M%S).tsv"
    fi
    [[ ${DAG_DRY} == 1 ]] && DAG_LEDGER=/dev/null
    echo "== ${DAG_STAGE}: tiers ${DAG_FROM}-${DAG_TO}${DAG_AFTER:+, after ${DAG_AFTER}}${DAG_DRY:+}" >&2
}

_dag_in_list() {  # _dag_in_list <id> <comma list>: exact id or id prefix match
    local id=$1 x; local IFS=,
    for x in $2; do [[ -n "$x" && ( "$id" == "$x" || "$id" == "$x".* || "$id" == "$x"* ) ]] && return 0; done
    return 1
}

_dag_selected() {
    local id=$1 tier=$((10#${1:0:2}))
    (( tier >= DAG_FROM && tier <= DAG_TO )) || return 1
    [[ -n "${DAG_ONLY}" ]] && ! _dag_in_list "$id" "${DAG_ONLY}" && return 1
    [[ -n "${DAG_SKIP}" ]] && _dag_in_list "$id" "${DAG_SKIP}" && return 1
    return 0
}

# job ids for a space-separated dependency list; returns 2 if any dependency is held
_dag_deps() {
    local d k out=()
    for d in $1; do
        for k in "${!DAG_HELD[@]}"; do
            [[ "$k" == "$d" || "$k" == "$d".* ]] && return 2
        done
        for k in "${!DAG_JOB[@]}"; do
            [[ "$k" == "$d" || "$k" == "$d".* ]] && out+=("${DAG_JOB[$k]}")
        done
    done
    local IFS=:; echo "${out[*]}"
}

step() {
    local id=$1 deps=$2; shift 2
    local sargs=()
    while [[ $# -gt 0 && "$1" == -* ]]; do sargs+=("$1"); shift; done
    [[ $# -gt 0 ]] || _dag_die "step ${id}: no script"
    local script=$1; shift
    [[ -f "${script}" ]] || _dag_die "step ${id}: missing ${script}"
    _dag_selected "${id}" || return 0
    local jobs rc=0
    jobs=$(_dag_deps "${deps}") || rc=$?
    if [[ ${rc} == 2 ]]; then
        DAG_HELD[$id]="waits on a held or login-node step"
        echo "   held ${id}  $(basename "${script}") $*" >&2
        return 0
    fi
    local d k
    for d in ${deps}; do     # marked here, not in _dag_deps: that runs in a subshell
        for k in "${!DAG_JOB[@]}"; do
            [[ "$k" == "$d" || "$k" == "$d".* ]] && DAG_USED[$k]=1
        done
    done
    [[ -z "${jobs}" ]] && jobs="${DAG_AFTER}"
    local cmd=(sbatch --parsable)
    [[ -n "${jobs}" ]] && cmd+=(--kill-on-invalid-dep=yes "--dependency=afterok:${jobs}")
    cmd+=(${sargs[@]+"${sargs[@]}"} "${script}" "$@")
    local jid
    if [[ ${DAG_DRY} == 1 ]]; then
        jid="dry${id}"
        echo "   ${id}  ${cmd[*]}" >&2
    else
        jid=$("${cmd[@]}") || _dag_die "sbatch failed for ${id}"
        jid=${jid%%;*}
        echo "   ${id}  ${jid}  $(basename "${script}") $*" >&2
    fi
    DAG_JOB[$id]=${jid}
    DAG_ALL+=("${jid}")
    printf '%s\t%s\t%s\n' "${id}" "$(basename "${script}") $*" "${jid}" >> "${DAG_LEDGER}"
}

login_step() {  # login_step <id> <deps> <script> [args]: needs network; never submitted
    local id=$1 deps=$2 script=$3; shift 3
    [[ -f "${script}" ]] || _dag_die "login step ${id}: missing ${script}"
    _dag_selected "${id}" || return 0
    if _dag_in_list "${id}" "${DAG_LOGIN_DONE}"; then
        echo "   done ${id}  (login step, reported done)" >&2
        return 0
    fi
    DAG_HELD[$id]="login-node step"
    echo "   LOGIN ${id}: after [${deps:-nothing}] finishes, run on a login node:" >&2
    echo "          bash ${script} $*" >&2
    echo "          then re-run this stage with --from <next tier> --login-done ${id}" >&2
}

manual_step() {  # manual_step <id> <text>: documented, never submitted
    _dag_selected "$1" || return 0
    echo "   manual $1: $2" >&2
}

dag_finish() {
    local n=${#DAG_ALL[@]} h=${#DAG_HELD[@]}
    echo "== ${DAG_STAGE}: ${n} submitted, ${h} held; ledger ${DAG_LEDGER}" >&2
    echo "STAGE_HELD=${h}"
    # Only the leaves: every other job is upstream of one of them, so a later stage that waits on
    # the leaves waits on the whole stage, with a dependency list short enough for sbatch.
    local id leaves=()
    for id in "${!DAG_JOB[@]}"; do
        [[ -z "${DAG_USED[$id]:-}" ]] && leaves+=("${DAG_JOB[$id]}")
    done
    local IFS=:
    echo "STAGE_JOBS=${leaves[*]}"
}
