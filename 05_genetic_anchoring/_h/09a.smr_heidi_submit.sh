#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=smr-heidi-submit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/%x-%j.log
## SMR + HEIDI, end to end for one QTL source, as one submittable unit.
##
## The 10a.smr_heidi.sh array size is the length of work_list.tsv, which does not exist until
## `smr_heidi --stage prep` has read the current nominations. So this job runs prep, sizes the
## array from the work list it just wrote, and submits the rest of the chain itself:
##
##   10a primary (5e-8, single-SNP, builds the BESD)  -> 11a meta
##   10a --peqtl-smr 1e-6 (reuses the BESD)          -> 11a meta --peqtl-smr 1e-6
##   10a --smr-multi      (reuses the BESD)          -> 11a meta --smr-multi
##
## Source is env-selected, the same variables 10a reads:
##   GTEx      sbatch --dependency=afterok:<coloc_signal meta --sqtl all> 05_genetic_anchoring/_h/09a.smr_heidi_submit.sh
##   BrainSEQ  sbatch --dependency=afterok:<coloc_brainseq meta>:<28 ea_only> \
##                    --export=ALL,SMR_QTL_SOURCE=brainseq,SMR_ARM=ea_only 05_genetic_anchoring/_h/09a.smr_heidi_submit.sh
##
## Set SMR_SUBMIT_LEDGER=<file> to append "<label>\t<jobid>" for every job this submits, so a
## wave ledger can account for jobs that no outer launcher saw.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 05_genetic_anchoring/_m/logs
command -v sbatch >/dev/null || { echo "ERROR: sbatch not on PATH in this job"; exit 1; }

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
H=05_genetic_anchoring/_h
SRC="${SMR_QTL_SOURCE:-gtex}"
ARM="${SMR_ARM:-}"
SRC_ARGS=(--qtl-source "${SRC}")
[[ -n "${ARM}" ]] && SRC_ARGS+=(--arm "${ARM}")
if [[ "${SRC}" == "gtex" ]]; then
    WORK=05_genetic_anchoring/_m/smr_heidi/gtex/work_list.tsv
    META_SRC=()
else
    [[ -n "${ARM}" ]] || { echo "ERROR: SMR_ARM is required for ${SRC}"; exit 1; }
    WORK="05_genetic_anchoring/_m/smr_heidi/${SRC}/${ARM}/work_list.tsv"
    META_SRC=(--qtl-source "${SRC}" --arm "${ARM}")
fi

log "**** smr_heidi --stage prep ${SRC_ARGS[*]} ****"
"${PY}" -u -m isograph_benchmark.real_data.smr_heidi --stage prep "${SRC_ARGS[@]}"

[[ -f "${WORK}" ]] || { echo "ERROR: prep wrote no ${WORK}"; exit 1; }
# header + data rows; 10a's own usage note sizes the array as (lines - 2)
N=$(( $(wc -l < "${WORK}") - 2 ))
(( N >= 0 )) || { echo "ERROR: empty work list ${WORK}"; exit 1; }
log "work list: $((N + 1)) tasks"

EXP="ALL,SMR_QTL_SOURCE=${SRC}${ARM:+,SMR_ARM=${ARM}}"
TAG="${SRC}${ARM:+_${ARM}}"
sub() {  # sub <label> <sbatch args...>; echoes the job id
    local label=$1; shift
    local jid; jid=$(sbatch --parsable --kill-on-invalid-dep=yes "$@"); jid=${jid%%;*}
    [[ -n "${SMR_SUBMIT_LEDGER:-}" ]] && printf '%s\t%s\n' "${label}" "${jid}" >> "${SMR_SUBMIT_LEDGER}"
    log "submitted ${label}: ${jid}" >&2
    echo "${jid}"
}

P=$(sub "smr_${TAG}" --array=0-${N} --export="${EXP}" "${H}/10a.smr_heidi.sh")
sub "smr_meta_${TAG}" --dependency=afterok:${P} "${H}/11a.smr_heidi_meta.sh" ${META_SRC[@]+"${META_SRC[@]}"} >/dev/null

R=$(sub "smr_peqtl1e-6_${TAG}" --dependency=afterok:${P} --array=0-${N} \
        --export="${EXP},SMR_PEQTL=1e-6,SMR_SKIP_BESD=1" "${H}/10a.smr_heidi.sh")
sub "smr_meta_peqtl1e-6_${TAG}" --dependency=afterok:${R} "${H}/11a.smr_heidi_meta.sh" \
    --peqtl-smr 1e-6 ${META_SRC[@]+"${META_SRC[@]}"} >/dev/null

U=$(sub "smr_multi_${TAG}" --dependency=afterok:${P} --array=0-${N} \
        --export="${EXP},SMR_MULTI=1,SMR_SKIP_BESD=1" "${H}/10a.smr_heidi.sh")
sub "smr_meta_multi_${TAG}" --dependency=afterok:${U} "${H}/11a.smr_heidi_meta.sh" \
    --smr-multi ${META_SRC[@]+"${META_SRC[@]}"} >/dev/null

log "**** SMR chain submitted for ${TAG} ****"
