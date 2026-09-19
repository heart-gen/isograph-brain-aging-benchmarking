#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ldsc-h2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/ldsc-h2-%j.log

## S-LDSC step 3 — munge one GWAS trait + joint partitioned heritability with the
## IsoGraph sQTL/eQTL/cis switch annotations on top of baselineLD v2.2. The joint
## model gives the conditional per-SNP heritability (tau) of the sQTL annotation
## beyond eQTL + cis + baseline = splicing-specific heritability in the switch layer.
##
## Two cases share this script:
##   disease : sbatch 05_genetic_anchoring/_h/03d.ldsc_munge_h2.sh scz brainseq-sczd
##   aging   : sbatch 05_genetic_anchoring/_h/03d.ldsc_munge_h2.sh ad  aging   (also pd / lbd / als / ftd)
##
## Trait column specs come from isograph_benchmark.real_data.gwas_traits (single source
## of truth). Munge is build-agnostic (HapMap3 rsID merge), so hg38 neurodegeneration
## sumstats need no liftover. Adapts ancestry-aging 14_gwas_ldsc steps 2 + 5 to PSC.
##
## Usage: sbatch 05_genetic_anchoring/_h/03d.ldsc_munge_h2.sh <trait> <annotation>
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

TRAIT="${1:-scz}"
ANNOT="${2:-brainseq-sczd}"
PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
mkdir -p 05_genetic_anchoring/_m/logs

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi

LDSC_DIR=/ocean/projects/bio250020p/shared/opt/ldsc
RES=/ocean/projects/bio250020p/shared/resources/ldsc
BASELINE="${RES}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD."
WEIGHTS="${RES}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
FRQ="${RES}/1000G_Phase3_frq/1000G.EUR.QC."
HM3_SNPLIST="${RES}/w_hm3.snplist"
WRAP="05_genetic_anchoring/_h/ldsc_wrapper.py"

MDIR="${PROJECT_ROOT}/05_genetic_anchoring/_m/ldsc/${ANNOT}"
SUMD="${MDIR}/sumstats/${TRAIT}"; RESD="${MDIR}/results"; mkdir -p "${SUMD}" "${RESD}"

## --- phase A: resolve munge args from the trait registry (isograph env) --------
module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph
log "Resolving munge args for trait '${TRAIT}'"
mapfile -t MUNGE_ARGS < <(python -m isograph_benchmark.real_data.gwas_traits \
    munge-args --trait "${TRAIT}" --clean-dir "${SUMD}")
[[ ${#MUNGE_ARGS[@]} -gt 0 ]] || { echo "ERROR: no munge args for ${TRAIT}"; exit 1; }
conda deactivate

## --- phase B: munge + partitioned h2 (ldsc genomics env) -----------------------
module load bedtools/2.30.0
conda activate /ocean/projects/bio250020p/shared/opt/env/genomics

if [[ ! -f "${SUMD}/${TRAIT}.sumstats.gz" ]]; then
    log "Munging ${TRAIT}"
    python "${WRAP}" "${LDSC_DIR}" munge_sumstats.py \
        "${MUNGE_ARGS[@]}" --out "${SUMD}/${TRAIT}" \
        --merge-alleles "${HM3_SNPLIST}" --chunksize 500000
fi

lp() { echo "${MDIR}/ldscores/${1}/${1}."; }
for A in sqtl_switch eqtl_switch cis_switch; do
    [[ -f "$(lp "${A}")1.l2.ldscore.gz" ]] || { echo "ERROR: LD scores missing for ${A} in ${ANNOT}; run step 2."; exit 1; }
done

## Four S-LDSC models, each on top of baselineLD v2.2:
##   sqtl_only / eqtl_only / cis_only — single-annotation Finucane enrichment + tau
##     conditional on baselineLD (the primary "is the switch-QTL layer genetically
##     anchored" test; size-robust, unlike MAGMA).
##   joint — baseline + sqtl_switch + eqtl_switch together, for the head-to-head
##     splicing-vs-expression conditional contrast. (cis_switch is the sqtl u eqtl
##     union, so it is never mixed WITH its components — that would be collinear.)
run_h2() {
    local model="$1" ref="$2"
    log "h2 model=${model} (${TRAIT}, annot=${ANNOT})"
    python "${WRAP}" "${LDSC_DIR}" ldsc.py \
        --h2 "${SUMD}/${TRAIT}.sumstats.gz" \
        --ref-ld-chr "${ref}" \
        --w-ld-chr "${WEIGHTS}" \
        --frqfile-chr "${FRQ}" \
        --overlap-annot --thin-annot --print-coefficients \
        --out "${RESD}/${TRAIT}_${model}"
}
run_h2 sqtl_only "${BASELINE},$(lp sqtl_switch)"
run_h2 eqtl_only "${BASELINE},$(lp eqtl_switch)"
run_h2 cis_only  "${BASELINE},$(lp cis_switch)"
run_h2 joint     "${BASELINE},$(lp sqtl_switch),$(lp eqtl_switch)"

conda deactivate
log "**** S-LDSC h2 done: ${RESD}/${TRAIT}_{sqtl_only,eqtl_only,cis_only,joint}.results ****"
