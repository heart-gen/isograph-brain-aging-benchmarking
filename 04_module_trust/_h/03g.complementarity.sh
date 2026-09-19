#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=trust-complementarity
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --array=1-12
#SBATCH --time=01:00:00
#SBATCH --output=04_module_trust/_m/logs/trust-complementarity-%A_%a.log

## Q4 complementarity (module_trust `complementarity` subcommand): for each trusted module,
## the biology an abundance method cannot see -- DTU-without-DGE gene fraction, overlap with
## age-associated classical WGCNA modules, and the structural switches of its driver
## transcripts.
##
## Had no committed launcher. It joins four upstream outputs, so run it only after all of:
##   04_module_trust/_h/02c.trust_gate.sh                                 (trusted sets)
##   03_module_characterization/_h/01a-01b.interpret_modules_*.sh          (drivers, structure_annotations)
##   03_module_characterization/_h/02b.characterize_composition_unique.sh (composition_unique genes)
## e.g.
##   sbatch --dependency=afterok:${gate}:${interp_b}:${interp_g}:${comp_unique} 04_module_trust/_h/03g.complementarity.sh
##
## composition_unique exists for BrainSEQ only; GTEx rows read an empty DTU-without-DGE set.
## Array = the six trust-funnel regions x {isograph, wgcna}. No model fits.
## Writes 04_module_trust/_m/stability/module_trust/module_complementarity__<cohort>__<region>__<method>.parquet
##
## Bridges memory is --cpus-per-task x 2000MB, so 4 cpus = 8G; do NOT pass --mem.
## Calls the env interpreter directly rather than `module load` + `conda activate`.
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 04_module_trust/_m/logs

PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python
REGIONS=(
    brainseq:caudate brainseq:dlpfc brainseq:hippocampus
    gtex:caudate_basal_ganglia gtex:frontal_cortex_ba9 gtex:hippocampus
)
METHODS=(isograph wgcna)
IDX=$((${SLURM_ARRAY_TASK_ID:-1} - 1))
PAIR="${REGIONS[$((IDX / 2))]}"
METHOD="${METHODS[$((IDX % 2))]}"
COHORT="${PAIR%%:*}"
REGION="${PAIR##*:}"

log "**** Q4 complementarity: ${COHORT}/${REGION}/${METHOD} ****"
"${PY}" -u -m isograph_benchmark.real_data.module_trust complementarity \
    --cohort "${COHORT}" --region "${REGION}" --method "${METHOD}" "$@"
log "**** Complete ****"
