#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=smr-heidi
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/smr-heidi-%A_%a.log
#
# SMR + HEIDI on the signal-level colocalization nominations (Analysis 5). One array task per
# row of work_list.tsv -- a (GTEx tissue, chromosome) cell -- which builds the ESD/BESD for
# the target genes there and runs SMR for every analysis with targets in that cell.
#
# This is orthogonal corroboration beneath coloc.susie, not a co-equal analysis. The reading
# rules live in isograph_benchmark/real_data/smr_heidi.py and in SMR_HEIDI.md: b_SMR is not
# causal direction, a non-rejected HEIDI is not proof of a shared variant, and a HEIDI
# rejection does not overrule a strong colocalization.
#
# Memory on PSC is --cpus-per-task x 2000MB; 8 cpus = 16 GB. One task reads one chromosome of
# one tissue's GTEx all-pairs by predicate pushdown for a handful of genes, and SMR holds one
# 2 Mb cis window of 1000G EUR genotypes at a time.
#
# Order (prep and meta are login-node safe):
#   python -m isograph_benchmark.real_data.smr_heidi --stage prep
#   N=$(( $(wc -l < 05_genetic_anchoring/_m/smr_heidi/gtex/work_list.tsv) - 2 ))
#   sbatch --array=0-${N} 05_genetic_anchoring/_h/27.smr_heidi.sh
#   python -m isograph_benchmark.real_data.smr_heidi --stage meta
#
# Depends on: coloc_signal_susie --stage meta --sqtl all (the nominations) and the per-locus
# GWAS from 08.coloc_prep.sh. Run prep again whenever the nominations change: targets, the
# work list and the .ma files are all derived from them.
#
# Instrument-threshold sensitivity arm (writes under sensitivity/peqtl_smr_<x>/; the BESD is
# shared, so run it after a primary array and skip the besd stage):
#   sbatch --array=0-${N} --export=ALL,SMR_PEQTL=1e-5,SMR_SKIP_BESD=1 05_genetic_anchoring/_h/27.smr_heidi.sh
#   python -m isograph_benchmark.real_data.smr_heidi --stage meta --peqtl-smr 1e-5
#
# BrainSEQ QTL source (EA-only; one task per (region, chr); refused unless both BrainSEQ QTL
# checks passed, and read against BrainSEQ's own coloc, so coloc_brainseq --stage meta first):
#   python -m isograph_benchmark.real_data.smr_heidi --stage prep --qtl-source brainseq --arm ea_only
#   N=$(( $(wc -l < 05_genetic_anchoring/_m/smr_heidi/brainseq/ea_only/work_list.tsv) - 2 ))
#   sbatch --array=0-${N} --export=ALL,SMR_QTL_SOURCE=brainseq,SMR_ARM=ea_only 05_genetic_anchoring/_h/27.smr_heidi.sh
#   python -m isograph_benchmark.real_data.smr_heidi --stage meta --qtl-source brainseq --arm ea_only
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
mkdir -p 05_genetic_anchoring/_m/logs

TASK="${SLURM_ARRAY_TASK_ID:-${1:-0}}"
PEQTL="${SMR_PEQTL:-5e-8}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
SRC="${SMR_QTL_SOURCE:-gtex}"
SRC_ARGS=(--qtl-source "${SRC}")
[[ -n "${SMR_ARM:-}" ]] && SRC_ARGS+=(--arm "${SMR_ARM}")

# Guard: some Bridges2 batch nodes start array tasks without Lmod initialised.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ocean/projects/bio260021p/shared/opt/envs/isograph

if [[ -z "${SMR_SKIP_BESD:-}" ]]; then
    log "**** SMR BESD: ${SRC} task ${TASK} ****"
    python -m isograph_benchmark.real_data.smr_heidi --stage besd --task "${TASK}" "${SRC_ARGS[@]}"
fi
log "**** SMR: ${SRC} task ${TASK}, --peqtl-smr ${PEQTL} ****"
python -m isograph_benchmark.real_data.smr_heidi --stage smr --task "${TASK}" \
    --peqtl-smr "${PEQTL}" --threads "${THREADS}" "${SRC_ARGS[@]}"
conda deactivate
log "**** Complete: task ${TASK} ****"
