#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=GPU-shared
#SBATCH --job-name=brainseq-swqtl
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=05_genetic_anchoring/_m/logs/brainseq-swqtl-%j.log
#
# BrainSEQ cis switch-QTL (S_g) and matched abundance-QTL (A_g) mapping.
#
# This is SAME-TISSUE GENETIC ANCHORING, not independent replication: BrainSEQ is the
# cohort the switches were discovered in. GTEx remains the external cohort.
#
# Two stages, two environments, because they need different toolchains:
#   phenotypes  (isograph env) -- recompute S_g / A_g on the expanded sample set with the
#               same `gene_feature_channels` the discovery fits use, and write tensorQTL
#               BED phenotypes + covariates keyed by BrNum.
#   map         (eqtl env)     -- tensorQTL cis nominal + permutation for BOTH modalities
#               under identical covariates, cis window and MAF floor.
#
# GPU-shared, because the permutation pass is the binding cost and tensorQTL is written
# for the GPU. Measured on a login-node CPU: chr22 alone (402 switch phenotypes) had not
# finished 10,000 permutations after 7 minutes. The full run is ~43x chr22 per modality
# and there are two modalities and three regions, so the CPU path is several days and the
# GPU path is a few hours. `--gres=gpu:1` with 5 cpus follows the convention already used
# by 01_synthetic_benchmark/01_synthetic/_h/run_batch_gpu.sh -- 5 is not a stylistic
# choice, it is the ceiling: GPU-shared refuses "cpus-per-gpu higher than maximum of
# 5/gpu" at submission.
#
# Permutations are seeded (brainseq_switch_qtl._SEED = 13). Note that seeding fixes the
# permutation draw, not bitwise cross-device equality: CPU and GPU reductions differ in
# floating-point association, so permutation p-values can differ in the last digits
# between devices. Run a whole arm on one device.
#
# Genotypes are controlled-access. The phenotype BEDs are per-donor expression matrices
# and the nominal cis output runs to GBs; both are gitignored and regenerated from here.
#
# Usage:
#   sbatch 05_genetic_anchoring/_h/01i.brainseq_switch_qtl.sh                  # all regions, both stages
#   sbatch 05_genetic_anchoring/_h/01i.brainseq_switch_qtl.sh caudate          # one region
#   sbatch --export=ALL,SWQTL_ARM=ea_only 05_genetic_anchoring/_h/01i.brainseq_switch_qtl.sh    # the coloc arm
#     (ea_only adjusts for genotype PCs computed WITHIN the EA panel -- `ARM_PCS` in
#     brainseq_switch_qtl.py; all_samples keeps the bundle's multi-ancestry PCs)
#   sbatch --export=ALL,SWQTL_STAGE=map   05_genetic_anchoring/_h/01i.brainseq_switch_qtl.sh    # re-map only
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export ISOGRAPH_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 05_genetic_anchoring/_m/logs

# `all_samples` is primary for the modality contrast; `ea_only` is required for anything
# colocalized against the EUR GWAS (LD in a ~50% AA cohort matches neither the EUR GWAS
# nor the 1000G EUR panel).
ARM="${SWQTL_ARM:-all_samples}"
STAGE="${SWQTL_STAGE:-both}"
REGIONS=("${@:-}")
if [[ -z "${REGIONS[0]:-}" ]]; then REGIONS=(caudate dlpfc hippocampus); fi
REGION_ARGS=(); for r in "${REGIONS[@]}"; do REGION_ARGS+=(--region "${r}"); done

# Guard: some Bridges2 batch nodes start without Lmod initialised.
if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"

ISO_ENV=/ocean/projects/bio260021p/shared/opt/envs/isograph
EQTL_ENV=/ocean/projects/bio260021p/shared/opt/envs/eqtl

if [[ "${STAGE}" == "both" || "${STAGE}" == "phenotypes" ]]; then
    log "**** phenotypes: arm=${ARM} regions=${REGIONS[*]} ****"
    conda activate "${ISO_ENV}"
    python -m isograph_benchmark.real_data.brainseq_switch_qtl \
        --stage phenotypes --arm "${ARM}" "${REGION_ARGS[@]}"
    conda deactivate
fi

if [[ "${STAGE}" == "both" || "${STAGE}" == "map" ]]; then
    log "**** map: arm=${ARM} regions=${REGIONS[*]} ****"
    conda activate "${EQTL_ENV}"
    python -m isograph_benchmark.real_data.brainseq_switch_qtl \
        --stage map --arm "${ARM}" "${REGION_ARGS[@]}"
    conda deactivate
fi

log "**** Complete: arm=${ARM} ****"
