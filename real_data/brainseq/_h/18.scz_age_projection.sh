#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=scz-age-projection
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=real_data/brainseq/_m/logs/scz-age-projection-%j.log

## Do SCZ-risk loci converge on age-sensitive isoform-switch programs disrupted in disease?
## Projects the age-sensitive SCZ-GWAS-enriched aging caudate modules (defined out-of-cohort in
## the independent GTEx caudate-basal-ganglia + BrainSeq caudate fits) onto the independent SCZ
## case/control cohort (BrainSeq caudate_sczd), testing convergence + module-level directional
## disruption (headline) and, per convergent module, its candidate RBP trans-regulators.
##
##   stage 1 (isograph env): prep genotype extract/keep lists for the SCZ coloc lead loci.
##   stage 2 (eqtl env):     plink2 --export A dosages for those loci in the caudate_sczd samples.
##   stage 3 (isograph env): scz_age_projection.py (A/A3 supp, B/C/D, convergence + mechanism).
##
## Depends on: coloc events (coloc_isoform_events_combined.parquet), MAGMA (magma_results_combined),
## the aging + disease IsoGraph fits, and the per-module RBP regulons (rbp_regulon --scope combined,
## from 16/17.rbp_regulon*.sh) for the mechanism layer. TOPMed LIBD genotypes are controlled-access;
## the extracted dosage table is a run-local intermediate and is NOT committed.
## Usage: sbatch real_data/brainseq/_h/18.scz_age_projection.sh
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
[[ -f .here && -d isograph_benchmark ]] || { echo "ERROR: submit from repo root."; exit 1; }
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p real_data/brainseq/_m/logs

ISO_ENV=/ocean/projects/bio260021p/shared/opt/envs/isograph
EQTL_ENV=/ocean/projects/bio260021p/shared/opt/envs/eqtl
PLINK2="${EQTL_ENV}/bin/plink2"
GENO=/ocean/projects/bio250020p/shared/resources/libd_data/genotypes/AA_EA/_m/TOPMed_LIBD
OUTDIR=real_data/_m/scz_age_projection/genotypes
mkdir -p "${OUTDIR}"

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module purge
module load anaconda3/2024.10-1
source "$(conda info --base)/etc/profile.d/conda.sh"

log "**** stage 0: regenerate per-module RBP regulons, all scopes (isograph env) ****"
## stage-2 only (motif counts from 16/17.rbp_regulon*.sh); picks up the module-gene-pool
## fallback so caudate_basal_ganglia + other GTEx regions are covered for the mechanism layer.
conda activate "${ISO_ENV}"
for scope in mature intronic combined; do
    python -m isograph_benchmark.real_data.rbp_regulon --scope "${scope}"
done

log "**** stage 1: genotype extract/keep lists (isograph env) ****"
python - "${GENO}" "${OUTDIR}" <<'PY'
import sys, pandas as pd
from isograph.io.artifacts import load_dataset_bundle
from isograph_benchmark.paths import rel
geno, outdir = sys.argv[1], sys.argv[2]
coloc = pd.read_parquet(rel("real_data", "coloc", "_m", "coloc_isoform_events_combined.parquet"))
rsids = set(coloc[coloc.trait.astype(str).str.lower().eq("scz")].best_rsid.dropna().astype(str))
ids = []
with open(f"{geno}.pvar") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        vid = f[2]; toks = vid.split("_")
        if toks and toks[-1].startswith("rs") and toks[-1] in rsids:
            ids.append(vid)
open(f"{outdir}/extract_ids.txt", "w").write("\n".join(ids) + "\n")
st = load_dataset_bundle(rel("inputs", "bundles", "brainseq_sczd", "caudate")).sample_table
brs = set(st["BrNum"].astype(str))
psam = pd.read_csv(f"{geno}.psam", sep="\t")
psam.columns = [c.lstrip("#") for c in psam.columns]
psam[psam["FID"].astype(str).isin(brs)][["FID", "IID"]].to_csv(
    f"{outdir}/keep.txt", sep="\t", index=False, header=False)
print(f"extract {len(ids)} variants; keep {psam['FID'].astype(str).isin(brs).sum()} samples")
PY
conda deactivate

log "**** stage 2: plink2 dosage export (eqtl env) ****"
conda activate "${EQTL_ENV}"
"${PLINK2}" --pfile "${GENO}" \
    --extract "${OUTDIR}/extract_ids.txt" --keep "${OUTDIR}/keep.txt" \
    --export A --out "${OUTDIR}/scz_loci_dosage"
conda deactivate

log "**** stage 3: SCZ age-projection (isograph env) ****"
conda activate "${ISO_ENV}"
python -m isograph_benchmark.real_data.scz_age_projection \
    --dosage "${OUTDIR}/scz_loci_dosage.raw" --n-perm 2000
conda deactivate
log "**** SCZ age-projection done ****"
