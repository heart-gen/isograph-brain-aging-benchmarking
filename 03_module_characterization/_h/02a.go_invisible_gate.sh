#!/usr/bin/env bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=go-invisible-gate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=03_module_characterization/_m/logs/go-invisible-gate-%j.log

## Calls the env interpreter directly rather than `module load` + `conda activate`, so the
## job cannot die on "module: command not found".
set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -f .here || ! -d isograph_benchmark ]]; then
    echo "ERROR: submit from the repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export PYTHONPATH="${PROJECT_ROOT}:/ocean/projects/bio260021p/kbenjamin/software/IsoGraph/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p 03_module_characterization/_m/logs

GTF=/ocean/projects/bio250020p/shared/resources/genomes/human/gencode-v47/gtf/gencode.v47.primary_assembly.annotation.gtf
SYM=03_module_characterization/_m/tmp/gene_id_symbol.parquet
PY=/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python

# Build the gene_id -> symbol cache once (awk over gene lines; the in-script python
# fallback reads the whole GTF and is slow on a login node).
if [[ ! -f "${SYM}" ]]; then
    log_message "building gene-symbol cache ..."
    mkdir -p 03_module_characterization/_m/tmp
    awk -F'\t' '$3=="gene"{
        match($9,/gene_id "[^."]+/); gid=substr($9,RSTART+9,RLENGTH-9);
        match($9,/gene_name "[^"]+/); nm=substr($9,RSTART+11,RLENGTH-11);
        print gid"\t"nm}' "${GTF}" > 03_module_characterization/_m/tmp/gene_id_symbol.tsv
    "${PY}" - <<'PY'
import pandas as pd
df = pd.read_csv("03_module_characterization/_m/tmp/gene_id_symbol.tsv", sep="\t",
                 names=["gene_id", "gene_name"]).drop_duplicates("gene_id")
df.to_parquet("03_module_characterization/_m/tmp/gene_id_symbol.parquet", index=False, compression="zstd")
PY
fi

# Default: the SCZD disease gate. Override the analysis/region via "$@", e.g.
#   sbatch 03_module_characterization/_h/02a.go_invisible_gate.sh --analysis brainseq-aging --region caudate
ARGS=("$@")
if [[ ${#ARGS[@]} -eq 0 ]]; then ARGS=(--analysis brainseq-sczd); fi

log_message "**** GO-invisible biology gate: ${ARGS[*]} ****"
"${PY}" -m isograph_benchmark.real_data.go_invisible_gate "${ARGS[@]}"
log_message "**** Complete ****"
