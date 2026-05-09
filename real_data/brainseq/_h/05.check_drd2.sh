#!/usr/bin/env bash
# Check whether DRD2 (ENSG00000149295) is present in the SCZD caudate IsoGraph modules.
# Exits 0 (pass) or 1 (fail). Run locally or on HPC after 02.run_isograph_sczd.sh completes.
set -euo pipefail

PROJECT_ROOT="${ISOGRAPH_BENCHMARK_ROOT:-${SLURM_SUBMIT_DIR:-${PWD}}}"
cd "${PROJECT_ROOT}"
if [[ ! -d isograph_benchmark ]]; then
    echo "ERROR: run from the isograph-brain-aging-benchmarking repo root or set ISOGRAPH_BENCHMARK_ROOT."
    exit 1
fi
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

python3 - <<'EOF'
import sys
import pandas as pd
from pathlib import Path

modules_path = Path("real_data/brainseq/caudate_sczd/_m/isograph_vae/modules.parquet")
tx_path = Path("inputs/bundles/brainseq_sczd/caudate/transcripts.parquet")

if not modules_path.exists():
    print(f"ERROR: modules file not found: {modules_path}")
    sys.exit(1)
if not tx_path.exists():
    print(f"ERROR: transcript table not found: {tx_path}")
    sys.exit(1)

tx = pd.read_parquet(tx_path)
drd2_mask = tx.get("transcript_name", pd.Series(dtype=str)).str.startswith("DRD2-", na=False)
drd2_gene_ids = tx.loc[drd2_mask, "gene_id"].unique().tolist()

if not drd2_gene_ids:
    print("FAIL: DRD2 gene_id not found in transcript table")
    sys.exit(1)

modules = pd.read_parquet(modules_path)
hits = modules[modules["gene_id"].isin(drd2_gene_ids)]

if hits.empty:
    print(f"FAIL: DRD2 ({drd2_gene_ids}) not found in any module")
    sys.exit(1)

module_ids = hits["module_id"].unique().tolist()
print(f"PASS: DRD2 ({drd2_gene_ids[0]}) found in module(s): {module_ids}")

roles_path = Path("real_data/brainseq/caudate_sczd/_m/isograph_vae/module_gene_roles.parquet")
if roles_path.exists():
    roles = pd.read_parquet(roles_path)
    role_hits = roles[roles["gene_id"].isin(drd2_gene_ids)]
    if not role_hits.empty:
        for _, r in role_hits.iterrows():
            print(f"  gene_id={r['gene_id']}  module_id={r.get('module_id','?')}  role={r.get('module_role','?')}")
EOF
