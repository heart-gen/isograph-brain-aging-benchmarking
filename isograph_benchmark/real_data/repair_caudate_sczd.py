"""Repair caudate_sczd artifacts in-place.

Runs Leiden at the target resolution on the already-saved edges.parquet and
rewrites modules / traits / diagnosis_assoc / module_gene_roles without
re-running the expensive VAE training step.

Also fixes the _eigengenes_to_sample_table pivot so that the diagnosis
association file carries M000/M001/... module IDs (not integers), which was
the root cause of the fallback_module_intersection selection failure.

Usage (from project root, with the isograph conda env active):
    python -m isograph_benchmark.real_data.repair_caudate_sczd
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import igraph as ig
import leidenalg
import numpy as np
import pandas as pd

from isograph.io.artifacts import load_dataset_bundle
from isograph.models.base import compute_module_gene_roles, compute_trait_associations
from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.run_models import diagnosis_association

DEFAULT_LEIDEN_RESOLUTION = 10.0
DEFAULT_SEED = 13
DEFAULT_MIN_MODULE_SIZE = 20

ARTIFACT_DIR = rel("real_data", "brainseq", "caudate_sczd", "_m", "isograph_vae")
BUNDLE_PATH = rel("inputs", "bundles", "brainseq_sczd", "caudate")

COVARIATE_COLS = [
    "Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
    "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5",
]
DIAGNOSIS_COVARIATE_COLS = ["Age"] + COVARIATE_COLS


def _build_module_table(
    edges: pd.DataFrame,
    all_gene_ids: list[str],
    leiden_resolution: float,
    seed: int,
    min_module_size: int,
) -> pd.DataFrame:
    """Run Leiden on positive edges; return module_table with M000/M001/... IDs."""
    pos_edges = edges[edges["weight"] > 0].copy()

    # Build node list from all genes in the network (include isolated nodes).
    nodes_list = sorted(all_gene_ids)
    node_to_idx = {n: i for i, n in enumerate(nodes_list)}

    # Only keep edges whose endpoints are in the node list.
    mask = pos_edges["source"].isin(node_to_idx) & pos_edges["target"].isin(node_to_idx)
    pos_edges = pos_edges[mask]

    ig_edges = [
        (node_to_idx[r["source"]], node_to_idx[r["target"]])
        for _, r in pos_edges.iterrows()
    ]
    g = ig.Graph(n=len(nodes_list), edges=ig_edges)

    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=leiden_resolution,
        seed=seed,
    )
    communities = sorted(partition, key=len, reverse=True)

    rows = []
    for module_index, community in enumerate(communities):
        nodes = {nodes_list[v] for v in community}
        if len(nodes) < min_module_size:
            continue
        for gene_id in sorted(nodes):
            rows.append({"gene_id": gene_id, "module_id": f"M{module_index:03d}"})
    return pd.DataFrame(rows)


def _pivot_eigengenes(eigengene_table: pd.DataFrame) -> pd.DataFrame:
    """Correct pivot: rows=samples, cols=[sample_id, M000, M001, ...]."""
    return (
        eigengene_table
        .set_index("module_id")
        .T
        .reset_index()
        .rename(columns={"index": "sample_id"})
        .rename_axis(None, axis=1)
    )


def repair(
    leiden_resolution: float = DEFAULT_LEIDEN_RESOLUTION,
    seed: int = DEFAULT_SEED,
    min_module_size: int = DEFAULT_MIN_MODULE_SIZE,
    dry_run: bool = False,
) -> None:
    t0 = time.time()

    # ── Load saved artifacts ────────────────────────────────────────────────────
    edges_path = ARTIFACT_DIR / "edges.parquet"
    fs_path = ARTIFACT_DIR / "feature_scores.parquet"
    if not edges_path.exists() or not fs_path.exists():
        sys.exit(f"ERROR: missing artifacts in {ARTIFACT_DIR}. Run the SCZD job first.")

    print(f"Loading edges from {edges_path} ...", flush=True)
    edges = pd.read_parquet(edges_path)
    print(f"  {len(edges):,} edges ({(edges['weight']>0).sum():,} positive)")

    print(f"Loading feature_scores from {fs_path} ...", flush=True)
    feature_scores = pd.read_parquet(fs_path)
    all_gene_ids = sorted(feature_scores["gene_id"].astype(str).unique().tolist())
    print(f"  {len(all_gene_ids):,} genes in feature_scores")

    print(f"Loading bundle from {BUNDLE_PATH} ...", flush=True)
    bundle = load_dataset_bundle(BUNDLE_PATH)

    # ── Run Leiden ──────────────────────────────────────────────────────────────
    print(
        f"\nRunning Leiden (resolution={leiden_resolution}, seed={seed}, "
        f"min_module_size={min_module_size}) ...",
        flush=True,
    )
    module_table = _build_module_table(
        edges, all_gene_ids, leiden_resolution, seed, min_module_size
    )
    module_sizes = module_table.groupby("module_id").size().sort_values(ascending=False)
    print(f"  {module_table['module_id'].nunique()} modules | "
          f"giant={module_sizes.iloc[0]} | "
          f"top5={module_sizes.head().tolist()}")

    drd2_gene = "ENSG00000149295"
    drd2_rows = module_table[module_table["gene_id"].str.startswith(drd2_gene)]
    if drd2_rows.empty:
        print("  WARNING: DRD2 not assigned to any module")
    else:
        drd2_mod = drd2_rows["module_id"].iloc[0]
        drd2_size = int(module_sizes[drd2_mod])
        print(f"  DRD2 → {drd2_mod} ({drd2_size} genes)")

    # ── Recompute eigengenes and trait associations ──────────────────────────────
    print("\nRecomputing trait associations and eigengenes ...", flush=True)
    trait_table, eigengene_table = compute_trait_associations(
        module_table, feature_scores, bundle.sample_table, trait_columns=["Dx"]
    )
    print(f"  {len(trait_table)} trait rows")

    # ── Recompute diagnosis associations (with covariates) ──────────────────────
    print("Recomputing diagnosis associations ...", flush=True)
    eg_pivot = _pivot_eigengenes(eigengene_table)
    diagnosis = diagnosis_association(
        eg_pivot,
        bundle.sample_table,
        covariate_cols=DIAGNOSIS_COVARIATE_COLS,
        dx_col="Dx",
        control_label="Control",
        case_label="SCZD",
    )
    if not diagnosis.empty and "fdr" in diagnosis.columns:
        sig = diagnosis[diagnosis["fdr"] <= 0.10]
        print(f"  {len(diagnosis)} module rows | {len(sig)} FDR≤10% modules")
        if not sig.empty:
            print(f"  Significant: {sig['module_id'].tolist()}")

    # ── Recompute module gene roles ─────────────────────────────────────────────
    print("Recomputing module gene roles ...", flush=True)
    module_gene_roles = compute_module_gene_roles(
        module_table, feature_scores, bundle.sample_table
    )
    print(f"  {len(module_gene_roles)} role rows")

    if drd2_rows.empty is False:
        drd2_role = module_gene_roles[module_gene_roles["gene_id"].str.startswith(drd2_gene)]
        if not drd2_role.empty:
            r = drd2_role.iloc[0]
            print(f"  DRD2 role: {r['module_role']} "
                  f"(switch_r={r['switch_r']:.3f}, abundance_r={r['abundance_r']:.3f})")

    if dry_run:
        print("\nDry-run — no files written.")
        return

    # ── Write repaired artifacts ────────────────────────────────────────────────
    print(f"\nWriting repaired artifacts to {ARTIFACT_DIR} ...", flush=True)
    ensure_dir(ARTIFACT_DIR)
    module_table.to_parquet(ARTIFACT_DIR / "modules.parquet", index=False, compression="zstd")
    trait_table.to_parquet(ARTIFACT_DIR / "traits.parquet", index=False, compression="zstd")
    if not diagnosis.empty:
        diagnosis.to_parquet(ARTIFACT_DIR / "diagnosis_assoc.parquet", index=False, compression="zstd")
    if not module_gene_roles.empty:
        module_gene_roles.to_parquet(ARTIFACT_DIR / "module_gene_roles.parquet", index=False, compression="zstd")

    # Invalidate stale interpretation outputs so interpret_modules re-runs.
    interp_dir = ARTIFACT_DIR / "module_interpret"
    if interp_dir.exists():
        stale = list(interp_dir.glob("*.parquet")) + list(interp_dir.glob("*.json"))
        for f in stale:
            f.unlink()
        print(f"  Removed {len(stale)} stale interpretation files from {interp_dir}")

    print(f"\nDone in {time.time() - t0:.0f}s", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Repair caudate SCZD IsoGraph module artifacts.")
    parser.add_argument("--resolution", type=float, default=DEFAULT_LEIDEN_RESOLUTION,
                        help=f"Leiden resolution parameter (default: {DEFAULT_LEIDEN_RESOLUTION})")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help=f"Random seed for Leiden (default: {DEFAULT_SEED})")
    parser.add_argument("--min-module-size", type=int, default=DEFAULT_MIN_MODULE_SIZE,
                        help=f"Minimum module size (default: {DEFAULT_MIN_MODULE_SIZE})")
    parser.add_argument("--dry-run", action="store_true",
                        help="Compute but do not write files")
    args = parser.parse_args()

    repair(
        leiden_resolution=args.resolution,
        seed=args.seed,
        min_module_size=args.min_module_size,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
