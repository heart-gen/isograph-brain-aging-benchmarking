"""Sweep Leiden resolution over existing saved IsoGraph artifacts.

Reuses saved edges.parquet + feature_scores.parquet without refitting the VAE.
Reports module statistics at each resolution and optionally writes the best result.

Usage (from project root, with the isograph conda env active):
    # Sweep SCZD caudate:
    python -m isograph_benchmark.real_data.sweep_leiden --analysis brainseq-sczd

    # Sweep aging regions (all or specific):
    python -m isograph_benchmark.real_data.sweep_leiden --analysis brainseq-aging
    python -m isograph_benchmark.real_data.sweep_leiden --analysis brainseq-aging --region caudate

    # Dry-run (compute, do not write):
    python -m isograph_benchmark.real_data.sweep_leiden --analysis brainseq-sczd --dry-run

    # Sweep and commit the best resolution to the artifact directory:
    python -m isograph_benchmark.real_data.sweep_leiden --analysis brainseq-sczd --write-best
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import igraph as ig
import leidenalg
import numpy as np
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score

from isograph.io.artifacts import load_dataset_bundle
from isograph.models.base import compute_module_gene_roles, compute_trait_associations
from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.go_enrichment import GoAnnotations, HAS_GOATOOLS
from isograph_benchmark.real_data.run_models import (
    diagnosis_association,
    linear_age_association,
    spline_age_association,
)

DRD2_GENE = "ENSG00000149295"

# Resolutions to sweep per analysis type
SCZD_RESOLUTIONS = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0]
# Caudate jumps from 11→33 modules between 2.0 and 3.0; include fine grid there.
AGING_RESOLUTIONS = [0.5, 1.0, 1.5, 2.0, 2.25, 2.5, 2.75, 3.0, 5.0]

# Hard constraint: reject resolutions where the giant module swamps the partition.
GIANT_FRACTION_MAX = 0.30

# Default GO annotation cache (relative to project root)
_DEFAULT_GO_CACHE = rel("inputs", "go_annotations")

SCZD_COVARIATE_COLS = [
    "Age", "Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
    "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5",
]
AGING_COVARIATE_COLS = [
    "Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
    "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5",
]
# GTEx v11 has no genotype PCs in the bundle; QC covariates mirror run_gtex_region.
# AGE is the continuous trait (years), handled by the same df=3 spline path as aging.
GTEX_COVARIATE_COLS = ["SEX", "SMRIN", "SMTSISCH", "SMMAPRT"]
GTEX_AGE_COL = "AGE"


def _compute_nmi(mt_a: pd.DataFrame, mt_b: pd.DataFrame) -> float:
    """NMI between two Leiden partitions over their common assigned genes.

    Returns the arithmetic-mean NMI ∈ [0, 1].  A value close to 1 means the
    two partitions are nearly identical (stable); a sudden drop between adjacent
    resolutions flags a fragmentation event.
    """
    common = set(mt_a["gene_id"]) & set(mt_b["gene_id"])
    if len(common) < 10:
        return np.nan
    genes = sorted(common)
    a = mt_a.set_index("gene_id")["module_id"].loc[genes].values
    b = mt_b.set_index("gene_id")["module_id"].loc[genes].values
    return float(normalized_mutual_info_score(a, b, average_method="arithmetic"))


def _select_best_resolution(results: pd.DataFrame) -> float:
    """Data-driven resolution selection.

    Primary criterion: go_n_enriched — the absolute number of modules with ≥1
    significant GO Biological Process term (BH FDR ≤ 0.05).  Maximising the
    count of biologically coherent modules avoids the resolution bias of the
    density (fraction) metric: density systematically favours low resolutions
    because larger modules have more genes and therefore more enrichment power
    per module.  Absolute counts instead reward the resolution that produces the
    most biologically interpretable modules overall.

    Hard constraint: giant_fraction ≤ GIANT_FRACTION_MAX (degenerate partitions
    excluded regardless of GO score).

    NMI (nmi_to_prev) is reported in the sweep table as a stability diagnostic
    but does not drive selection — a drop in NMI alone does not indicate a worse
    partition, only a different one.

    Fallback when GO data are unavailable: maximise n_sig_spline_fdr10 (aging)
    or n_sig_linear_fdr10 (SCZD), both of which are trait-driven and still
    independent of WGCNA.
    """
    df = results[results["giant_fraction"] <= GIANT_FRACTION_MAX].copy()
    if df.empty:
        df = results.copy()

    if "go_n_enriched" in df.columns and df["go_n_enriched"].gt(0).any():
        # Break ties by n_modules descending: prefer the most resolved partition
        # when GO evidence is equal (avoids collapsing to the coarsest tied option).
        best_idx = (
            df.sort_values(["go_n_enriched", "n_modules"], ascending=[False, False])
            .index[0]
        )
    elif "n_sig_spline_fdr10" in df.columns:
        best_idx = df["n_sig_spline_fdr10"].idxmax()
    else:
        best_idx = df["n_sig_linear_fdr10"].idxmax()

    return float(df.loc[best_idx, "leiden_resolution"])


def _artifact_dir(analysis: str, region: str | None, variant: str = "standard") -> Path:
    # variant "with-abundance" sweeps the abundance-channel refit (separate dir),
    # so both IsoGraph variants can be resolution-selected by the same GO criterion.
    subdir = "isograph_vae_with_abundance" if variant == "with-abundance" else "isograph_vae"
    if analysis == "brainseq-sczd":
        return rel("real_data", "brainseq", "caudate_sczd", "_m", subdir)
    if analysis == "brainseq-aging":
        assert region is not None
        return rel("real_data", "brainseq", region, "_m", subdir)
    if analysis == "gtex-aging":
        assert region is not None
        return rel("real_data", "gtex", region, "_m", subdir)
    raise ValueError(f"Unknown analysis: {analysis!r}")


def _bundle_path(analysis: str, region: str | None) -> Path:
    if analysis == "brainseq-sczd":
        return rel("inputs", "bundles", "brainseq_sczd", "caudate")
    if analysis == "brainseq-aging":
        assert region is not None
        return rel("inputs", "bundles", "brainseq_v1", region)
    if analysis == "gtex-aging":
        assert region is not None
        return rel("inputs", "bundles", "gtex_v11_brain", region)
    raise ValueError(f"Unknown analysis: {analysis!r}")


def _build_module_table(
    edges: pd.DataFrame,
    all_gene_ids: list[str],
    leiden_resolution: float,
    seed: int = 13,
    min_module_size: int = 20,
) -> pd.DataFrame:
    """Run Leiden on positive edges; return module_table with M000/M001/... IDs."""
    pos_edges = edges[edges["weight"] > 0].copy()
    nodes_list = sorted(all_gene_ids)
    node_to_idx = {n: i for i, n in enumerate(nodes_list)}
    mask = pos_edges["source"].isin(node_to_idx) & pos_edges["target"].isin(node_to_idx)
    pos_edges = pos_edges[mask]
    ig_edges = [
        (node_to_idx[r["source"]], node_to_idx[r["target"]])
        for _, r in pos_edges.iterrows()
    ]
    g = ig.Graph(n=len(nodes_list), edges=ig_edges)
    partition = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=leiden_resolution, seed=seed,
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
    return (
        eigengene_table
        .set_index("module_id").T.reset_index()
        .rename(columns={"index": "sample_id"})
        .rename_axis(None, axis=1)
    )


def _module_stats(module_table: pd.DataFrame) -> dict:
    if module_table.empty:
        return {"n_modules": 0, "n_assigned": 0, "giant_size": 0, "giant_fraction": 0.0, "top5_sizes": []}
    sizes = module_table.groupby("module_id").size().sort_values(ascending=False)
    n_assigned = len(module_table)
    giant = int(sizes.iloc[0])
    return {
        "n_modules": int(len(sizes)),
        "n_assigned": n_assigned,
        "giant_size": giant,
        "giant_fraction": round(giant / n_assigned, 4) if n_assigned > 0 else 0.0,
        "top5_sizes": sizes.head(5).tolist(),
    }


def _drd2_info(module_table: pd.DataFrame) -> tuple[str | None, int | None]:
    hits = module_table[module_table["gene_id"].str.startswith(DRD2_GENE)]
    if hits.empty:
        return None, None
    mod = hits["module_id"].iloc[0]
    size = int(module_table.groupby("module_id").size()[mod])
    return mod, size


def sweep_one(
    analysis: str,
    region: str | None,
    resolutions: list[float],
    seed: int = 13,
    min_module_size: int = 20,
    dry_run: bool = False,
    write_best: bool = False,
    go_cache_dir: Path | None = None,
    skip_go: bool = False,
    variant: str = "standard",
) -> pd.DataFrame:
    artifact_dir = _artifact_dir(analysis, region, variant)
    label = f"{analysis}/{region}" if region else analysis
    if variant != "standard":
        label = f"{label}[{variant}]"

    edges_path = artifact_dir / "edges.parquet"
    fs_path = artifact_dir / "feature_scores.parquet"
    if not edges_path.exists() or not fs_path.exists():
        print(f"[{label}] ERROR: missing artifacts in {artifact_dir}. Run the model first.")
        return pd.DataFrame()

    print(f"[{label}] Loading edges ({edges_path}) ...", flush=True)
    t0 = time.time()
    edges = pd.read_parquet(edges_path)
    print(f"  {len(edges):,} edges ({(edges['weight'] > 0).sum():,} positive) in {time.time()-t0:.0f}s")

    print(f"[{label}] Loading feature_scores ...", flush=True)
    feature_scores = pd.read_parquet(fs_path)
    all_gene_ids = sorted(feature_scores["gene_id"].astype(str).unique().tolist())
    print(f"  {len(all_gene_ids):,} genes")

    print(f"[{label}] Loading bundle ...", flush=True)
    bundle = load_dataset_bundle(_bundle_path(analysis, region))
    sample_table = bundle.sample_table

    # Initialise GO enrichment helper (shared across all resolutions in this sweep).
    go_helper: GoAnnotations | None = None
    if not skip_go and HAS_GOATOOLS:
        try:
            cache = go_cache_dir or _DEFAULT_GO_CACHE
            go_helper = GoAnnotations(cache)
            go_helper.prepare(all_gene_ids)
        except Exception as exc:
            print(f"[{label}] WARNING: GO enrichment disabled ({exc})", flush=True)
            go_helper = None
    elif not skip_go and not HAS_GOATOOLS:
        print(f"[{label}] WARNING: goatools not installed; skipping GO enrichment.", flush=True)

    rows = []
    prev_module_table: pd.DataFrame | None = None

    for res in resolutions:
        t1 = time.time()
        print(f"\n[{label}] Leiden resolution={res} ...", flush=True)
        module_table = _build_module_table(edges, all_gene_ids, res, seed, min_module_size)
        stats = _module_stats(module_table)
        # DRD2 is striatal — only meaningful in the SCZD caudate dataset
        drd2_mod, drd2_size = _drd2_info(module_table) if analysis == "brainseq-sczd" else (None, None)

        # Partition stability: NMI between this resolution and the previous one.
        nmi_to_prev = (
            _compute_nmi(prev_module_table, module_table)
            if prev_module_table is not None
            else np.nan
        )
        prev_module_table = module_table

        # GO enrichment density
        go_density = np.nan
        go_n_enriched = 0
        if go_helper is not None and not module_table.empty:
            try:
                go_density, go_n_enriched = go_helper.enrichment_density(module_table)
            except Exception as exc:
                print(f"  WARNING: GO enrichment failed: {exc}")

        # Trait associations
        n_sig_linear = 0
        n_sig_spline = 0
        if not module_table.empty:
            try:
                trait_col = "Dx" if analysis == "brainseq-sczd" else "Age"
                _, eigengene_table = compute_trait_associations(
                    module_table, feature_scores, sample_table,
                    trait_columns=[trait_col],
                )
                eg_pivot = _pivot_eigengenes(eigengene_table)
                if analysis == "brainseq-sczd":
                    assoc = diagnosis_association(
                        eg_pivot, sample_table,
                        covariate_cols=SCZD_COVARIATE_COLS,
                    )
                    if not assoc.empty and "fdr" in assoc.columns:
                        n_sig_linear = int((assoc["fdr"] <= 0.10).sum())
                else:
                    linear = linear_age_association(eg_pivot, sample_table, age_col="Age")
                    spline = spline_age_association(
                        eg_pivot, sample_table,
                        covariate_cols=AGING_COVARIATE_COLS, age_col="Age",
                    )
                    if not linear.empty and "fdr" in linear.columns:
                        n_sig_linear = int((linear["fdr"] <= 0.10).sum())
                    if not spline.empty and "fdr_ftest" in spline.columns:
                        n_sig_spline = int(
                            spline[spline["fdr_ftest"] <= 0.10]["module_id"].nunique()
                        )
            except Exception as exc:
                print(f"  WARNING: association failed: {exc}")

        row = {
            "analysis": analysis,
            "region": region or "caudate_sczd",
            "leiden_resolution": res,
            **stats,
            "drd2_module": drd2_mod,
            "drd2_size": drd2_size,
            "nmi_to_prev": round(nmi_to_prev, 4) if np.isfinite(nmi_to_prev) else None,
            "go_enrichment_density": round(go_density, 4) if np.isfinite(go_density) else None,
            "go_n_enriched": go_n_enriched,
            "n_sig_linear_fdr10": n_sig_linear,
            "n_sig_spline_fdr10": n_sig_spline,
            "elapsed_s": round(time.time() - t1, 1),
        }
        rows.append(row)
        go_str = (
            f"go_density={go_density:.2f}({go_n_enriched}/{stats['n_modules']}) | "
            if np.isfinite(go_density) else ""
        )
        print(
            f"  n_modules={stats['n_modules']} | giant={stats['giant_size']} "
            f"({stats['giant_fraction']:.1%}) | "
            + (f"DRD2={drd2_mod}({drd2_size}) | " if drd2_mod else "")
            + (f"nmi={nmi_to_prev:.3f} | " if np.isfinite(nmi_to_prev) else "")
            + go_str
            + f"n_sig_linear={n_sig_linear} n_sig_spline={n_sig_spline} | {row['elapsed_s']:.0f}s"
        )

    results = pd.DataFrame(rows)
    print(f"\n[{label}] Sweep summary:")
    summary_cols = [
        "leiden_resolution", "n_modules", "giant_fraction",
        "nmi_to_prev", "go_enrichment_density", "go_n_enriched",
        "n_sig_linear_fdr10", "n_sig_spline_fdr10",
    ]
    print(results[[c for c in summary_cols if c in results.columns]].to_string(index=False))

    if not dry_run:
        out_path = artifact_dir / "leiden_sweep_results.parquet"
        results.to_parquet(out_path, index=False, compression="zstd")
        print(f"[{label}] Sweep results written to {out_path}")

    if write_best and not dry_run:
        best_resolution = _select_best_resolution(results)
        _r = results.loc[results["leiden_resolution"] == best_resolution].iloc[0]
        print(
            f"[{label}] Selected resolution={best_resolution} "
            f"(GO n_enriched={_r['go_n_enriched']}, density={_r['go_enrichment_density']})",
            flush=True,
        )
        best_module_table = _build_module_table(
            edges, all_gene_ids, best_resolution, seed, min_module_size
        )
        print(f"[{label}] Writing best artifacts ...", flush=True)
        _write_best_artifacts(
            analysis, region, best_module_table, feature_scores,
            bundle, artifact_dir, best_resolution, seed, dry_run=False,
        )

    return results


def _write_best_artifacts(
    analysis: str,
    region: str | None,
    module_table: pd.DataFrame,
    feature_scores: pd.DataFrame,
    bundle,
    artifact_dir: Path,
    resolution: float,
    seed: int,
    dry_run: bool,
) -> None:
    label = f"{analysis}/{region}" if region else analysis
    sample_table = bundle.sample_table

    trait_col = "Dx" if analysis == "brainseq-sczd" else "Age"
    trait_table, eigengene_table = compute_trait_associations(
        module_table, feature_scores, sample_table, trait_columns=[trait_col],
    )
    eg_pivot = _pivot_eigengenes(eigengene_table)

    if analysis == "brainseq-sczd":
        assoc = diagnosis_association(eg_pivot, sample_table, covariate_cols=SCZD_COVARIATE_COLS)
        n_sig = int((assoc["fdr"] <= 0.10).sum()) if not assoc.empty and "fdr" in assoc.columns else 0
        # DRD2 is striatal — only report for SCZD caudate
        drd2_mod, drd2_size = _drd2_info(module_table)
        drd2_str = f" | DRD2={drd2_mod}({drd2_size})" if drd2_mod else ""
    else:
        linear = linear_age_association(eg_pivot, sample_table, age_col="Age")
        spline = spline_age_association(
            eg_pivot, sample_table, covariate_cols=AGING_COVARIATE_COLS, age_col="Age",
        )
        n_sig_linear = int((linear["fdr"] <= 0.10).sum()) if not linear.empty and "fdr" in linear.columns else 0
        n_sig_spline = int(
            spline[spline["fdr_ftest"] <= 0.10]["module_id"].nunique()
        ) if not spline.empty and "fdr_ftest" in spline.columns else 0
        n_sig = n_sig_linear  # use linear for summary line; spline printed separately
        drd2_str = ""

    module_gene_roles = compute_module_gene_roles(module_table, feature_scores, sample_table)

    sizes = module_table.groupby("module_id").size().sort_values(ascending=False)
    print(f"  {module_table['module_id'].nunique()} modules | "
          f"giant={sizes.iloc[0]} | n_sig_fdr10={n_sig}{drd2_str}")
    if analysis != "brainseq-sczd":
        print(f"  n_sig_linear={n_sig_linear} | n_sig_spline={n_sig_spline}")

    if dry_run:
        return

    module_table.to_parquet(artifact_dir / "modules.parquet", index=False, compression="zstd")
    trait_table.to_parquet(artifact_dir / "traits.parquet", index=False, compression="zstd")
    if analysis == "brainseq-sczd":
        if not assoc.empty:
            assoc.to_parquet(artifact_dir / "diagnosis_assoc.parquet", index=False, compression="zstd")
    else:
        if not linear.empty:
            linear.to_parquet(artifact_dir / "age_linear.parquet", index=False, compression="zstd")
        if not spline.empty:
            spline.to_parquet(artifact_dir / "age_spline.parquet", index=False, compression="zstd")
    if not module_gene_roles.empty:
        module_gene_roles.to_parquet(artifact_dir / "module_gene_roles.parquet", index=False, compression="zstd")

    # Invalidate stale interpretation outputs
    interp_dir = artifact_dir / "module_interpret"
    if interp_dir.exists():
        stale = list(interp_dir.glob("*.parquet")) + list(interp_dir.glob("*.json"))
        for f in stale:
            f.unlink()
        print(f"  Removed {len(stale)} stale interpretation files")

    print(f"[{label}] Best artifacts written (resolution={resolution})")


def main() -> None:
    parser = argparse.ArgumentParser(description="Sweep Leiden resolution on saved IsoGraph artifacts.")
    parser.add_argument(
        "analysis", nargs="?", default="brainseq-sczd",
        choices=["brainseq-sczd", "brainseq-aging"],
        help="Which analysis to sweep (default: brainseq-sczd)",
    )
    parser.add_argument(
        "--region", action="append", dest="regions",
        help="For brainseq-aging: which region(s) to sweep. Repeatable. Default: all 3.",
    )
    parser.add_argument(
        "--resolutions", type=float, nargs="+",
        help="Override the default resolution grid (space-separated floats).",
    )
    parser.add_argument(
        "--seed", type=int, default=13,
        help="Leiden random seed (default: 13)",
    )
    parser.add_argument(
        "--min-module-size", type=int, default=20,
        help="Minimum module size (default: 20)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Compute but do not write any files.",
    )
    parser.add_argument(
        "--write-best", action="store_true",
        help="After sweeping, write the best-resolution artifacts (selected by GO density).",
    )
    parser.add_argument(
        "--go-cache-dir", type=Path, default=None,
        help=(
            "Directory for cached GO annotation files (go-basic.obo, gene2go.gz, "
            "Homo_sapiens.gene_info.gz). Default: inputs/go_annotations/"
        ),
    )
    parser.add_argument(
        "--no-go", action="store_true",
        help="Skip GO enrichment (falls back to trait-association criterion for --write-best).",
    )
    parser.add_argument(
        "--variant", choices=["standard", "with-abundance"], default="standard",
        help=(
            "Which saved artifact set to sweep: 'standard' (isograph_vae/) or "
            "'with-abundance' (isograph_vae_with_abundance/). Default: standard."
        ),
    )
    args = parser.parse_args()

    go_kws = dict(go_cache_dir=args.go_cache_dir, skip_go=args.no_go)

    if args.analysis == "brainseq-sczd":
        resolutions = args.resolutions or SCZD_RESOLUTIONS
        sweep_one(
            "brainseq-sczd", None, resolutions,
            seed=args.seed, min_module_size=args.min_module_size,
            dry_run=args.dry_run, write_best=args.write_best,
            variant=args.variant,
            **go_kws,
        )
    elif args.analysis == "brainseq-aging":
        resolutions = args.resolutions or AGING_RESOLUTIONS
        regions = args.regions or ["caudate", "hippocampus", "dlpfc"]
        for region in regions:
            sweep_one(
                "brainseq-aging", region, resolutions,
                seed=args.seed, min_module_size=args.min_module_size,
                dry_run=args.dry_run, write_best=args.write_best,
                variant=args.variant,
                **go_kws,
            )


if __name__ == "__main__":
    main()
