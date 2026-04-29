from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel


DEFAULT_MINUTES = {
    "isograph_baseline": 2.0,
    "isograph_latent": 4.0,
    "isograph_graph": 5.0,
    "isograph_gpu_latent": 4.0,
    "isograph_vae": 12.0,
    "wgcna_gene": 15.0,
}


def estimate_minutes(row: pd.Series) -> float:
    base = DEFAULT_MINUTES[str(row["method"])]
    n_genes = max(float(row["n_genes"]), 1.0)
    n_samples = max(float(row["n_samples"]), 1.0)
    size_factor = (n_genes / 400.0) * math.sqrt(n_samples / 160.0)
    if row["method"] == "wgcna_gene":
        size_factor = (n_genes / 400.0) ** 2 * math.sqrt(n_samples / 160.0)
    if row["method"] == "isograph_vae":
        size_factor = (n_genes / 400.0) * (n_samples / 160.0)
    return max(0.5, base * size_factor)


def load_runtime_estimates(path: Path | None) -> pd.DataFrame | None:
    if path is None or not path.exists():
        return None
    estimates = pd.read_parquet(path)
    required = {"method", "scenario", "n_genes", "estimated_minutes"}
    missing = required - set(estimates.columns)
    if missing:
        raise ValueError(f"Runtime estimate file missing columns: {sorted(missing)}")
    return estimates


def apply_runtime_estimates(grid: pd.DataFrame, estimates: pd.DataFrame | None) -> pd.DataFrame:
    grid = grid.copy()
    grid["estimated_minutes"] = grid.apply(estimate_minutes, axis=1)
    if estimates is None:
        return grid
    keys = ["method", "scenario", "n_genes"]
    merged = grid.merge(estimates[keys + ["estimated_minutes"]], on=keys, how="left", suffixes=("", "_observed"))
    observed = merged["estimated_minutes_observed"].notna()
    merged.loc[observed, "estimated_minutes"] = merged.loc[observed, "estimated_minutes_observed"]
    return merged.drop(columns=["estimated_minutes_observed"])


def make_batches(grid: pd.DataFrame, target_fraction: float = 0.85) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    batch_number = 0
    for resource_class, group in grid.groupby("resource_class", sort=True):
        target = float(group["target_minutes"].iloc[0]) * target_fraction
        current: list[str] = []
        current_minutes = 0.0
        for row in group.sort_values(["scenario", "method", "estimated_minutes"], ascending=[True, True, False]).itertuples():
            est = float(row.estimated_minutes)
            if current and current_minutes + est > target:
                batch_number += 1
                rows.append(
                    {
                        "batch_id": f"B{batch_number:05d}",
                        "resource_class": resource_class,
                        "run_ids": ",".join(current),
                        "n_runs": len(current),
                        "estimated_minutes": current_minutes,
                    }
                )
                current = []
                current_minutes = 0.0
            current.append(str(row.run_id))
            current_minutes += est
        if current:
            batch_number += 1
            rows.append(
                {
                    "batch_id": f"B{batch_number:05d}",
                    "resource_class": resource_class,
                    "run_ids": ",".join(current),
                    "n_runs": len(current),
                    "estimated_minutes": current_minutes,
                }
            )
    return pd.DataFrame(rows)


def write_class_tsvs(batches: pd.DataFrame, out_dir: Path) -> None:
    ensure_dir(out_dir)
    for resource_class, group in batches.groupby("resource_class", sort=True):
        group.reset_index(drop=True).to_csv(out_dir / f"batches_{resource_class}.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", default="benchmark/00_design/_m/synthetic_run_grid.parquet")
    parser.add_argument("--runtime-estimates", default="benchmark/00_design/_m/runtime_estimates.parquet")
    parser.add_argument("--out", default="benchmark/01_synthetic/_m/synthetic_batches.parquet")
    args = parser.parse_args()

    grid_path = rel(args.grid)
    runtime_path = rel(args.runtime_estimates)
    grid = pd.read_parquet(grid_path)
    grid = apply_runtime_estimates(grid, load_runtime_estimates(runtime_path))
    batches = make_batches(grid)
    out = rel(args.out)
    ensure_dir(out.parent)
    batches.to_parquet(out, index=False, compression="zstd")
    write_class_tsvs(batches, out.parent)
    print(f"Wrote {len(batches):,} batches to {out}")
    for resource_class, count in batches["resource_class"].value_counts().sort_index().items():
        print(f"{resource_class}: {count} array tasks")


if __name__ == "__main__":
    main()
