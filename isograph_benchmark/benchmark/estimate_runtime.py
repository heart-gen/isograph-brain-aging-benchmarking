from __future__ import annotations

import argparse

import pandas as pd

from isograph_benchmark.benchmark.collect_results import collect
from isograph_benchmark.paths import ensure_dir, rel


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", default="benchmark/01_synthetic/_o/runs")
    parser.add_argument("--out", default="benchmark/00_design/_m/runtime_estimates.parquet")
    parser.add_argument("--multiplier", type=float, default=1.25)
    args = parser.parse_args()

    results = collect(rel(args.run_root))
    if results.empty:
        raise SystemExit("No completed telemetry found. Run a small pilot batch first.")
    completed = results.loc[results["status"] == "completed"].copy()
    completed["estimated_minutes"] = completed["measurement_elapsed_sec"] / 60.0 * args.multiplier
    estimates = (
        completed.groupby(["run_method", "run_scenario", "run_n_genes"], as_index=False)["estimated_minutes"]
        .median()
        .rename(
            columns={
                "run_method": "method",
                "run_scenario": "scenario",
                "run_n_genes": "n_genes",
            }
        )
    )
    out = rel(args.out)
    ensure_dir(out.parent)
    estimates.to_parquet(out, index=False, compression="zstd")
    print(f"Wrote {len(estimates):,} runtime estimates to {out}")


if __name__ == "__main__":
    main()
