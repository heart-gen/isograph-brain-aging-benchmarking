"""List declared synthetic-grid runs that have no ``done.json``, minus the deliberate skips.

The grid declares 14,590 runs. 1,180 of them -- ``isograph_vae_residual`` on the six
unconfounded switch scenarios -- must never be executed: their cached ``samples.parquet``
predates the covariate columns, so ``build_design_matrix`` receives nothing and
residualization silently reduces to plain ``isograph_vae``.  Running them would produce a
"residualization costs nothing" result manufactured by the method not running, and the
sample tables cannot be repaired because 1,204 cached datasets no longer regenerate from the
current generator (``01_synthetic_benchmark/01_synthetic/_m/SAMPLE_TABLE_REFRESH.md``).  That question is
answered instead by ``residual_cost.py`` on fresh, generator-pinned datasets.

Everything else missing is a genuine coverage gap and is emitted here.  The skip list is
encoded rather than left to whoever reads the grid next, because "missing" and "must not be
run" look identical from the outside.

Usage:
    python -m isograph_benchmark.benchmark.list_missing_runs
    python -m isograph_benchmark.benchmark.list_missing_runs --include-skipped   # audit
"""
from __future__ import annotations

import argparse

import pandas as pd

from isograph_benchmark.paths import rel

GRID = ("benchmark", "00_design", "_m", "synthetic_run_grid.parquet")
RUN_ROOT = ("benchmark", "01_synthetic", "_o", "runs")
OUT = ("benchmark", "01_synthetic", "_m", "missing_runs.tsv")

# isograph_vae_residual on these scenarios is a silent no-op; see module docstring.
_NO_OP_RESIDUAL_SCENARIOS = frozenset({
    "idealized_switching", "noise_stress", "feature_space_interactions",
    "unequal_isoform_abundance", "non_switching_background", "negative_control_noise",
})


def missing(include_skipped: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    grid = pd.read_parquet(rel(*GRID))
    root = rel(*RUN_ROOT)
    grid["done"] = [(root / str(r) / "done.json").exists() for r in grid["run_id"]]
    miss = grid.loc[~grid["done"]].copy()
    skip = ((miss["method"] == "isograph_vae_residual")
            & (miss["scenario"].isin(_NO_OP_RESIDUAL_SCENARIOS)))
    return (miss if include_skipped else miss.loc[~skip]), miss.loc[skip]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--include-skipped", action="store_true",
                   help="also emit the deliberate no-op rows (for auditing, not for running)")
    args = p.parse_args()

    run, skipped = missing(args.include_skipped)
    out = rel(*OUT)
    run[["run_id", "scenario", "method", "resource_class", "dataset_id"]].to_csv(
        out, sep="\t", index=False)
    print(f"{len(run)} runs -> {out}")
    if len(run):
        print(run.groupby(["scenario", "method"]).size().to_string())
    if not args.include_skipped and len(skipped):
        print(f"\ndeliberately skipped: {len(skipped)} isograph_vae_residual rows on "
              f"{skipped.scenario.nunique()} unconfounded scenarios (silent no-ops; see "
              f"configs/synthetic_grid.yaml)")


if __name__ == "__main__":
    main()
