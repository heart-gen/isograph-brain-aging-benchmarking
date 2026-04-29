from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel


def batch_row(batch_file: Path, batch_index: int | None, batch_id: str | None) -> pd.Series:
    batches = pd.read_csv(batch_file, sep="\t") if batch_file.suffix == ".tsv" else pd.read_parquet(batch_file)
    if batch_id is not None:
        matches = batches.loc[batches["batch_id"].astype(str) == batch_id]
    else:
        if batch_index is None:
            raise ValueError("Provide --batch-index or --batch-id")
        matches = batches.iloc[[batch_index - 1]]
    if len(matches) != 1:
        raise ValueError(f"Expected one batch, found {len(matches)}")
    return matches.iloc[0]


def run_ids(row: pd.Series) -> list[str]:
    return [value for value in str(row["run_ids"]).split(",") if value]


def run_one(run_id: str, grid: str, log_dir: Path, force: bool = False) -> int:
    ensure_dir(log_dir)
    cmd = [
        sys.executable,
        "-m",
        "isograph_benchmark.benchmark.run_one",
        "--grid",
        grid,
        "--run-id",
        run_id,
    ]
    if force:
        cmd.append("--force")
    stdout = log_dir / f"{run_id}.stdout"
    stderr = log_dir / f"{run_id}.stderr"
    time_file = log_dir / f"{run_id}.time.txt"
    if shutil.which("/usr/bin/time"):
        cmd = ["/usr/bin/time", "-v", "-o", str(time_file), *cmd]
    with stdout.open("w") as out, stderr.open("w") as err:
        proc = subprocess.run(cmd, stdout=out, stderr=err, text=True, check=False)
    return int(proc.returncode)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-file", default="benchmark/01_synthetic/_m/synthetic_batches.parquet")
    parser.add_argument("--batch-index", type=int)
    parser.add_argument("--batch-id")
    parser.add_argument("--grid", default="benchmark/00_design/_m/synthetic_run_grid.parquet")
    parser.add_argument("--log-dir", default="benchmark/01_synthetic/_m/logs")
    parser.add_argument("--stop-margin-min", type=float, default=10.0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    batch_file = rel(args.batch_file)
    row = batch_row(batch_file, args.batch_index, args.batch_id)
    ids = run_ids(row)
    log_dir = rel(args.log_dir) / str(row["batch_id"])
    start = time.monotonic()
    failures = 0
    for idx, run_id in enumerate(ids, start=1):
        elapsed_min = (time.monotonic() - start) / 60.0
        if elapsed_min + args.stop_margin_min >= float(row["estimated_minutes"]):
            print(f"Stopping batch {row['batch_id']} before run {idx}; safety margin reached.")
            break
        code = run_one(run_id, args.grid, log_dir, force=args.force)
        if code != 0:
            failures += 1
            print(f"Run {run_id} failed with exit code {code}")
    if failures:
        raise SystemExit(f"{failures} runs failed in batch {row['batch_id']}")


if __name__ == "__main__":
    main()
