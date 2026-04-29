from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel


def flatten(prefix: str, value: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, item in value.items():
        name = f"{prefix}_{key}" if prefix else key
        if isinstance(item, dict):
            out.update(flatten(name, item))
        else:
            out[name] = item
    return out


def collect(run_root: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for telemetry_path in sorted(run_root.glob("*/telemetry.json")):
        data = json.loads(telemetry_path.read_text())
        row: dict[str, Any] = {"telemetry_path": str(telemetry_path)}
        for section in ["run", "metrics", "measurement", "hardware", "slurm", "software", "torch", "wgcna"]:
            value = data.get(section, {})
            if isinstance(value, dict):
                row.update(flatten(section, value))
        row["status"] = data.get("status")
        error = data.get("error")
        if isinstance(error, dict):
            row["error_type"] = error.get("type")
            row["error_message"] = error.get("message")
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", default="benchmark/01_synthetic/_o/runs")
    parser.add_argument("--out", default="benchmark/01_synthetic/_m/synthetic_results.parquet")
    args = parser.parse_args()
    results = collect(rel(args.run_root))
    out = rel(args.out)
    ensure_dir(out.parent)
    results.to_parquet(out, index=False, compression="zstd")
    print(f"Wrote {len(results):,} result rows to {out}")


if __name__ == "__main__":
    main()
