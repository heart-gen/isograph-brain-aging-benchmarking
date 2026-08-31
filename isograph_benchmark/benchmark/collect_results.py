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


def add_compute_fields(row: dict[str, Any]) -> None:
    requested_gpus = row.get("run_requested_gpus")
    slurm_gpus = row.get("slurm_slurm_gpus") or row.get("slurm_slurm_job_gpus")
    cuda_available = row.get("torch_cuda_available")
    gpu_query = row.get("hardware_nvidia_smi_gpu_query")

    has_gpu = False
    if requested_gpus not in (None, "", 0, "0"):
        has_gpu = True
    if slurm_gpus not in (None, "", "0"):
        has_gpu = True
    if cuda_available is True or str(cuda_available).lower() == "true":
        has_gpu = True
    if gpu_query not in (None, ""):
        has_gpu = True

    row["compute_backend"] = "GPU" if has_gpu else "CPU"
    row["accelerator_name"] = None
    if gpu_query:
        row["accelerator_name"] = str(gpu_query).split(",", 1)[0].strip()
    row["cpu_model"] = row.get("hardware_cpu_model")
    row["requested_cpus"] = row.get("run_requested_cpus")
    row["requested_gpus"] = row.get("run_requested_gpus")
    row["requested_mem_gb"] = row.get("run_requested_mem_gb")
    row["max_rss_mb"] = row.get("measurement_max_rss_mb")
    row["gpu_peak_allocated_mb"] = row.get("torch_gpu_peak_allocated_mb")
    row["gpu_peak_reserved_mb"] = row.get("torch_gpu_peak_reserved_mb")


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
        add_compute_fields(row)
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", default="01_synthetic_benchmark/01_synthetic/_o/runs")
    parser.add_argument("--out", default="01_synthetic_benchmark/01_synthetic/_m/synthetic_results.parquet")
    args = parser.parse_args()
    results = collect(rel(args.run_root))
    out = rel(args.out)
    ensure_dir(out.parent)
    results.to_parquet(out, index=False, compression="zstd")
    print(f"Wrote {len(results):,} result rows to {out}")


if __name__ == "__main__":
    main()
