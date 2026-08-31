"""Backfill the partition metrics (reviewer items 1 and 4) onto completed benchmark runs.

The 13k runs under ``01_synthetic_benchmark/01_synthetic/_o/runs`` each keep their fitted
``modules.parquet``, and every dataset keeps its ``truth_modules.parquet`` — so ARI, AMI,
the homogeneity/completeness decomposition and the module-count-preserving Jaccard null can
all be recomputed without re-fitting a single model.

New keys are merged into each run's ``telemetry.json`` / ``done.json`` ``metrics`` block, so
``collect_results.py`` picks them up unchanged.  Existing keys are never touched, and every
run's published ``module_recovery`` is re-derived and checked against the stored value, so a
silent convention drift between the null and the metric it calibrates cannot pass.

Usage (SLURM array):
    python -m isograph_benchmark.benchmark.backfill_metrics --shard "$SLURM_ARRAY_TASK_ID" \
        --n-shards 16
"""
from __future__ import annotations

import argparse
import json
import os
import zlib
from pathlib import Path
from typing import Any

import pandas as pd

from isograph_benchmark.benchmark.partition_metrics import (
    DEFAULT_N_PERM,
    DEFAULT_SEED,
    METRIC_KEYS,
    partition_metrics,
)
from isograph_benchmark.paths import ensure_dir, rel

# Recomputed vs stored best-match Jaccard must agree to this tolerance on every run.
MRS_TOLERANCE = 1e-9


def _run_dirs(run_root: Path) -> list[Path]:
    # Ocean FS: iterdir(), not glob() (AGENTS.md).
    return sorted(p for p in run_root.iterdir() if p.is_dir())


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    """Write via a temp file + os.replace so an interrupted job cannot truncate a run dir."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, default=str))
    os.replace(tmp, path)


class _TruthCache:
    """Truth tables keyed by dataset dir — 2,290 datasets serve 13,060 runs."""

    def __init__(self, max_entries: int = 64) -> None:
        self._cache: dict[str, tuple[pd.DataFrame, list[str]]] = {}
        self._max = max_entries

    def get(self, dataset_path: Path) -> tuple[pd.DataFrame, list[str]] | None:
        key = str(dataset_path)
        if key in self._cache:
            return self._cache[key]
        modules = dataset_path / "truth_modules.parquet"
        switch = dataset_path / "truth_switch.parquet"
        if not modules.exists():
            return None
        truth = pd.read_parquet(modules)
        universe: list[str] = []
        if switch.exists():
            universe = pd.read_parquet(switch, columns=["gene_id"])["gene_id"].astype(str).tolist()
        if len(self._cache) >= self._max:
            self._cache.pop(next(iter(self._cache)))
        self._cache[key] = (truth, universe)
        return self._cache[key]


def _seed_for(run_id: str, seed: int) -> int:
    """Run-specific but deterministic seed, so shards and reruns agree.

    ``hash()`` is salted per process (PYTHONHASHSEED), so it cannot be used here — CRC32 is
    stable across processes and machines.
    """
    return (seed * 1_000_003 + zlib.crc32(run_id.encode())) % (2**31 - 1)


def backfill_run(
    run_dir: Path,
    truth_cache: _TruthCache,
    n_perm: int,
    seed: int,
    force: bool,
) -> dict[str, Any]:
    """Compute and persist the new metrics for one run; returns a QC row."""
    run_id = run_dir.name
    row: dict[str, Any] = {"run_id": run_id, "partition_status": "ok", "mrs_delta": None}

    telemetry_path = run_dir / "telemetry.json"
    telemetry = _read_json(telemetry_path)
    if telemetry is None:
        row["partition_status"] = "missing_telemetry"
        return row

    metrics = telemetry.get("metrics") or {}
    row["dataset_id"] = (telemetry.get("run") or {}).get("dataset_id")
    row["scenario"] = (telemetry.get("run") or {}).get("scenario")
    row["method"] = (telemetry.get("run") or {}).get("method")
    if not force and "ari_planted" in metrics:
        row["partition_status"] = "skipped_present"
        return row

    dataset_path = telemetry.get("dataset_path")
    if not dataset_path:
        row["partition_status"] = "missing_dataset_path"
        return row
    truth_entry = truth_cache.get(Path(dataset_path))
    if truth_entry is None:
        row["partition_status"] = "missing_truth"
        return row
    truth, universe = truth_entry

    modules_path = run_dir / "modules.parquet"
    if not modules_path.exists():
        row["partition_status"] = "missing_modules"
        predicted = pd.DataFrame(columns=["gene_id", "module_id"])
    else:
        predicted = pd.read_parquet(modules_path)
        if predicted.empty or "module_id" not in predicted.columns:
            row["partition_status"] = "empty_modules"
            predicted = pd.DataFrame(columns=["gene_id", "module_id"])

    new_metrics = partition_metrics(
        predicted,
        truth,
        universe=universe or None,
        seed=_seed_for(run_id, seed),
        n_perm=n_perm,
    )

    # Guard: the null is only a calibration of the published metric if it is the same
    # computation.  Compare against the stored value whenever both are defined.
    stored = metrics.get("module_recovery")
    recomputed = new_metrics.get("module_recovery_recomputed")
    if stored is not None and recomputed is not None:
        row["mrs_delta"] = float(recomputed) - float(stored)
        if abs(row["mrs_delta"]) > MRS_TOLERANCE:
            row["partition_status"] = "mrs_mismatch"

    metrics.update(new_metrics)
    telemetry["metrics"] = metrics
    _write_json_atomic(telemetry_path, telemetry)

    done_path = run_dir / "done.json"
    done = _read_json(done_path)
    if done is not None:
        done_metrics = done.get("metrics") or {}
        done_metrics.update(new_metrics)
        done["metrics"] = done_metrics
        _write_json_atomic(done_path, done)

    row.update({k: new_metrics.get(k) for k in METRIC_KEYS})
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-root", default="01_synthetic_benchmark/01_synthetic/_o/runs")
    ap.add_argument("--out-dir", default="01_synthetic_benchmark/01_synthetic/_m/partition_backfill")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--n-perm", type=int, default=DEFAULT_N_PERM)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--limit", type=int, default=None, help="debug: process at most N runs")
    ap.add_argument("--force", action="store_true", help="recompute runs that already carry the metrics")
    ap.add_argument("--dry-run", action="store_true", help="compute and report, but write nothing")
    args = ap.parse_args()

    run_root = rel(args.run_root)
    runs = _run_dirs(run_root)
    runs = [r for i, r in enumerate(runs) if i % args.n_shards == args.shard]
    if args.limit:
        runs = runs[: args.limit]
    print(f"shard {args.shard}/{args.n_shards}: {len(runs):,} runs, n_perm={args.n_perm}", flush=True)

    truth_cache = _TruthCache()
    rows = []
    for i, run_dir in enumerate(runs, 1):
        if args.dry_run:
            telemetry = _read_json(run_dir / "telemetry.json") or {}
            entry = truth_cache.get(Path(telemetry.get("dataset_path", "")))
            rows.append({"run_id": run_dir.name, "partition_status": "dry_run", "resolvable": entry is not None})
        else:
            rows.append(backfill_run(run_dir, truth_cache, args.n_perm, args.seed, args.force))
        if i % 250 == 0:
            print(f"  {i:,}/{len(runs):,}", flush=True)

    qc = pd.DataFrame(rows)
    out_dir = ensure_dir(rel(args.out_dir))
    out = out_dir / f"partition_backfill_shard{args.shard:03d}.parquet"
    qc.to_parquet(out, index=False, compression="zstd")

    counts = qc["partition_status"].value_counts().to_dict()
    print(f"status: {counts}", flush=True)
    if "mrs_mismatch" in counts:
        raise SystemExit(
            f"ERROR: {counts['mrs_mismatch']} runs disagree with the stored module_recovery "
            f"beyond {MRS_TOLERANCE:g} — the null does not calibrate the published metric."
        )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
