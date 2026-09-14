"""Is the ``isograph_vae_gpu`` arm reproducible, and from what?

The PI review called the GPU arm "preserved but not reproducible from config". The config half
is wrong -- ``run_one.build_model`` has an ``isograph_vae_gpu`` case identical to
``isograph_vae`` except ``device="cuda"`` -- but three separate questions sit behind the label,
and a full re-run answers none of them on its own:

1. **Run-to-run determinism on GPU.** Does the same run, same code, same seed, give the same
   result twice? If not, no re-run can make the arm reproducible; only deterministic kernels can.
2. **Code drift.** The stored rows were produced by IsoGraph 0.1.2 (1,204) and 0.1.4 (126);
   current is newer. Does current code reproduce the stored values -- for the GPU arm *and* for
   its CPU twin? If the CPU twin drifts too, the issue is grid-wide, not a GPU property.
3. **Device equivalence.** Documentation called the two backends "numerically matched"; the
   stored pairs are bit-identical in only 334 of 1,330 module-recovery values.

This probe re-runs a scenario-stratified sample into isolated output roots -- GPU twice, CPU once
-- and compares everything against the stored telemetry. It never touches the stored runs and
never regenerates a dataset (``select`` refuses any run whose cached dataset is missing).

  python -m isograph_benchmark.benchmark.gpu_reproducibility_probe select
  sbatch 01_synthetic_benchmark/01_synthetic/_h/run_gpu_repro_probe.sh      # array over rows
  python -m isograph_benchmark.benchmark.gpu_reproducibility_probe compare
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel

GRID = ("01_synthetic_benchmark", "00_design", "_m", "synthetic_run_grid.parquet")
DATASETS = ("01_synthetic_benchmark", "01_synthetic", "_m", "datasets")
STORED_RUNS = ("01_synthetic_benchmark", "01_synthetic", "_o", "runs")
PROBE_RUNS = ("01_synthetic_benchmark", "01_synthetic", "_o", "gpu_repro_probe")
OUT = ("01_synthetic_benchmark", "03_metrics", "_m", "gpu_repro_probe")
ARMS = ("gpu_rep1", "gpu_rep2", "cpu_rep1")
METRICS = ("module_recovery", "ari_planted", "n_predicted_modules", "n_edges",
           "switch_gene_detection_rate")
_KEY = ["dataset_id", "seed", "replicate"]


def select(per_scenario: int, seed: int) -> Path:
    grid = pd.read_parquet(rel(*GRID))
    gpu = grid[grid["method"] == "isograph_vae_gpu"]
    cpu = grid[grid["method"] == "isograph_vae"][_KEY + ["run_id"]].rename(
        columns={"run_id": "run_id_cpu"})
    pairs = gpu.merge(cpu, on=_KEY, how="inner").rename(columns={"run_id": "run_id_gpu"})
    ds_root = rel(*DATASETS)
    pairs = pairs[[(ds_root / d / "manifest.json").exists() for d in pairs["dataset_id"]]]
    rng = np.random.default_rng(seed)
    picks = []
    for _, g in pairs.groupby("scenario", sort=True):
        take = g.iloc[rng.permutation(len(g))[:per_scenario]]
        picks.append(take)
    out = pd.concat(picks, ignore_index=True)[
        ["scenario", "dataset_id", "seed", "replicate", "run_id_gpu", "run_id_cpu"]]
    out.insert(0, "probe_index", np.arange(1, len(out) + 1))
    dest = ensure_dir(rel(*OUT)) / "probe_runs.tsv"
    out.to_csv(dest, sep="\t", index=False)
    print(f"{len(out)} probe rows over {out['scenario'].nunique()} scenarios -> {dest}")
    return dest


def _metrics(path: Path) -> dict:
    if not path.exists():
        return {}
    t = json.loads(path.read_text())
    m = t.get("metrics") or {}
    return {k: m.get(k) for k in METRICS} | {
        "status": t.get("status"),
        "isograph_version": (t.get("software") or {}).get("isograph"),
    }


def compare() -> pd.DataFrame:
    probe = pd.read_csv(rel(*OUT) / "probe_runs.tsv", sep="\t")
    stored, probe_root = rel(*STORED_RUNS), rel(*PROBE_RUNS)
    rows = []
    for r in probe.itertuples(index=False):
        rec = {"probe_index": r.probe_index, "scenario": r.scenario}
        sources = {
            "stored_gpu": stored / r.run_id_gpu / "telemetry.json",
            "stored_cpu": stored / r.run_id_cpu / "telemetry.json",
            "gpu_rep1": probe_root / "gpu_rep1" / r.run_id_gpu / "telemetry.json",
            "gpu_rep2": probe_root / "gpu_rep2" / r.run_id_gpu / "telemetry.json",
            "cpu_rep1": probe_root / "cpu_rep1" / r.run_id_cpu / "telemetry.json",
        }
        for name, path in sources.items():
            for k, v in _metrics(path).items():
                rec[f"{name}__{k}"] = v
        rows.append(rec)
    df = pd.DataFrame(rows)
    ensure_dir(rel(*OUT))
    df.to_parquet(rel(*OUT) / "probe_comparison.parquet", index=False)

    contrasts = {
        "GPU run-to-run (rep1 vs rep2, current code)": ("gpu_rep1", "gpu_rep2"),
        "GPU current code vs stored GPU": ("gpu_rep1", "stored_gpu"),
        "CPU current code vs stored CPU": ("cpu_rep1", "stored_cpu"),
        "GPU vs CPU, current code": ("gpu_rep1", "cpu_rep1"),
        "GPU vs CPU, stored": ("stored_gpu", "stored_cpu"),
    }
    summary = []
    for label, (a, b) in contrasts.items():
        for m in METRICS:
            x = pd.to_numeric(df.get(f"{a}__{m}"), errors="coerce")
            y = pd.to_numeric(df.get(f"{b}__{m}"), errors="coerce")
            ok = x.notna() & y.notna()
            d = (x[ok] - y[ok]).abs()
            summary.append({
                "contrast": label, "metric": m, "n": int(ok.sum()),
                "n_identical": int((d == 0).sum()),
                "max_abs_diff": float(d.max()) if len(d) else np.nan,
                "median_abs_diff": float(d.median()) if len(d) else np.nan,
            })
    s = pd.DataFrame(summary)
    s.to_parquet(rel(*OUT) / "probe_summary.parquet", index=False)
    _report(df, s)
    print(s.to_string(index=False))
    return s


def _report(df: pd.DataFrame, s: pd.DataFrame) -> None:
    cols = list(s.columns)
    table = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in s.itertuples(index=False):
        table.append("| " + " | ".join(
            f"{v:.4g}" if isinstance(v, float) else str(v) for v in r) + " |")
    versions = {c: sorted(df[c].dropna().astype(str).unique().tolist())
                for c in df.columns if c.endswith("__isograph_version")}
    lines = [
        "# GPU VAE reproducibility probe",
        "",
        f"`gpu_reproducibility_probe.py` over {len(df)} runs stratified by scenario. Each run "
        "was re-executed from the current grid config into an isolated root — GPU twice, its "
        "CPU twin once — and compared with the stored telemetry. Stored runs were not touched "
        "and no dataset was regenerated.",
        "",
        f"IsoGraph versions seen: `{json.dumps(versions)}`.",
        "",
        "## How to read it",
        "",
        "- **GPU run-to-run** separates nondeterminism from everything else: if rep1 and rep2 "
        "differ, re-running cannot make the arm reproducible.",
        "- **Current vs stored, per device** separates code drift from device effects: if the "
        "CPU twin also fails to reproduce its stored value, the drift is grid-wide.",
        "- **GPU vs CPU** tests the documented claim that the backends are numerically matched.",
        "",
        *table,
        "",
    ]
    (rel(*OUT) / "GPU_REPRO_PROBE.md").write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--per-scenario", type=int, default=3)
    s.add_argument("--seed", type=int, default=13)
    sub.add_parser("compare")
    args = ap.parse_args()
    if args.cmd == "select":
        select(args.per_scenario, args.seed)
    else:
        compare()


if __name__ == "__main__":
    main()
