"""Rewrite cached synthetic ``samples.parquet`` files -- and nothing else.

Why this exists
---------------
``isograph_vae_residual`` asks for the covariates ``RIN``, ``neuron_frac``, ``batch`` and
``library_size`` (``run_one.py:116``).  1,730 of the cached datasets were generated before
``synthetic_data.py`` learned to emit those columns and carry only
``sample_id, Dx, Age, Sex``, so ``build_design_matrix`` receives nothing and residualization
is a **silent no-op** -- the residual arm reduces to plain ``isograph_vae``.

That matters because the open benchmark gap is whether residualization is *free* on
unconfounded data.  Running that ablation against a stale sample table would return "costs
exactly nothing" as an artifact of the method never running.  A clean-looking confirmation
of the hypothesis, manufactured by a bug.

``ensure_dataset`` cannot fix this: it short-circuits on ``manifest.json`` and so never
revisits a materialised dataset.

What it does
------------
The sample table is not hand-writable -- ``library_size`` is an RNG draw interleaved with
the expression construction (``synthetic_data.py:231``) -- so the only correct source is the
generator itself.  For each dataset this rebuilds the full bundle from its grid row, then:

1. fingerprints every cached artifact **except** ``samples.parquet`` (count matrices,
   feature tables, truth tables),
2. fingerprints the same artifacts on the rebuilt bundle,
3. writes ``samples.parquet`` **only if every one of those fingerprints matches**.

The regenerated matrices are then discarded.  13,060 completed runs and every published
number key off the cached expression data, so overwriting it would mean that a single
non-reproducing dataset silently invalidates results already in the figures.  Discarding
makes that failure mode structurally impossible: the worst case is a loud mismatch.

A mismatch is therefore a finding, not an error to route around -- it says the cached grid
is no longer reproducible from the committed config, which is a bigger problem than a stale
sample table and must be surfaced rather than patched over.

Usage
-----
    # probe: rebuild ONE stale dataset and report whether it reproduces (seconds, no writes)
    python -m isograph_benchmark.benchmark.refresh_sample_tables --probe

    # full pass (SLURM); aborts before any write if the built-in probe fails
    python -m isograph_benchmark.benchmark.refresh_sample_tables --shard 0 --n-shards 8

Outputs (under ``01_synthetic_benchmark/01_synthetic/_m/``):
    sample_table_refresh__shard{i}.parquet   one row per dataset, with both fingerprints
    SAMPLE_TABLE_REFRESH.md                  written by --collect once all shards land
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from isograph.io.artifacts import load_dense_matrix
from isograph.validation import DatasetManifest
from isograph_benchmark.benchmark.synthetic_data import build_synthetic_bundle, dataset_dir
from isograph_benchmark.paths import ensure_dir, rel

GRID_PATH = ("benchmark", "00_design", "_m", "synthetic_run_grid.parquet")
DATASET_ROOT = ("benchmark", "01_synthetic", "_m", "datasets")
OUT_ROOT = ("benchmark", "01_synthetic", "_m")

# The covariates isograph_vae_residual regresses out; a table missing any of them makes
# residualization a no-op.  median_tin is deliberately excluded -- synthetic_data only emits
# it under active 3' degradation, so its absence is correct for every other scenario.
REQUIRED_COVARIATES = ("RIN", "neuron_frac", "batch", "library_size")


# --------------------------------------------------------------------------- #
# Fingerprints
# --------------------------------------------------------------------------- #
def _hash_array(arr: np.ndarray) -> str:
    arr = np.ascontiguousarray(arr)
    h = hashlib.sha256()
    h.update(str(arr.dtype).encode())
    h.update(str(arr.shape).encode())
    h.update(arr.tobytes())
    return h.hexdigest()


def _hash_frame(df: pd.DataFrame) -> str:
    """Value-level digest, independent of parquet encoding or row index."""
    h = hashlib.sha256()
    for col in df.columns:
        h.update(str(col).encode())
        h.update(str(df[col].dtype).encode())
        h.update(pd.util.hash_pandas_object(df[col], index=False).to_numpy().tobytes())
    return h.hexdigest()


def _cached_fingerprint(ddir: Path, manifest: DatasetManifest) -> dict[str, str]:
    """Everything in the dataset dir EXCEPT the sample table."""
    parts: dict[str, str] = {}
    for spec in manifest.matrices:
        parts[f"m:{spec.assay_name}"] = _hash_array(load_dense_matrix(ddir / spec.filename))
    for spec in manifest.feature_tables:
        parts[f"f:{spec.kind}"] = _hash_frame(pd.read_parquet(ddir / spec.filename))
    for name in manifest.truth_tables:
        parts[f"t:{name}"] = _hash_frame(pd.read_parquet(ddir / name))
    return parts


def _built_fingerprint(bundle) -> dict[str, str]:
    parts: dict[str, str] = {}
    for name, arr in bundle.matrices.items():
        parts[f"m:{name}"] = _hash_array(arr)
    for kind, tbl in bundle.feature_tables.items():
        parts[f"f:{kind}"] = _hash_frame(tbl)
    for name, tbl in bundle.truth_tables.items():
        parts[f"t:{name}"] = _hash_frame(tbl)
    return parts


def _digest(parts: dict[str, str]) -> str:
    """Order-independent roll-up of a fingerprint dict, for compact reporting."""
    h = hashlib.sha256()
    for key in sorted(parts):
        h.update(key.encode())
        h.update(parts[key].encode())
    return h.hexdigest()[:16]


# --------------------------------------------------------------------------- #
# Per-dataset work
# --------------------------------------------------------------------------- #
def _unique_datasets() -> pd.DataFrame:
    """One row per dataset_id.  The grid repeats a dataset once per method; the dataset
    parameters are identical across those rows, so the first is representative."""
    grid = pd.read_parquet(rel(*GRID_PATH))
    return grid.drop_duplicates(subset=["dataset_id"], keep="first").reset_index(drop=True)


def _is_stale(sample_table: pd.DataFrame) -> bool:
    return any(c not in sample_table.columns for c in REQUIRED_COVARIATES)


def process(row: pd.Series, root: Path, write: bool) -> dict:
    ddir = dataset_dir(row, root)
    rec = {
        "dataset_id": str(row["dataset_id"]),
        "scenario": str(row["scenario"]),
        "seed": int(row["seed"]),
        "status": "",
        "stale_before": None,
        "cached_digest": None,
        "built_digest": None,
        "n_mismatched_artifacts": 0,
        "mismatched": "",
        "n_new_artifacts": 0,
        "new_artifacts": "",
    }
    if not (ddir / "manifest.json").exists():
        rec["status"] = "absent"
        return rec

    manifest = DatasetManifest.model_validate_json((ddir / "manifest.json").read_text())
    cached_samples = pd.read_parquet(ddir / manifest.sample_table)
    rec["stale_before"] = _is_stale(cached_samples)

    cached = _cached_fingerprint(ddir, manifest)
    bundle = build_synthetic_bundle(row)
    built = _built_fingerprint(bundle)

    rec["cached_digest"] = _digest(cached)
    rec["built_digest"] = _digest(built)

    # Three ways the two fingerprints can disagree, and only two of them are faults.
    #   shared-but-different -> the cached data is NOT reproducible.  Hard stop.
    #   cached-only          -> the generator no longer emits something the runs consumed.
    #                           Also a hard stop: the cache cannot be regenerated.
    #   built-only           -> the generator gained an output after this dataset was
    #                           materialised (e.g. truth_switch_event), exactly like the
    #                           sample-table columns this tool exists to add.  Benign:
    #                           nothing cached changed, so record it and carry on.
    bad = sorted({k for k in set(cached) & set(built) if cached[k] != built[k]}
                 | (set(cached) - set(built)))
    new = sorted(set(built) - set(cached))
    rec["n_new_artifacts"] = len(new)
    rec["new_artifacts"] = ",".join(new[:8])
    if bad:
        # Do NOT write.  The cached expression data is what every published number was
        # computed from; a mismatch means it is not reproducible from the committed config.
        rec["status"] = "mismatch"
        rec["n_mismatched_artifacts"] = len(bad)
        rec["mismatched"] = ",".join(bad[:8])
        return rec

    if _hash_frame(cached_samples) == _hash_frame(bundle.sample_table):
        rec["status"] = "unchanged"
        return rec

    if write:
        _atomic_parquet(bundle.sample_table, ddir / manifest.sample_table)
        rec["status"] = "refreshed"
    else:
        rec["status"] = "would_refresh"
    return rec


def _atomic_parquet(df: pd.DataFrame, dest: Path) -> None:
    tmp = dest.with_suffix(dest.suffix + f".tmp{os.getpid()}")
    df.to_parquet(tmp, index=False)
    os.replace(tmp, dest)


# --------------------------------------------------------------------------- #
# Probe
# --------------------------------------------------------------------------- #
def probe(root: Path, datasets: pd.DataFrame) -> dict | None:
    """Rebuild the first stale dataset and report whether it reproduces.  Returns the
    record, or None when no stale dataset is materialised (nothing to do)."""
    for row in datasets.itertuples(index=False):
        row = pd.Series(row._asdict())
        ddir = dataset_dir(row, root)
        if not (ddir / "manifest.json").exists():
            continue
        manifest = DatasetManifest.model_validate_json((ddir / "manifest.json").read_text())
        if not _is_stale(pd.read_parquet(ddir / manifest.sample_table)):
            continue
        return process(row, root, write=False)
    return None


def _report_probe(rec: dict | None) -> bool:
    if rec is None:
        print("probe: no stale dataset found on disk -- nothing to refresh.")
        return True
    print(f"probe dataset {rec['dataset_id']} ({rec['scenario']}, seed {rec['seed']})")
    print(f"  cached digest {rec['cached_digest']}   built digest {rec['built_digest']}")
    if rec["status"] == "mismatch":
        print(f"  MISMATCH in {rec['n_mismatched_artifacts']} shared artifact(s): "
              f"{rec['mismatched']}")
        print("  The cached expression data does not reproduce from the committed config.")
        print("  Refusing to refresh. This is a reproducibility finding, not a transient "
              "error -- do not work around it.")
        return False
    print(f"  status: {rec['status']} -- every shared artifact reproduces exactly.")
    if rec["n_new_artifacts"]:
        print(f"  note: the generator now also emits {rec['n_new_artifacts']} artifact(s) "
              f"absent from this cache ({rec['new_artifacts']}). Nothing cached changed; "
              f"this tool writes only the sample table and leaves those alone.")
    return True


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def run(shard: int, n_shards: int, root: Path, write: bool, limit: int | None) -> int:
    datasets = _unique_datasets()

    # The probe gates every write in this process: if the generator does not reproduce the
    # cached data, nothing downstream is trustworthy and no file should be touched.
    if not _report_probe(probe(root, datasets)):
        return 2
    print()

    mine = datasets.iloc[shard::n_shards] if n_shards > 1 else datasets
    if limit:
        mine = mine.head(limit)

    records = []
    for i, row in enumerate(mine.itertuples(index=False), start=1):
        records.append(process(pd.Series(row._asdict()), root, write))
        if i % 100 == 0:
            print(f"  {i}/{len(mine)} datasets", flush=True)

    out = pd.DataFrame(records)
    outdir = ensure_dir(rel(*OUT_ROOT))
    out.to_parquet(outdir / f"sample_table_refresh__shard{shard}.parquet", index=False)

    counts = out["status"].value_counts().to_dict()
    print(f"\nshard {shard}/{n_shards}: {len(out)} datasets -> {counts}")

    n_bad = int((out["status"] == "mismatch").sum())
    if n_bad:
        print(f"ERROR: {n_bad} dataset(s) did not reproduce; their sample tables were left "
              f"untouched. See sample_table_refresh__shard{shard}.parquet.")
        return 1
    return 0


def collect() -> int:
    outdir = ensure_dir(rel(*OUT_ROOT))
    frames = [pd.read_parquet(p) for p in sorted(outdir.iterdir())
              if p.name.startswith("sample_table_refresh__shard") and p.suffix == ".parquet"]
    if not frames:
        print("no shard outputs found")
        return 1
    out = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["dataset_id"])
    out.to_parquet(outdir / "sample_table_refresh.parquet", index=False)

    counts = out["status"].value_counts()
    bad = out[out["status"] == "mismatch"]
    lines = [
        "# Synthetic sample-table refresh", "",
        "> **Archival decision (2026-08-14).** All benchmark datasets will be archived and "
        "deposited alongside the manuscript. Reproducibility of the synthetic benchmark is "
        "therefore provided **by the archive**, not by regeneration from the generator: the "
        "`mismatch` datasets below are a valid realization of the documented generative "
        "process that current code no longer reproduces, and they are kept as-is rather "
        "than being regenerated. Anything generated after this decision carries a generator "
        "pin (cf. `residual_cost_generator_pin.json`) so the same drift is detectable "
        "rather than silent.", "",
        f"{len(out)} unique datasets in the run grid. Only `samples.parquet` is rewritten; "
        "the regenerated count matrices, feature tables and truth tables are fingerprinted "
        "against the cached copies and then discarded, so a dataset that fails to reproduce "
        "is reported rather than overwritten.", "",
        "| status | n | meaning |", "|---|---|---|",
    ]
    meaning = {
        "refreshed": "sample table rewritten; all other artifacts reproduced exactly",
        "unchanged": "sample table already matched the generator",
        "would_refresh": "dry run; would have rewritten the sample table",
        "mismatch": "**a shared artifact did not reproduce; left untouched**",
        "absent": "dataset never materialised on disk",
    }
    for status, n in counts.items():
        lines.append(f"| `{status}` | {n} | {meaning.get(status, '')} |")

    stale = out["stale_before"].fillna(False)
    lines += ["", f"Datasets whose cached sample table was missing at least one of "
                  f"`{'`, `'.join(REQUIRED_COVARIATES)}` before this pass: "
                  f"**{int(stale.sum())}**. Those are the datasets on which "
                  "`isograph_vae_residual` was silently reducing to `isograph_vae`.", ""]
    if len(bad):
        lines += ["## Datasets that did not reproduce", "",
                  "These block the refresh: their cached expression data is not recoverable "
                  "from the committed config. Per the archival decision above they are "
                  "**retained as archived artifacts**, not regenerated -- regenerating would "
                  "replace the data every published number was computed on with a different "
                  "draw. Their sample tables keep the old schema, so `isograph_vae_residual` "
                  "must not be run on them (it would be a silent no-op); see the comment at "
                  "the top of `configs/synthetic_grid.yaml`.", "",
                  "| dataset_id | scenario | seed | n mismatched | artifacts |",
                  "|---|---|---|---|---|"]
        for r in bad.head(50).itertuples():
            lines.append(f"| `{r.dataset_id}` | {r.scenario} | {r.seed} | "
                         f"{r.n_mismatched_artifacts} | `{r.mismatched}` |")
    else:
        lines.append("Every materialised dataset reproduced all of its shared artifacts "
                     "exactly -- count matrices included -- so the cached grid is "
                     "reproducible from the committed config and the refresh is safe.")
    newart = out["n_new_artifacts"].fillna(0)
    if (newart > 0).any():
        lines += ["", f"**Generator additions:** {int((newart > 0).sum())} dataset(s) "
                      "predate one or more artifacts the generator now emits (for example "
                      "`truth_switch_event`). Nothing cached differs; those artifacts are "
                      "simply absent. This tool writes only `samples.parquet` and does not "
                      "backfill them -- that is a separate decision, since downstream "
                      "metrics may or may not require them."]
    (outdir / "SAMPLE_TABLE_REFRESH.md").write_text("\n".join(lines) + "\n")
    print(f"{len(out)} datasets -> {counts.to_dict()}")
    return 1 if len(bad) else 0


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--probe", action="store_true",
                   help="rebuild one stale dataset, report reproducibility, write nothing")
    p.add_argument("--collect", action="store_true",
                   help="merge shard outputs into the report")
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--n-shards", type=int, default=1)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--dry-run", action="store_true", help="verify only; write no files")
    p.add_argument("--dataset-root", default="/".join(DATASET_ROOT))
    args = p.parse_args()

    root = rel(*args.dataset_root.split("/"))
    if args.collect:
        sys.exit(collect())
    if args.probe:
        sys.exit(0 if _report_probe(probe(root, _unique_datasets())) else 2)
    sys.exit(run(args.shard, args.n_shards, root, not args.dry_run, args.limit))


if __name__ == "__main__":
    main()
