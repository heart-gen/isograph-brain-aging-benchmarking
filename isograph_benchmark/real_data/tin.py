"""Per-transcript TIN (RSeQC transcript integrity number) extraction for the
real-data stability test.

RSeQC writes one ``<RNum>_Aligned.sortedByCoord.out.tin.xls`` per sample, each a
tab-separated table of per-transcript TIN (column ``geneID`` actually holds the
transcript id; column ``TIN`` the score). This module collapses those per-sample
files into a single transcript x sample matrix aligned to a dataset bundle, cached
as parquet, plus a per-sample median-TIN sidecar.

Consumers:
  * ``load_tin_aligned`` returns the matrix reindexed to a requested transcript and
    sample order (missing entries filled with the per-sample median = neutral), for
    IsoGraph's ``transcript_tin`` argument (source='tin_differential').
  * the per-sample median sidecar can be added as a residualization covariate.

TIN==0 marks transcripts RSeQC could not score (insufficient coverage); we treat it
as missing rather than "fully degraded" so it contributes no spurious differential.

Heavy parse (hundreds of files x ~388k rows) -> run ``extract`` on SLURM.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from isograph.io.artifacts import load_dataset_bundle

from isograph_benchmark.real_data.stability import COHORTS, ensure_dir, rel

# Region -> directory of <RNum>_*.tin.xls files. Extend as more regions arrive.
TIN_DIRS = {
    ("brainseq", "caudate"): Path(
        "/ocean/projects/bio260021p/shared/resources/processed-data/log-files/caudate/tin"
    ),
}

_RNUM = re.compile(r"(R\d+)_")


def _cache_path(cohort: str, region: str) -> Path:
    return ensure_dir(rel("inputs", "tin")) / f"{cohort}__{region}__transcript_tin.parquet"


def _median_path(cohort: str, region: str) -> Path:
    return ensure_dir(rel("inputs", "tin")) / f"{cohort}__{region}__sample_median_tin.csv"


def _tin_files_by_rnum(tin_dir: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for p in tin_dir.iterdir():  # iterdir (not glob) on Ocean FS
        if not p.name.endswith(".tin.xls"):
            continue
        m = _RNUM.match(p.name)
        if m:
            out.setdefault(m.group(1), p)  # first file wins on the rare duplicate RNum
    return out


def extract(cohort: str, region: str) -> None:
    """Parse per-sample TIN files into a bundle-aligned transcript x sample matrix."""
    key = (cohort, region)
    if key not in TIN_DIRS:
        raise SystemExit(f"no TIN directory registered for {key}")
    spec = COHORTS[cohort]
    bundle = load_dataset_bundle(rel(*spec["bundle_root"], region))
    tt = bundle.feature_tables["transcript"]
    transcript_ids = tt["transcript_id"].astype(str).to_numpy()
    sample_ids = bundle.sample_table["sample_id"].astype(str).to_numpy()
    tx_index = {tid: i for i, tid in enumerate(transcript_ids)}

    files = _tin_files_by_rnum(TIN_DIRS[key])
    present = [s for s in sample_ids if s in files]
    missing = [s for s in sample_ids if s not in files]
    print(f"[{cohort}/{region}] bundle: {len(sample_ids)} samples, "
          f"{len(transcript_ids)} transcripts | TIN files matched: {len(present)} "
          f"| missing: {len(missing)}", flush=True)
    if missing:
        print(f"  WARNING missing TIN for {len(missing)} samples: "
              f"{missing[:8]}{'...' if len(missing) > 8 else ''}", flush=True)

    mat = np.full((len(transcript_ids), len(sample_ids)), np.nan, dtype=np.float32)
    col_of = {s: j for j, s in enumerate(sample_ids)}
    for n, s in enumerate(present, 1):
        df = pd.read_csv(files[s], sep="\t", usecols=["geneID", "TIN"])
        tin = df["TIN"].to_numpy(dtype=np.float32)
        tin[tin == 0.0] = np.nan  # 0 = unscored (insufficient coverage) -> missing
        rows = df["geneID"].astype(str).map(tx_index)
        keep = rows.notna().to_numpy()
        mat[rows[keep].to_numpy(dtype=int), col_of[s]] = tin[keep]
        if n % 25 == 0 or n == len(present):
            print(f"  parsed {n}/{len(present)} ({s})", flush=True)

    # per-sample median over scored transcripts (neutral fill value for downstream)
    sample_median = np.nanmedian(mat, axis=0)
    sample_median = np.where(np.isfinite(sample_median), sample_median,
                             float(np.nanmedian(sample_median)))

    out = pd.DataFrame(mat, columns=sample_ids)
    out.insert(0, "transcript_id", transcript_ids)
    out.to_parquet(_cache_path(cohort, region), index=False, compression="zstd")
    pd.DataFrame({"sample_id": sample_ids, "median_tin": sample_median}).to_csv(
        _median_path(cohort, region), index=False)
    cov = float(np.mean(np.isfinite(mat)))
    print(f"  wrote {_cache_path(cohort, region)}  (scored-entry coverage {cov:.3f})",
          flush=True)
    print(f"  wrote {_median_path(cohort, region)}", flush=True)


def load_tin_aligned(
    cohort: str, region: str, transcript_ids, sample_ids
) -> np.ndarray:
    """Return the cached TIN matrix reindexed to the given transcript/sample order.

    Missing/unscored entries are filled with the per-sample median (neutral: adds no
    differential-degradation signal), so the row order may be any transcript subset
    (e.g. after expression filtering) and the columns any sample subset (a split-half).
    """
    cache = _cache_path(cohort, region)
    if not cache.exists():
        raise SystemExit(f"TIN cache missing: {cache}; run `tin extract --cohort "
                         f"{cohort} --region {region}` first")
    df = pd.read_parquet(cache).set_index("transcript_id")
    df = df.reindex(index=[str(t) for t in transcript_ids],
                    columns=[str(s) for s in sample_ids])
    mat = df.to_numpy(dtype=np.float32)
    col_median = np.nanmedian(mat, axis=0)
    col_median = np.where(np.isfinite(col_median), col_median,
                          float(np.nanmedian(col_median)) if np.isfinite(np.nanmedian(col_median)) else 0.0)
    inds = np.where(np.isnan(mat))
    mat[inds] = np.take(col_median, inds[1])
    return mat


def load_sample_median_tin(cohort: str, region: str, sample_ids) -> np.ndarray:
    """Return per-sample median TIN reindexed to ``sample_ids`` (neutral fill = global
    median), for use as a residualization covariate."""
    path = _median_path(cohort, region)
    if not path.exists():
        raise SystemExit(f"median-TIN sidecar missing: {path}; run `tin extract` first")
    s = pd.read_csv(path).set_index("sample_id")["median_tin"]
    vals = s.reindex([str(x) for x in sample_ids]).to_numpy(dtype=float)
    if np.isnan(vals).any():
        vals = np.nan_to_num(vals, nan=float(np.nanmedian(vals)))
    return vals


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    ex = sub.add_parser("extract", help="parse per-sample TIN files -> bundle-aligned parquet")
    ex.add_argument("--cohort", required=True, choices=sorted({c for c, _ in TIN_DIRS}))
    ex.add_argument("--region", required=True)
    args = ap.parse_args()
    if args.cmd == "extract":
        extract(args.cohort, args.region)


if __name__ == "__main__":
    main()
