#!/usr/bin/env python
"""Audit every committed (fit, module_enrichment) pair for the stale-join failure.

Run this before quoting any module-set number, and before submission. It is the
check that would have caught the 2026-06-29 scramble in seconds: `qtl_anchoring`
was re-run against a refreshed fit but the pre-refresh enrichment table, and
because Leiden module ids are re-assigned on every fit, the GO-invisible /
GO-visible partition it produced was a relabelling rather than a stale copy.

Keys strictly on (cohort, region). An earlier ad-hoc version of this scan keyed on
region name alone and silently merged BrainSEQ hippocampus (50 modules) with GTEx
hippocampus (41), reporting 54 false positives that all traced to the 9-module
difference. Region names are NOT unique across cohorts; do not "simplify" this.

Exit status is 1 if anything is stale, so it can gate CI or a submission script.

  python scripts/audit_partition_provenance.py
  python scripts/audit_partition_provenance.py --verbose
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Run directly (scripts/ is on sys.path, the repo root is not), so locate the root
# by its .here marker the same way isograph_benchmark.paths does.
_d = Path(__file__).resolve().parent
while not (_d / ".here").exists() and _d != _d.parent:
    _d = _d.parent
sys.path.insert(0, str(_d))

from isograph_benchmark.paths import COHORTS, cohort_dir, region_store  # noqa: E402
from isograph_benchmark.real_data.partition_provenance import (  # noqa: E402
    ENRICH_PREFIX_FIT_DIR,
    StalePartitionError,
    check_partition,
    read_fingerprint,
)

import pandas as pd  # noqa: E402


def iter_pairs():
    """Yield (cohort, region, prefix, enrich_path, fit_path) for every stored pair."""
    for cohort in COHORTS:
        base = cohort_dir(cohort)
        if not base.exists():
            continue
        # Ocean FS: iterdir, never a recursive glob over the artifact store.
        for region_dir in sorted(p for p in base.iterdir() if p.is_dir()):
            region = region_dir.name
            enrich_dir = region_store(cohort, region, "module_enrichment")
            if not enrich_dir.exists():
                continue
            for prefix, fit_subdir in ENRICH_PREFIX_FIT_DIR.items():
                enrich_path = enrich_dir / f"{prefix}_modules.parquet"
                fit_path = region_store(cohort, region, fit_subdir, "modules.parquet")
                if enrich_path.exists() and fit_path.exists():
                    yield cohort, region, prefix, enrich_path, fit_path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verbose", action="store_true", help="list every pair, not just failures")
    args = ap.parse_args()

    checked = stale = unstamped = 0
    failures: list[str] = []
    for cohort, region, prefix, enrich_path, fit_path in iter_pairs():
        checked += 1
        label = f"{cohort}/{region} [{prefix}]"
        stamped = read_fingerprint(enrich_path) is not None
        unstamped += not stamped
        try:
            check_partition(
                pd.read_parquet(fit_path),
                pd.read_parquet(enrich_path),
                context=label,
                enrich_path=enrich_path,
            )
        except StalePartitionError as exc:
            stale += 1
            failures.append(f"{label}: {exc}")
            print(f"STALE       {label}")
        else:
            if args.verbose:
                tier = "fingerprint" if stamped else "structural only"
                print(f"ok          {label}  ({tier})")

    print(f"\n{checked} (fit, enrichment) pairs checked; {stale} stale.")
    if unstamped:
        print(
            f"{unstamped} table(s) carry no fit fingerprint and were checked "
            f"structurally only — re-run module_enrichment to stamp them, which "
            f"upgrades the check from 'group sizes agree' to 'same partition'."
        )
    if failures:
        print("\nDetail:")
        for f in failures:
            print(f"  - {f}")
        print("\nFix: re-run module_enrichment for the region(s) above, then re-run "
              "every analysis that consumes the enrichment table (qtl_anchoring, "
              "module_genetic_anchoring, sqtl_concordance, coloc_prep, "
              "go_invisible_gate, rbp_regulon, rbp_target_panel, replication_*).")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
