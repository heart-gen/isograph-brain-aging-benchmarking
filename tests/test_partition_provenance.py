"""The guard must fire on the exact failure that scrambled the 2026-06-29 GO partition.

A refit re-assigns Leiden module ids, so a stale enrichment table joined on
module_id hands each module some other module's pheno_fdr / n_go_terms. Without a
guard that is silent. These tests pin both tiers of detection and, importantly,
pin that a CORRECT pairing still passes — a guard that rejects valid inputs would
just be routed around.
"""
from __future__ import annotations

import pandas as pd
import pytest

from isograph_benchmark.real_data.partition_provenance import (
    StalePartitionError,
    check_partition,
    load_enrichment,
    partition_fingerprint,
    read_fingerprint,
    write_with_fingerprint,
)


def _modules(sizes: dict[str, int]) -> pd.DataFrame:
    rows = []
    g = 0
    for mid, n in sizes.items():
        for _ in range(n):
            rows.append({"gene_id": f"ENSG{g:08d}", "module_id": mid})
            g += 1
    return pd.DataFrame(rows)


def _enrich(modules: pd.DataFrame, **cols) -> pd.DataFrame:
    counts = modules["module_id"].value_counts()
    ids = sorted(counts.index)
    return pd.DataFrame({
        "module_id": ids,
        "n_genes": [int(counts[m]) for m in ids],
        "pheno_fdr": [0.01] * len(ids),
        "n_go_terms": [0] * len(ids),
        **cols,
    })


def test_matching_pair_passes():
    m = _modules({"M000": 10, "M001": 5, "M002": 3})
    check_partition(m, _enrich(m), context="test")


def test_relabelled_fit_is_rejected_by_gene_counts():
    """The real failure: same ids, different gene sets after a refit."""
    old = _modules({"M000": 10, "M001": 5, "M002": 3})
    new = _modules({"M000": 12, "M001": 4, "M002": 2})
    with pytest.raises(StalePartitionError, match="disagree on gene count"):
        check_partition(new, _enrich(old), context="test")


def test_module_added_by_refit_is_rejected():
    """A refit that finds more modules leaves orphans the stale table cannot describe."""
    old = _modules({"M000": 10, "M001": 5})
    new = _modules({"M000": 10, "M001": 5, "M002": 4})
    with pytest.raises(StalePartitionError, match="module-id sets differ"):
        check_partition(new, _enrich(old), context="test")


def test_fingerprint_is_order_and_index_independent():
    m = _modules({"M000": 4, "M001": 3})
    shuffled = m.sample(frac=1.0, random_state=0).reset_index(drop=True)
    assert partition_fingerprint(m) == partition_fingerprint(shuffled)


def test_fingerprint_changes_when_assignment_changes():
    m = _modules({"M000": 4, "M001": 3})
    moved = m.copy()
    moved.loc[0, "module_id"] = "M001"
    assert partition_fingerprint(m) != partition_fingerprint(moved)


def test_fingerprint_catches_a_size_preserving_relabel(tmp_path):
    """Tier 2 cannot see a permutation that preserves every group size. Tier 1 can.

    This is the case that motivates stamping the fingerprint rather than relying on
    the structural check alone.
    """
    m = _modules({"M000": 5, "M001": 5})
    swapped = m.copy()
    swapped["module_id"] = swapped["module_id"].map({"M000": "M001", "M001": "M000"})

    path = tmp_path / "isograph_modules.parquet"
    write_with_fingerprint(_enrich(m), path, m)

    check_partition(m, pd.read_parquet(path), context="test", enrich_path=path)
    with pytest.raises(StalePartitionError, match="fingerprint mismatch"):
        check_partition(swapped, pd.read_parquet(path), context="test", enrich_path=path)


def test_fingerprint_round_trips_through_parquet(tmp_path):
    m = _modules({"M000": 3, "M001": 2})
    path = tmp_path / "e.parquet"
    sha = write_with_fingerprint(_enrich(m), path, m)
    assert read_fingerprint(path) == sha == partition_fingerprint(m)
    assert list(pd.read_parquet(path)["module_id"]) == ["M000", "M001"]


def test_legacy_table_without_fingerprint_still_structurally_checked(tmp_path):
    """Tables written before the guard carry no stamp; tier 2 must still apply."""
    old = _modules({"M000": 10, "M001": 5})
    new = _modules({"M000": 11, "M001": 4})
    path = tmp_path / "legacy.parquet"
    _enrich(old).to_parquet(path, index=False)
    assert read_fingerprint(path) is None
    with pytest.raises(StalePartitionError):
        load_enrichment(path, new, context="test")


def test_missing_table_returns_none_unless_required(tmp_path):
    m = _modules({"M000": 2})
    assert load_enrichment(tmp_path / "absent.parquet", m, context="test") is None
    with pytest.raises(FileNotFoundError):
        load_enrichment(tmp_path / "absent.parquet", m, context="test", required=True)


def test_error_message_names_the_caller():
    old = _modules({"M000": 10})
    new = _modules({"M000": 11})
    with pytest.raises(StalePartitionError, match=r"qtl_anchoring gtex-aging/cortex"):
        check_partition(new, _enrich(old), context="qtl_anchoring gtex-aging/cortex")
