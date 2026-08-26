from __future__ import annotations

import gzip

import pandas as pd
import pytest

from isograph_benchmark.real_data.neuronal_clip_fetch import (
    _bed12_event_stats,
    _candidate_identity_sha256,
    _guard_candidate_identity,
    _quantas_event_stats,
)


def _candidates(module_id: str = "M008") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "region": ["frontal_cortex_ba9", "frontal_cortex_ba9"],
            "module_id": [module_id, "M010"],
            "rbp": ["TARDBP", "NOVA2"],
            "gene": ["ENSG1", "ENSG2"],
            "transcript_id_1": ["ENST1", "ENST3"],
            "transcript_id_2": ["ENST2", "ENST4"],
        }
    )


def test_candidate_identity_is_order_invariant() -> None:
    frame = _candidates()
    assert _candidate_identity_sha256(frame) == _candidate_identity_sha256(
        frame.iloc[::-1].reset_index(drop=True)
    )


def test_candidate_identity_changes_when_membership_changes_but_count_does_not() -> None:
    # The failure mode the count guard cannot see: same number of nominations,
    # different nominations.
    before = _candidates("M008")
    after = _candidates("M014")
    assert len(before) == len(after)
    assert _candidate_identity_sha256(before) != _candidate_identity_sha256(after)


def test_guard_raises_on_changed_identity() -> None:
    observed = _candidate_identity_sha256(_candidates())
    with pytest.raises(RuntimeError, match="candidate identity changed"):
        _guard_candidate_identity(
            {"expected_candidate_identity_sha256": "0" * 64}, observed
        )


def test_guard_passes_when_pinned_identity_matches() -> None:
    observed = _candidate_identity_sha256(_candidates())
    _guard_candidate_identity(
        {"expected_candidate_identity_sha256": observed}, observed
    )


def test_guard_is_permissive_when_unpinned() -> None:
    _guard_candidate_identity({}, _candidate_identity_sha256(_candidates()))


def test_quantas_event_stats_handles_crlf_and_counts_unique_ids(tmp_path) -> None:
    path = tmp_path / "events.txt.gz"
    with gzip.open(path, "wt", newline="") as handle:
        handle.write(
            "AS.name\tFDR\tdelta.exon.inclusion.rate(dI)\r\n"
            "event_1\t0.01\t0.2\r\n"
            "event_1\t0.02\t0.3\r\n"
            "event_2\t0.03\t-0.4\r\n"
        )
    assert _quantas_event_stats(path) == {
        "event_rows": 3,
        "unique_event_ids": 2,
    }


def test_bed12_event_stats_allows_multiple_coordinates_per_event(tmp_path) -> None:
    path = tmp_path / "events.bed.gz"
    rows = [
        "chr1\t10\t30\tevent_1\t0\t+\t10\t30\t0\t2\t5,5\t0,15",
        "chr1\t40\t60\tevent_1\t0\t+\t40\t60\t0\t2\t5,5\t0,15",
        "chr2\t70\t90\tevent_2\t0\t-\t70\t90\t0\t2\t5,5\t0,15",
    ]
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(rows) + "\n")
    assert _bed12_event_stats(path) == {
        "coordinate_rows": 3,
        "coordinate_event_ids": 2,
    }
