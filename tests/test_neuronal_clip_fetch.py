from __future__ import annotations

import gzip

from isograph_benchmark.real_data.neuronal_clip_fetch import (
    _bed12_event_stats,
    _quantas_event_stats,
)


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
