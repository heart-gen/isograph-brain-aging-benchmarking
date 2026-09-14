"""Run telemetry must identify the code that produced a result, not just a package version."""
import re
from pathlib import Path

from isograph_benchmark.benchmark import telemetry


def test_software_versions_records_both_commits():
    v = telemetry.software_versions()
    for key in ("isograph_commit", "benchmark_commit"):
        assert key in v
        assert v[key] is None or re.fullmatch(r"[0-9a-f]{40}", v[key])
    assert v["benchmark_commit"] is not None  # the test suite runs inside this repository


def test_git_state_outside_a_repository_is_none(tmp_path: Path):
    assert telemetry._git_state(tmp_path) == (None, None)
