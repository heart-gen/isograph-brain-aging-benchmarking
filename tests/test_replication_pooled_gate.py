"""Tests for the pooled cross-cohort replication arm's gate and its empty-result behaviour.

Both parquets shipped at 0 rows for months and no consumer noticed, because a 0-row result
table is indistinguishable from "nothing replicates". These pin the two decisions that fixed
it: the gate no longer selects on module granularity, and a zero result is never written as
a readable table.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest

from isograph_benchmark.real_data import module_trust as mt


def _argparse_default(cmd: str, flag: str):
    """Default for `flag` on subcommand `cmd`, read from the real parser."""
    import argparse
    import contextlib
    import io

    parser_defaults = {}
    real_add_parser = argparse._SubParsersAction.add_parser

    def spy(self, name, **kw):
        p = real_add_parser(self, name, **kw)
        parser_defaults[name] = p
        return p

    argparse._SubParsersAction.add_parser = spy
    try:
        with contextlib.redirect_stderr(io.StringIO()), pytest.raises(SystemExit):
            mt.main()  # no argv -> exits after the parser is built
    finally:
        argparse._SubParsersAction.add_parser = real_add_parser
    p = parser_defaults[cmd]
    (action,) = [a for a in p._actions if flag in a.option_strings]
    return action.default


def test_pooled_gate_does_not_select_on_granularity():
    """The 0.25 default emptied the arm by demanding a Jaccard IsoGraph's partition cannot reach.

    Cross-cohort module Jaccard rises with module size, so a threshold gate admits WGCNA's few
    giant modules (max 0.65 on these data) and excludes IsoGraph's fine ones (max 0.12) --
    selecting on granularity rather than on replication.
    """
    assert _argparse_default("replication-pooled", "--min-jaccard") == 0.0


def test_pooled_gate_still_filters_when_asked(tmp_path, monkeypatch):
    """The old 0.25 behaviour must remain reachable as a sensitivity arm, not be deleted."""
    monkeypatch.setattr(mt, "_out_dir", lambda: tmp_path)
    rows = []
    for i in range(8):
        rows.append({
            "pair": "caudate", "bs_module": f"M{i:03d}", "gtex_match": f"M{i:03d}",
            # half above the old gate, half below it
            "gene_jaccard": 0.40 if i % 2 else 0.10,
            "driver_jaccard": 0.0, "match_trusted": True,
            "age_effect_bs": 0.5 + i * 0.1, "age_p_bs": 0.01,
            "age_effect_gtex": 0.4 + i * 0.1, "age_p_gtex": 0.02,
            "sign_match": True, "both_sig": True, "replicates": True,
        })
    monkeypatch.setattr(mt, "_crosscohort_rows", lambda *a, **k: rows)

    mt.replication_pooled("isograph", k=5, sig=0.05, min_jaccard=0.25, n_perm=10, seed=13)

    out = pd.read_parquet(tmp_path / "module_aging_replication_pooled__isograph.parquet")
    assert len(out) == 4 * len(mt.REGION_PAIRS)
    assert (out["gene_jaccard"] >= 0.25).all()


def test_zero_result_writes_no_parquet_but_records_the_gate(tmp_path, monkeypatch):
    """An impossible gate must leave no readable result table, and must say why in the json."""
    monkeypatch.setattr(mt, "_out_dir", lambda: tmp_path)

    rows = [{
        "pair": "caudate", "bs_module": "M001", "gtex_match": "M002",
        "gene_jaccard": 0.10, "driver_jaccard": 0.0, "match_trusted": True,
        "age_effect_bs": 0.5, "age_p_bs": 0.01,
        "age_effect_gtex": 0.4, "age_p_gtex": 0.02,
        "sign_match": True, "both_sig": True, "replicates": True,
    }]
    monkeypatch.setattr(mt, "_crosscohort_rows", lambda *a, **k: rows)

    stale = tmp_path / "module_aging_replication_pooled__isograph.parquet"
    pd.DataFrame(rows).to_parquet(stale, index=False)

    with pytest.raises(SystemExit):
        mt.replication_pooled("isograph", k=5, sig=0.05, min_jaccard=0.99,
                              n_perm=10, seed=13)

    assert not stale.exists(), "a stale result table survived a zero run"
    stats = json.loads(
        (tmp_path / "module_aging_replication_pooled__isograph__stats.json").read_text())
    assert stats["n_pairs_reproducible"] == 0
    assert stats["min_jaccard"] == 0.99
    assert stats["max_gene_jaccard"] == pytest.approx(0.10)


def test_nonzero_result_writes_both_files(tmp_path, monkeypatch):
    monkeypatch.setattr(mt, "_out_dir", lambda: tmp_path)
    rows = []
    for i in range(8):
        rows.append({
            "pair": "caudate", "bs_module": f"M{i:03d}", "gtex_match": f"M{i:03d}",
            "gene_jaccard": 0.10, "driver_jaccard": 0.0, "match_trusted": True,
            # vary the effects so the Spearman arm is defined
            "age_effect_bs": 0.5 + i * 0.1, "age_p_bs": 0.01,
            "age_effect_gtex": 0.4 + i * 0.1, "age_p_gtex": 0.02,
            "sign_match": True, "both_sig": True, "replicates": True,
        })
    monkeypatch.setattr(mt, "_crosscohort_rows", lambda *a, **k: rows)

    mt.replication_pooled("isograph", k=5, sig=0.05, min_jaccard=0.0, n_perm=10, seed=13)

    out = pd.read_parquet(tmp_path / "module_aging_replication_pooled__isograph.parquet")
    # pooling concatenates every region pair, which is the whole point of the arm
    n_expected = 8 * len(mt.REGION_PAIRS)
    assert len(out) == n_expected
    stats = json.loads(
        (tmp_path / "module_aging_replication_pooled__isograph__stats.json").read_text())
    assert stats["n_pairs_reproducible"] == n_expected
