"""Composition-adjustment rollup: gene-level turnover, not a ratio of the two counts.

Adjustment both drops and adds switch-unique genes, so the adjusted set is not nested in
the unadjusted one. The rollup therefore reports `n_overlap`/`n_new` and no retained
fraction; this pins that down against a case where a count ratio would mislead.
"""
import json

import pandas as pd
import pytest

cc = pytest.importorskip("isograph_benchmark.real_data.celltype_composition")


def _write_arm(art, subdir, unique, other):
    d = art / subdir
    d.mkdir(parents=True)
    genes = list(unique) + list(other)
    cats = ["composition_unique"] * len(unique) + ["neither"] * len(other)
    pd.DataFrame({"gene_id": genes, "category": cats}).to_parquet(d / "gene_level.parquet")
    (d / "summary.json").write_text(json.dumps({"gene_level": {
        "n_tested": len(genes), "composition_unique": len(unique),
        "abundance_unique": 0, "both": 0, "neither": len(other)}}))


def test_overlap_and_new_replace_a_count_ratio(tmp_path, monkeypatch):
    art = tmp_path / "isograph_vae"
    # base {A,B,C} -> adj {B,C,D}: equal counts, but only two genes persist and one is new.
    _write_arm(art, "incremental_association", ["A", "B", "C"], ["X", "Y"])
    _write_arm(art, "incremental_association_composition", ["B", "C", "D"], ["X", "Y"])
    monkeypatch.setattr(cc, "_artifact_dir", lambda *a, **k: art)

    row = cc._contrast_rows([("gtex-aging", "cortex", "GTEx cortex")], "standard").iloc[0]
    assert row["comp_unique_base"] == 3 and row["comp_unique_adj"] == 3
    assert row["n_overlap"] == 2 and row["n_new"] == 1
    assert "retained_frac" not in row.index


def test_missing_arm_is_skipped(tmp_path, monkeypatch):
    art = tmp_path / "isograph_vae"
    _write_arm(art, "incremental_association", ["A"], ["X"])
    monkeypatch.setattr(cc, "_artifact_dir", lambda *a, **k: art)
    assert cc._contrast_rows([("gtex-aging", "cortex", "GTEx cortex")], "standard").empty
