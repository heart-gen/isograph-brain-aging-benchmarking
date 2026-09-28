"""Marker-test rollup behind the composition supplement's marker panel."""
import pandas as pd
import pytest

cac = pytest.importorskip("isograph_benchmark.real_data.composition_age_coupling")


def _store(root, cohort, region, fdr_enrich, fdr_deplete):
    d = root / cohort / region / "_m" / "isograph_vae" / "celltype_composition"
    d.mkdir(parents=True)
    pd.DataFrame({"fdr_enrich": fdr_enrich, "fdr_deplete": fdr_deplete}).to_parquet(
        d / "marker_enrichment.parquet")


def test_marker_counts_carry_their_own_denominator(tmp_path):
    # Two regions with the same single significant test but very different test counts.
    # The rollup has to keep the number of tests, or 1 of 3 and 1 of 5 read as the same.
    _store(tmp_path, "brainseq", "caudate", [0.01, 0.5, 0.9], [1.0, 1.0, 1.0])
    _store(tmp_path, "gtex", "amygdala", [0.9, 0.9, 0.9, 0.9, 0.01], [1.0, 0.02, 1.0, 1.0, 1.0])
    out = cac._marker_tests("default", root=tmp_path).set_index("label")

    assert out.loc["aging caudate", "n_marker_tests"] == 3
    assert out.loc["aging caudate", "marker_enriched"] == 1
    assert out.loc["GTEx amygdala", "n_marker_tests"] == 5
    assert out.loc["GTEx amygdala", "marker_depleted"] == 1


def test_labels_join_the_composition_rollup(tmp_path):
    # The marker table is read beside the composition ledgers, so it has to use the same
    # analysis labels, including the hand-named BrainSEQ ones and the disease analysis.
    _store(tmp_path, "brainseq", "caudate_sczd", [0.9], [0.9])
    _store(tmp_path, "brainseq", "dlpfc", [0.9], [0.9])
    out = cac._marker_tests("default", root=tmp_path)
    assert set(out["label"]) == {"SCZD (caudate)", "aging DLPFC"}
    assert set(out["class"]) == {"Disease (SCZD)", "Cortical"}
    assert set(out["cohort"]) == {"BrainSEQ"}


def test_no_marker_files_gives_an_empty_table_with_its_columns(tmp_path):
    out = cac._marker_tests("default", root=tmp_path)
    assert out.empty
    assert list(out.columns) == ["label", "cohort", "class", "n_marker_tests",
                                 "marker_enriched", "marker_depleted"]
