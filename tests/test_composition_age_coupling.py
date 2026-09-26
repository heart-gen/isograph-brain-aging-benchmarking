"""Composition-vs-age coupling: Spearman per analysis x cell type, BH within analysis."""
import numpy as np
import pandas as pd
import pytest

cac = pytest.importorskip("isograph_benchmark.real_data.composition_age_coupling")


def _tidy(label, cohort, cell, prop, age):
    return pd.DataFrame({"label": label, "cohort": cohort, "cell_type": cell,
                         "proportion": prop, "age": age})


def test_monotone_cell_type_is_detected_and_flat_one_is_not():
    age = np.arange(20, 80, dtype=float)
    tidy = pd.concat([
        _tidy("GTEx cortex", "GTEx", "Astro", age / 200, age),      # perfectly monotone
        _tidy("GTEx cortex", "GTEx", "Oligo", np.tile([0.1, 0.2], 30), age),
    ], ignore_index=True)
    out = cac._coupling(tidy).set_index("cell_type")

    assert out.loc["Astro", "rho_age"] == pytest.approx(1.0)
    assert out.loc["Astro", "fdr_age"] < 0.05
    assert abs(out.loc["Oligo", "rho_age"]) < 0.3
    assert out.loc["Astro", "n_samples"] == len(age)


def test_zero_variance_cell_type_reports_missing_rather_than_dropping_the_row():
    age = np.arange(20, 50, dtype=float)
    tidy = pd.concat([
        _tidy("GTEx cortex", "GTEx", "Astro", age / 200, age),
        _tidy("GTEx cortex", "GTEx", "Mural", np.zeros(len(age)), age),
    ], ignore_index=True)
    out = cac._coupling(tidy).set_index("cell_type")

    # The cell type was still adjusted for, so it must stay visible in the panel.
    assert "Mural" in out.index
    assert np.isnan(out.loc["Mural", "rho_age"]) and np.isnan(out.loc["Mural", "fdr_age"])
    # ... and must not enter the BH family, which would deflate the other test.
    assert out.loc["Astro", "fdr_age"] == pytest.approx(out.loc["Astro", "p_age"])


def test_fdr_is_computed_within_analysis_not_across_them():
    age = np.arange(20, 80, dtype=float)
    rng = np.random.default_rng(13)
    noise = [_tidy("GTEx cortex", "GTEx", f"C{i}", rng.random(len(age)), age)
             for i in range(8)]
    tidy = pd.concat(noise + [_tidy("GTEx amygdala", "GTEx", "Astro", age / 200, age)],
                     ignore_index=True)
    out = cac._coupling(tidy).set_index(["label", "cell_type"])

    # The lone perfect correlation is its own family of one, so BH cannot inflate it.
    row = out.loc[("GTEx amygdala", "Astro")]
    assert row["fdr_age"] == pytest.approx(row["p_age"])
