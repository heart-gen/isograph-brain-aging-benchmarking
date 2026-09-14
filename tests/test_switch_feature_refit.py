"""Setting grid and partition comparison for the full-refit sensitivity."""
import numpy as np
import pandas as pd

from isograph_benchmark.real_data import switch_feature_refit as sfr


def test_grid_has_the_published_setting_once_and_first():
    # published + 4 pseudocounts + 6 transcript-filter alternatives + 3 minor-isoform thresholds
    for cohort, n in (("brainseq", 14), ("gtex", 14)):
        grid = sfr.setting_grid(cohort)
        assert grid[0]["axis"] == "published"
        assert len(grid) == n
        assert len({sfr.slug(e) for e in grid}) == n


def test_counts_for_setting_shifts_the_pseudocount_and_skips_an_absent_filter():
    counts = np.array([[0.0, 10.0], [5.0, 0.0]])
    table = pd.DataFrame({"gene_id": ["g", "g"], "transcript_id": ["t1", "t2"]})
    setting = {"pseudocount": 1.0, "transcript_filter": "none", "min_tx_prop": 0.0,
               "min_tx_fraction": 0.0, "min_usage": 0.0}
    tc, tt = sfr.counts_for_setting(counts, table, setting)
    assert tc.shape == counts.shape
    assert np.allclose(tc, counts + 0.5)


def test_identical_partitions_agree_perfectly_and_retain_every_association():
    mods = pd.DataFrame({"gene_id": [f"g{i}" for i in range(6)],
                         "module_id": ["M000"] * 3 + ["M001"] * 3})
    age = pd.DataFrame({"module_id": ["M000", "M001"], "effect": [0.4, -0.3], "fdr": [0.01, 0.2]})
    out = sfr.compare_partitions(mods, mods, age, age)
    assert out["ari"] == 1.0 and out["nmi"] == 1.0
    assert out["n_age_sig_published"] == 1 and out["n_age_sig_retained"] == 1


def test_a_relabelled_refit_is_matched_by_content_not_id():
    pub = pd.DataFrame({"gene_id": [f"g{i}" for i in range(6)],
                        "module_id": ["M000"] * 3 + ["M001"] * 3})
    ref = pub.assign(module_id=pub["module_id"].map({"M000": "M001", "M001": "M000"}))
    pub_age = pd.DataFrame({"module_id": ["M000", "M001"], "effect": [0.4, -0.3], "fdr": [0.01, 0.2]})
    ref_age = pd.DataFrame({"module_id": ["M001", "M000"], "effect": [0.4, -0.3], "fdr": [0.01, 0.2]})
    out = sfr.compare_partitions(pub, ref, pub_age, ref_age)
    assert out["ari"] == 1.0
    assert out["n_age_sig_retained"] == 1
