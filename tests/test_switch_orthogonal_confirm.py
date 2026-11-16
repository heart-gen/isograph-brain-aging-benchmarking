"""The long-read anchored check can read either coloc layer without mixing them."""
import pandas as pd
import pytest

from isograph_benchmark.real_data.switch_orthogonal_confirm import (
    EVIDENCE_COL,
    anchored_pairs,
    out_path,
)


def _events(ev_col: str, support: list[tuple[float, str]]) -> pd.DataFrame:
    return pd.DataFrame([
        {"analysis": f"aging__{trait}", "trait": trait, "gene": "ENSG1", "gene_name": "G1",
         "kind": "sQTL", "tissue": "Brain_Cortex", "iso_region": "cortex", "best_rsid": None,
         "junction": "chr1:10-20(+)", "go_invisible": False,
         "switch_pair": "ENST1.1 | ENST2.3 | ENST3.1", "junction_transcripts": "ENST1.2",
         ev_col: value}
        for value, trait in support
    ])


def test_signal_layer_writes_its_own_subdirectory():
    clpp, signal = out_path("clpp"), out_path("signal")
    assert signal != clpp and signal.parent == clpp
    with pytest.raises(SystemExit):
        out_path("abf")


def test_signal_layer_scores_pairs_by_its_own_posterior():
    assert EVIDENCE_COL == {"clpp": "clpp", "signal": "PP4_sQTL"}
    p = anchored_pairs(_events("PP4_sQTL", [(0.85, "pd"), (0.97, "lbd")]), ev_col="PP4_sQTL")
    assert "clpp" not in p.columns
    assert set(p["anchored_tx"]) == {"ENST1"}          # pair members that carry the junction
    assert sorted(p["partner_tx"]) == ["ENST2", "ENST3"]
    assert (p["PP4_sQTL"] == 0.97).all()               # strongest supporting event kept
    assert (p["traits"] == "lbd,pd").all() and (p["n_events"] == 2).all()


def test_clpp_layer_is_the_default():
    p = anchored_pairs(_events("clpp", [(0.2, "ad")]))
    assert "clpp" in p.columns and len(p) == 2


def test_exception_events_swap_in_the_cross_tissue_pair_and_skip_matched_rows():
    from isograph_benchmark.real_data.switch_orthogonal_confirm import exception_events
    ev = pd.concat([_events("PP4_sQTL", [(0.96, "als")]), _events("PP4_sQTL", [(0.90, "pd")])],
                   ignore_index=True)
    ev["junction_in_switch_pair"] = [False, True]
    ev["cross_tissue_in_switch_pair"] = [True, False]
    ev["cross_tissue_switch_pair"] = ["ENST1.1 | ENST7.1", ""]
    ex = exception_events(ev)
    assert list(ex["trait"]) == ["als"]            # a tissue-matched row is never an exception
    assert ex["switch_pair"].iloc[0] == "ENST1.1 | ENST7.1"
    p = anchored_pairs(ex, ev_col="PP4_sQTL")
    assert set(p["anchored_tx"]) == {"ENST1"} and set(p["partner_tx"]) == {"ENST7"}


def test_exception_events_is_empty_without_the_exception_columns():
    from isograph_benchmark.real_data.switch_orthogonal_confirm import exception_events
    assert exception_events(_events("clpp", [(0.2, "ad")])).empty


# --------------------------------------------------------------------------- #
# Joint (multi-key) matched null
# --------------------------------------------------------------------------- #
def _null_frame(n: int, seed: int) -> pd.DataFrame:
    import numpy as np
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "max_mean_if": rng.uniform(0, 1, n),
        "detect_rate": rng.uniform(0, 1, n),
        "n_iso": rng.integers(2, 12, n).astype(float),
        "usage_sd": rng.uniform(0, 0.3, n),
        "switch_like": rng.uniform(0, 1, n) < 0.4,
        "usage_spearman": rng.uniform(-1, 1, n),
    })


def test_default_matching_is_the_abundance_decile_null():
    from isograph_benchmark.real_data.longread_null import matched_null
    focal, bg = _null_frame(40, 1), _null_frame(400, 2)
    a = matched_null(focal, bg, "switch_like_rate", 200, 13)
    b = matched_null(focal, bg, "switch_like_rate", 200, 13, keys=("max_mean_if",), bins=(10,))
    assert a == b
    assert "match_keys" not in a          # abundance-only summaries keep their shape


def test_joint_null_draws_only_from_the_matching_joint_stratum():
    from isograph_benchmark.real_data.longread_null import (
        JOINT_BINS,
        JOINT_KEYS,
        matched_null,
    )
    bg = _null_frame(3000, 3)
    # Background pairs are switch-like exactly when they are well detected; abundance
    # matching alone cannot see that, joint matching on detection rate must.
    bg["switch_like"] = bg["detect_rate"] > 2 / 3
    focal = bg[bg["detect_rate"] > 0.9].head(50).copy()
    abundance = matched_null(focal, bg, "switch_like_rate", 300, 13)
    joint = matched_null(focal, bg, "switch_like_rate", 300, 13, keys=JOINT_KEYS,
                         bins=JOINT_BINS)
    assert abundance["null_mean"] < 0.5
    assert joint["null_mean"] == pytest.approx(1.0)
    assert joint["match_keys"] == list(JOINT_KEYS) and joint["match_bins"] == [5, 3, 3, 3]
    assert joint["n_focal"] + joint["n_focal_dropped_no_stratum"] == len(focal)


def test_joint_null_excludes_pairs_with_a_missing_key():
    import numpy as np

    from isograph_benchmark.real_data.longread_null import (
        JOINT_BINS,
        JOINT_KEYS,
        matched_null,
    )
    focal, bg = _null_frame(30, 4), _null_frame(500, 5)
    focal.loc[:4, "usage_sd"] = np.nan
    r = matched_null(focal, bg, "switch_like_rate", 50, 13, keys=JOINT_KEYS, bins=JOINT_BINS)
    assert r["n_focal"] + r["n_focal_dropped_no_stratum"] == 25
    with pytest.raises(ValueError):
        matched_null(focal, bg, "switch_like_rate", 5, 13, keys=JOINT_KEYS, bins=(5, 3))


def test_add_match_keys_from_the_long_read_matrix():
    import numpy as np

    from isograph_benchmark.real_data.switch_orthogonal_confirm import add_match_keys
    tx = pd.DataFrame({"gene": ["G1", "G1", "G1", "G2", "G2"],
                       "tx": ["A", "B", "C", "D", "E"],
                       "n_detected_samples": [4, 2, 0, 4, 1]})
    frac = pd.DataFrame([[0.5, 0.7, 0.6, 0.8], [0.5, 0.3, 0.4, 0.2], [0.0] * 4,
                         [1.0, 0.0, 0.5, np.nan], [0.0, 1.0, 0.5, np.nan]],
                        columns=list("wxyz"))
    pairs = pd.DataFrame({"gene": ["G1", "G2", "G1"], "t1": ["A", "D", "A"],
                          "t2": ["B", "E", "Z"]})
    out = add_match_keys(pairs, tx, frac)
    assert out["detect_rate"].tolist()[:2] == [0.75, 0.625]
    assert np.isnan(out["detect_rate"].iloc[2])        # partner absent from long read
    assert out["n_iso"].tolist() == [3.0, 2.0, 3.0]
    assert out["usage_sd"].iloc[0] == pytest.approx(pd.Series([0.5, 0.7, 0.6, 0.8]).std())
    assert out["usage_sd"].iloc[1] == pytest.approx(pd.Series([1.0, 0.0, 0.5]).std())


def test_distinct_pairs_carry_the_joint_keys():
    from isograph_benchmark.real_data.switch_orthogonal_confirm import distinct_pairs
    scored = pd.DataFrame({
        "gene": ["G1", "G1"], "gene_name": ["g1", "g1"], "t1": ["A", "B"], "t2": ["B", "A"],
        "pair_detected": [True, True], "switch_like": [True, True],
        "usage_spearman": [-0.5, -0.5], "max_mean_if": [0.6, 0.6],
        "detect_rate": [0.75, 0.75], "n_iso": [3.0, 3.0], "usage_sd": [0.1, 0.2],
        "anchored_usable": [False, True],
    })
    dp = distinct_pairs(scored)
    assert len(dp) == 1 and dp["n_orientations"].iloc[0] == 2
    assert dp["usage_sd"].iloc[0] == pytest.approx(0.15)   # mean of both anchored SDs
    assert bool(dp["anchored_usable"].iloc[0])
    scored.loc[1, "detect_rate"] = 0.5
    with pytest.raises(ValueError):
        distinct_pairs(scored)
