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
