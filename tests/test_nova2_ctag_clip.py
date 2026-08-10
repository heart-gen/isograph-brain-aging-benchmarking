from __future__ import annotations

import gzip

import pandas as pd

from isograph_benchmark.real_data.nova2_ctag_clip import (
    _bedgraph_calls,
    _cell_type_interaction,
    _effect_row,
    _map_chain,
)


def test_chain_mapping_preserves_plus_and_flips_minus(tmp_path) -> None:
    chain = tmp_path / "test.chain.gz"
    content = "\n".join(
        [
            "chain 1000 chr1 1000 + 100 300 chrM 1000 + 500 700 1",
            "200",
            "",
            "chain 900 chr2 1000 + 100 200 chrN 1000 - 300 400 2",
            "100",
            "",
        ]
    )
    with gzip.open(chain, "wt") as handle:
        handle.write(content)
    intervals = pd.DataFrame(
        {
            "interval_id": ["plus", "minus", "gap"],
            "chrom": ["chr1", "chr2", "chr1"],
            "start": [120, 120, 290],
            "end": [150, 150, 310],
            "strand": ["+", "+", "+"],
        }
    )
    mapped = _map_chain(chain, intervals).set_index("interval_id")
    assert mapped.loc["plus", "mapped_chrom"] == "chrM"
    assert mapped.loc["plus", "mapped_start"] == 520
    assert mapped.loc["plus", "mapped_end"] == 550
    assert mapped.loc["plus", "mapped_strand"] == "+"
    assert mapped.loc["minus", "mapped_start"] == 650
    assert mapped.loc["minus", "mapped_end"] == 680
    assert mapped.loc["minus", "mapped_strand"] == "-"
    assert mapped.loc["gap", "mapping_status"] == "unmapped_or_crosses_alignment_gap"


def test_bedgraph_calls_use_reciprocal_mapping_and_coverage(tmp_path) -> None:
    bedgraph = tmp_path / "signal.bedgraph.gz"
    with gzip.open(bedgraph, "wt") as handle:
        handle.write('track type=bedGraph name="test"\n')
        handle.write("chrM\t110\t130\t2\n")
        handle.write("chrM\t150\t160\t1\n")
    windows = pd.DataFrame(
        {
            "window_id": ["bound", "unbound", "unmapped"],
            "start": [0, 0, 0],
            "end": [20, 20, 20],
            "reciprocal_mapped": [True, True, False],
            "forward_mapped_chrom": ["chrM", "chrM", pd.NA],
            "forward_mapped_start": [100, 200, pd.NA],
            "forward_mapped_end": [120, 220, pd.NA],
        }
    )
    calls = _bedgraph_calls(windows, bedgraph, "ctx", 1).set_index("window_id")
    assert bool(calls.loc["bound", "callable"])
    assert bool(calls.loc["bound", "bound"])
    assert calls.loc["bound", "covered_nt"] == 10
    assert not bool(calls.loc["unbound", "bound"])
    assert not bool(calls.loc["unmapped", "callable"])


def test_effect_row_uses_gene_level_discordance() -> None:
    frame = pd.DataFrame(
        {
            "case_callable": [True] * 5,
            "control_callable": [True] * 5,
            "case_bound": [True, True, True, False, True],
            "control_bound": [False, False, True, True, False],
        }
    )
    effect = _effect_row(frame, "gene")
    assert effect["callable_units"] == 5
    assert effect["case_only_discordant"] == 3
    assert effect["control_only_discordant"] == 1
    assert effect["informative_discordant"] == 4
    assert effect["matched_odds_ratio_haldane"] > 1


def test_cell_type_interaction_excludes_non_callable_genes() -> None:
    rows = []
    for context_id in ["emx1", "gad2"]:
        rows.extend(
            [
                {
                    "context_id": context_id,
                    "tree": "brainseq",
                    "region": "hippocampus",
                    "module_id": "M001",
                    "gene": "callable",
                    "window_width": 100,
                    "case_callable": True,
                    "control_callable": True,
                    "case_bound": context_id == "emx1",
                    "control_bound": False,
                },
                {
                    "context_id": context_id,
                    "tree": "brainseq",
                    "region": "hippocampus",
                    "module_id": "M001",
                    "gene": "not_callable",
                    "window_width": 100,
                    "case_callable": False,
                    "control_callable": False,
                    "case_bound": False,
                    "control_bound": False,
                },
            ]
        )
    stage = {"primary_width": 100, "primary_contexts": ["emx1", "gad2"]}
    result = _cell_type_interaction(pd.DataFrame(rows), stage).iloc[0]
    assert result["jointly_callable_genes"] == 1
    assert result["positive_difference"] == 1
    assert result["negative_difference"] == 0
