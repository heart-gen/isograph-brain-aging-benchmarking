from __future__ import annotations

import math

import pandas as pd

from isograph_benchmark.real_data.neuronal_clip_motif_qc import _cluster_instances
from isograph_benchmark.real_data.neuronal_clip_overlap import (
    _interval_calls,
    _matched_effects,
)


def test_nova_cluster_rule_is_coordinate_preserving() -> None:
    instances = pd.DataFrame(
        {
            "instance_id": ["i1", "i2", "i3", "i4"],
            "transcript_id": ["tx1"] * 4,
            "gene_id": ["g1"] * 4,
            "chrom": ["chr1"] * 4,
            "start": [100, 108, 120, 170],
            "end": [104, 112, 124, 174],
            "strand": ["+"] * 4,
            "flank_id": ["f1"] * 4,
            "structural_class": ["donor_flank"] * 4,
            "intron_index": [1] * 4,
            "sense_offset": [0, 8, 20, 70],
            "motif": ["TCAT", "CCAC", "TCAC", "CCAT"],
        }
    )
    clusters = _cluster_instances(instances, min_instances=3, max_span=30)
    assert len(clusters) == 1
    cluster = clusters.iloc[0]
    assert cluster["start"] == 100
    assert cluster["end"] == 124
    assert cluster["n_instances"] == 3
    assert cluster["span_nt"] == 24
    assert cluster["evidence_scope"] == "exploratory_qc_only"


def test_interval_calls_require_strand_and_qvalue() -> None:
    windows = pd.DataFrame(
        {
            "window_id": ["w1", "w2"],
            "chrom": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [200, 300],
            "strand": ["+", "-"],
            "gene_id": ["ENSG1.1", "ENSG2.1"],
        }
    )
    tested = pd.DataFrame(
        {
            "chr": ["chr1", "chr1", "chr1"],
            "start": [150, 220, 220],
            "end": [250, 260, 260],
            "strand": ["+", "+", "-"],
            "input": [2, 3, 4],
            "clip": [10, 11, 12],
            "qvalue": [0.01, 0.01, 0.20],
        }
    )
    calls = _interval_calls(windows, tested, q_max=0.05, require_same_strand=True)
    first = calls.set_index("window_id").loc["w1"]
    second = calls.set_index("window_id").loc["w2"]
    assert bool(first["callable"])
    assert bool(first["bound"])
    assert bool(second["callable"])
    assert not bool(second["bound"])

    gene_callable = _interval_calls(
        windows,
        tested,
        q_max=0.05,
        require_same_strand=True,
        callable_gene_ids={"ENSG1"},
    ).set_index("window_id")
    assert bool(gene_callable.loc["w1", "callable"])
    assert not bool(gene_callable.loc["w2", "callable"])


def test_matched_effects_use_only_jointly_callable_sets() -> None:
    windows = pd.DataFrame(
        {
            "window_id": ["c1", "n1", "c2", "n2", "c3", "n3"],
            "candidate_id": ["a", "a", "b", "b", "d", "d"],
            "rbp": ["TARDBP"] * 6,
            "match_id": ["m1", "m1", "m2", "m2", "m3", "m3"],
            "window_role": ["case", "control"] * 3,
        }
    )
    consensus = pd.DataFrame(
        {
            "window_id": windows["window_id"],
            "context_callable": [True, True, True, True, True, False],
            "context_bound": [True, False, False, True, True, False],
        }
    )
    context = {
        "context_id": "ctx",
        "rbp": "TARDBP",
        "analysis_scope": "primary",
        "evidence_tier": 1,
        "qc_type": "enriched_windows",
    }
    result = _matched_effects(consensus, windows, context, min_discordant=1)
    assert result["callable_matched_sets"] == 2
    assert result["case_only_discordant"] == 1
    assert result["control_only_discordant"] == 1
    assert result["matched_odds_ratio"] == 1
    assert math.isclose(result["confirmatory_pvalue"], 0.75)
