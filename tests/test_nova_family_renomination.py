from __future__ import annotations

import numpy as np
import pandas as pd

from isograph_benchmark.real_data.nova_family_renomination import (
    _candidate_identity_sha256,
    _cluster_offsets,
    _fit_regulons,
    _safe_exp,
    _unique_transcript_opportunity,
)


def test_cluster_offsets_require_three_instances_within_span() -> None:
    assert _cluster_offsets([0, 8, 20, 70], 4, 3, 30) == [(0, 2)]
    assert _cluster_offsets([0, 15, 29], 4, 3, 30) == []


def test_safe_exp_preserves_extreme_finite_estimates() -> None:
    assert np.isinf(_safe_exp(1_000.0))
    assert _safe_exp(-1_000.0) == 0.0
    assert _safe_exp(0.0) == 1.0


def test_candidate_identity_hash_is_order_invariant_and_membership_sensitive() -> None:
    frame = pd.DataFrame(
        [
            {
                "tree": "brainseq",
                "region": "hippocampus",
                "module_id": "M001",
                "gene": "GENE1",
                "transcript_id_1": "TX1",
                "transcript_id_2": "TX2",
                "cluster_positive_transcript": "TX1",
                "cluster_negative_transcript": "TX2",
            },
            {
                "tree": "gtex",
                "region": "hippocampus",
                "module_id": "M002",
                "gene": "GENE2",
                "transcript_id_1": "TX3",
                "transcript_id_2": "TX4",
                "cluster_positive_transcript": "TX4",
                "cluster_negative_transcript": "TX3",
            },
        ]
    )
    baseline = _candidate_identity_sha256(frame)
    assert _candidate_identity_sha256(frame.iloc[::-1]) == baseline
    changed = frame.copy()
    changed.loc[0, "gene"] = "GENE3"
    assert _candidate_identity_sha256(changed) != baseline


def test_gene_opportunity_sums_each_eligible_transcript_once() -> None:
    pairs = pd.DataFrame(
        {
            "transcript_id_1": ["TX1", "TX1"],
            "intronic_opportunity_nt_1": [100, 100],
            "transcript_id_2": ["TX2", "TX3"],
            "intronic_opportunity_nt_2": [200, 300],
        }
    )
    assert _unique_transcript_opportunity(pairs) == (3, 600.0)


def test_adjusted_regulon_model_preserves_tree_region_namespace() -> None:
    rows = []
    for tree in ["brainseq", "gtex"]:
        for index in range(20):
            rows.append(
                {
                    "tree": tree,
                    "region": "hippocampus",
                    "gene": f"{tree}_{index}",
                    "module_id": "M1" if index < 10 else "M2",
                    "go_invisible": True,
                    "pool_source": "module_genes",
                    "nova_family_switched": bool(index % 3 == 0),
                    "log_total_intronic_opportunity": 5.0 + index / 100,
                    "log_eligible_pair_count": 1.0,
                    "log_transcript_count": 1.5,
                }
            )
    stage = {
        "model_covariates": [
            "log_total_intronic_opportunity",
            "log_eligible_pair_count",
            "log_transcript_count",
        ],
        "min_module_genes": 3,
    }
    result = _fit_regulons(pd.DataFrame(rows), stage)
    assert len(result) == 4
    assert set(result["tree"]) == {"brainseq", "gtex"}
    assert set(result["region"]) == {"hippocampus"}
