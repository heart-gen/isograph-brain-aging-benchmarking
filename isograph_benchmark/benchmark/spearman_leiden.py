"""Spearman-rank + Leiden community detection baseline method."""

from __future__ import annotations

from dataclasses import dataclass, field

import networkx as nx
import numpy as np
import pandas as pd

from isograph.features.channels import gene_feature_channels, make_feature_scores
from isograph.models.base import FitArtifacts, NetworkModel


@dataclass
class SpearmanLeidenConfig:
    # Minimum Spearman r (positive correlation) to include an edge.
    min_r: float = 0.30
    leiden_resolution: float = 1.0
    min_module_size: int = 2
    trait_columns: list = field(default_factory=list)
    residualize_covariates: list = field(default_factory=list)


@dataclass
class SpearmanLeidenModel(NetworkModel):
    """Spearman correlation network on PSI switch coordinates with Leiden clustering.

    Uses the switch (PSI first-PC) feature channel — the same signal IsoGraph
    models — then applies plain Spearman r thresholding and Leiden community
    detection. This ablates IsoGraph's partial-correlation and significance-test
    components while keeping the same feature input, making it a direct measure
    of whether those components add value beyond simple correlation.
    """

    config: SpearmanLeidenConfig

    def fit(
        self,
        transcript_counts: np.ndarray,
        transcript_table: pd.DataFrame,
        sample_table: pd.DataFrame,
        gene_counts: np.ndarray | None = None,
        gene_table: pd.DataFrame | None = None,
    ) -> FitArtifacts:
        from scipy.stats import rankdata

        all_matrix, feature_info = gene_feature_channels(
            transcript_counts, transcript_table, gene_counts, gene_table
        )
        # Switch (PSI first-PC) features carry the isoform-switching signal.
        # Fall back to abundance if no switch features are present.
        switch_mask = feature_info["feature_type"] == "switch"
        if switch_mask.any():
            feat_matrix = all_matrix[switch_mask.values]
            feat_info = feature_info[switch_mask].reset_index(drop=True)
        else:
            feat_matrix = all_matrix
            feat_info = feature_info.reset_index(drop=True)

        gene_ids = feat_info["gene_id"].tolist()
        n_genes = len(gene_ids)

        # Spearman rank correlation: rank each gene's PSI values across samples, then Pearson
        ranks = np.apply_along_axis(rankdata, 1, feat_matrix)  # (n_genes, n_samples)
        corr = np.corrcoef(ranks) if n_genes > 1 else np.array([[1.0]])
        np.fill_diagonal(corr, 0.0)

        # Build network from positive correlation edges (r >= min_r)
        min_r = self.config.min_r
        i_idx, j_idx = np.where(np.triu(corr >= min_r, k=1))
        graph = nx.Graph()
        graph.add_nodes_from(gene_ids)
        edge_rows = []
        for i, j in zip(i_idx.tolist(), j_idx.tolist()):
            r = float(corr[i, j])
            graph.add_edge(gene_ids[i], gene_ids[j], weight=r)
            edge_rows.append({"gene_id_a": gene_ids[i], "gene_id_b": gene_ids[j], "weight": r})

        module_table = self._module_table(graph)
        feature_scores = make_feature_scores(feat_matrix, feat_info, sample_table)

        return FitArtifacts(
            module_table=module_table,
            edge_table=pd.DataFrame(edge_rows),
            trait_table=pd.DataFrame(columns=["module_id", "trait", "effect", "pvalue"]),
            feature_scores=feature_scores,
        )
