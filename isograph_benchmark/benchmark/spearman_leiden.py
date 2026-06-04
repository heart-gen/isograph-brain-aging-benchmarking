"""Spearman-rank + Leiden community detection baseline method."""

from __future__ import annotations

from dataclasses import dataclass, field

import networkx as nx
import numpy as np
import pandas as pd

from isograph.features.channels import gene_feature_channels, make_feature_scores
from isograph.features.residualize import build_design_matrix, residualize_rows
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
    """Modern comparator: Spearman correlation network + Leiden clustering.

    Consumes the identical ``gene_feature_channels`` input used by WGCNA and the
    IsoGraph models (abundance + switch feature rows, optionally residualized for
    covariates), builds a feature-feature graph from positive Spearman r >= min_r,
    runs Leiden community detection, then maps feature communities back to genes
    exactly as the WGCNA backend does.

    Holding the input features and the feature->gene mapping identical to WGCNA,
    this isolates the clustering choice: Spearman r + Leiden vs. WGCNA's signed
    topological-overlap matrix + dynamic tree cut.
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

        feat_matrix, feature_info = gene_feature_channels(
            transcript_counts, transcript_table, gene_counts, gene_table
        )
        if feat_matrix.size:
            design = build_design_matrix(sample_table, self.config.residualize_covariates)
            feat_matrix = residualize_rows(feat_matrix, design)

        feature_ids = feature_info["feature_id"].tolist()
        feature_to_gene = feature_info.set_index("feature_id")["gene_id"].to_dict()
        n_features = len(feature_ids)

        # Spearman rank correlation: rank each feature's values across samples, then Pearson
        if n_features > 1:
            ranks = np.apply_along_axis(rankdata, 1, feat_matrix)
            corr = np.corrcoef(ranks)
        else:
            corr = np.array([[1.0]])
        np.fill_diagonal(corr, 0.0)

        # Feature-feature graph from positive correlation edges (r >= min_r)
        min_r = self.config.min_r
        i_idx, j_idx = np.where(np.triu(corr >= min_r, k=1))
        graph = nx.Graph()
        graph.add_nodes_from(feature_ids)
        for i, j in zip(i_idx.tolist(), j_idx.tolist()):
            graph.add_edge(feature_ids[i], feature_ids[j], weight=float(corr[i, j]))

        # Leiden community detection on the feature graph, then map features -> genes.
        communities = self._detect_communities(graph)
        rows = []
        for module_index, feature_nodes in enumerate(communities):
            genes = sorted({feature_to_gene[f] for f in feature_nodes if f in feature_to_gene})
            if len(genes) < self.config.min_module_size:
                continue
            for gene_id in genes:
                rows.append({"gene_id": gene_id, "module_id": f"M{module_index:03d}"})
        module_table = pd.DataFrame(rows, columns=["gene_id", "module_id"])

        # Gene-level edge table (map feature edges to gene pairs, drop self-loops)
        edge_rows = []
        for i, j in zip(i_idx.tolist(), j_idx.tolist()):
            ga = feature_to_gene.get(feature_ids[i])
            gb = feature_to_gene.get(feature_ids[j])
            if ga is not None and gb is not None and ga != gb:
                edge_rows.append({"source": ga, "target": gb, "weight": float(corr[i, j])})
        edge_table = pd.DataFrame(edge_rows, columns=["source", "target", "weight"])

        feature_scores = make_feature_scores(feat_matrix, feature_info, sample_table)

        return FitArtifacts(
            module_table=module_table,
            edge_table=edge_table,
            trait_table=pd.DataFrame(columns=["module_id", "trait", "effect", "pvalue"]),
            feature_scores=feature_scores,
        )
