from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from isograph.explain.annotation import load_annotation_table
from isograph.explain.config import ExplainConfig
from isograph.explain.core import explain_module
from isograph.explain.structure import annotate_switch_pairs
from isograph.io.artifacts import DatasetBundle, load_dataset_bundle
from isograph_benchmark.paths import ensure_dir


def artifact_sample_columns(artifact_dir: Path | str) -> list[str]:
    """Return the sample columns present in an IsoGraph feature score table."""
    feature_scores = pd.read_parquet(Path(artifact_dir) / "feature_scores.parquet")
    return [str(c) for c in feature_scores.columns if c != "gene_id"]


def _sample_index_for_artifact(bundle: DatasetBundle, artifact_dir: Path | str) -> list[str]:
    bundle_sample_ids = bundle.sample_table["sample_id"].astype(str).tolist()
    score_columns = artifact_sample_columns(artifact_dir)
    if len(score_columns) == len(bundle_sample_ids):
        return score_columns
    return bundle_sample_ids


def build_transcript_usage_inputs(
    bundle: DatasetBundle,
    artifact_dir: Path | str,
    pseudocount: float = 0.5,
    restrict_genes: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build samples-by-transcript usage features and metadata for explain-module."""
    transcript_table = bundle.feature_tables["transcript"].reset_index(drop=True).copy()
    transcript_counts_full = np.asarray(bundle.matrices["transcript_counts"], dtype=np.float32)
    if transcript_counts_full.shape[0] != len(transcript_table):
        raise ValueError(
            "transcript_counts row count does not match transcript metadata: "
            f"{transcript_counts_full.shape[0]} != {len(transcript_table)}"
        )
    transcript_table["_matrix_idx"] = np.arange(len(transcript_table))
    if restrict_genes:
        transcript_table = transcript_table.loc[
            transcript_table["gene_id"].astype(str).isin(restrict_genes)
        ].copy()
    matrix_idx = transcript_table["_matrix_idx"].to_numpy(dtype=int)
    transcript_counts = transcript_counts_full[matrix_idx]
    transcript_table = transcript_table.drop(columns=["_matrix_idx"]).reset_index(drop=True)

    usage = np.empty_like(transcript_counts, dtype=np.float32)
    for _, frame in transcript_table.groupby("gene_id", sort=False):
        idx = frame.index.to_numpy()
        counts = transcript_counts[idx].astype(np.float32, copy=False) + np.float32(pseudocount)
        totals = counts.sum(axis=0, keepdims=True)
        usage[idx] = counts / np.clip(totals, 1e-12, None)

    sample_index = _sample_index_for_artifact(bundle, artifact_dir)
    if len(sample_index) != usage.shape[1]:
        raise ValueError(
            "artifact/sample count does not match transcript_counts columns: "
            f"{len(sample_index)} != {usage.shape[1]}"
        )

    transcript_ids = transcript_table["transcript_id"].astype(str).tolist()
    feature_table = pd.DataFrame(usage.T, index=pd.Index(sample_index, name="sample_id"), columns=transcript_ids)

    feature_meta = pd.DataFrame(
        {
            "feature_id": transcript_ids,
            "gene_id": transcript_table["gene_id"].astype(str).to_numpy(),
            "feature_type": "transcript_usage",
            "transcript_id": transcript_ids,
        }
    )
    for col in ["gene_name", "transcript_name", "transcript_type", "Length", "EffectiveLength", "length"]:
        if col in transcript_table.columns and col not in feature_meta.columns:
            feature_meta[col] = transcript_table[col].to_numpy()
    return feature_table, feature_meta


def _module_genes(artifact_dir: Path, module_ids: list[str] | None) -> set[str]:
    modules = pd.read_parquet(artifact_dir / "modules.parquet")
    if modules.empty or "module_id" not in modules.columns:
        return set()
    if module_ids is not None:
        modules = modules.loc[modules["module_id"].astype(str).isin(set(map(str, module_ids)))]
    return set(modules["gene_id"].astype(str))


def build_switch_pairs(
    feature_meta: pd.DataFrame,
    module_genes: set[str],
    feature_table: pd.DataFrame | None = None,
    max_transcripts_per_gene: int = 6,
    max_pairs_per_gene: int = 15,
) -> pd.DataFrame:
    """Build bounded transcript pairs for GTF-based structural annotation."""
    if not module_genes:
        return pd.DataFrame(columns=["gene_id", "transcript_id_1", "transcript_id_2"])

    meta = feature_meta.loc[feature_meta["gene_id"].astype(str).isin(module_genes)].copy()
    if meta.empty:
        return pd.DataFrame(columns=["gene_id", "transcript_id_1", "transcript_id_2"])

    variance: pd.Series | None = None
    if feature_table is not None:
        available = [c for c in meta["feature_id"].astype(str) if c in feature_table.columns]
        if available:
            variance = feature_table[available].var(axis=0)

    rows: list[dict[str, str]] = []
    for gene_id, frame in meta.groupby("gene_id", sort=False):
        frame = frame.drop_duplicates("transcript_id").copy()
        if len(frame) < 2:
            continue
        if variance is not None:
            frame["_variance"] = frame["feature_id"].astype(str).map(variance).fillna(0.0)
            frame = frame.sort_values("_variance", ascending=False)
        transcripts = frame["transcript_id"].astype(str).head(max_transcripts_per_gene).tolist()
        n_pairs = 0
        for i, tx1 in enumerate(transcripts):
            for tx2 in transcripts[i + 1 :]:
                rows.append(
                    {
                        "gene_id": str(gene_id),
                        "transcript_id_1": tx1,
                        "transcript_id_2": tx2,
                    }
                )
                n_pairs += 1
                if n_pairs >= max_pairs_per_gene:
                    break
            if n_pairs >= max_pairs_per_gene:
                break
    return pd.DataFrame(rows, columns=["gene_id", "transcript_id_1", "transcript_id_2"])


def build_gtf_annotation_table(
    artifact_dir: Path | str,
    feature_table: pd.DataFrame,
    feature_meta: pd.DataFrame,
    module_ids: list[str] | None,
    gtf_path: Path | str,
    gtf_cache: Path | str | None = None,
    max_transcripts_per_gene: int = 6,
    max_pairs_per_gene: int = 15,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build structural transcript annotations for selected modules from a GTF."""
    artifact_dir = Path(artifact_dir)
    module_genes = _module_genes(artifact_dir, module_ids)
    switch_pairs = build_switch_pairs(
        feature_meta,
        module_genes=module_genes,
        feature_table=feature_table,
        max_transcripts_per_gene=max_transcripts_per_gene,
        max_pairs_per_gene=max_pairs_per_gene,
    )
    if switch_pairs.empty:
        annotations = pd.DataFrame(columns=["transcript_id", "gene_id"])
    else:
        if gtf_cache is not None:
            ensure_dir(Path(gtf_cache).parent)
        annotations = annotate_switch_pairs(switch_pairs, gtf_path, gtf_cache=gtf_cache)
    return annotations, switch_pairs


def explain_artifact_modules(
    artifact_dir: Path | str,
    bundle_path: Path | str,
    output_dir: Path | str,
    module_ids: list[str] | None = None,
    annotation_table: Path | str | pd.DataFrame | None = None,
    gtf_path: Path | str | None = None,
    gtf_cache: Path | str | None = None,
    max_transcripts_per_gene: int = 6,
    max_pairs_per_gene: int = 15,
    split_percentile: float = 50.0,
    min_complete_pairs: int = 3,
    fdr_method: str = "bh",
    plot: bool = False,
    output_format: str | list[str] = "png",
    pseudocount: float = 0.5,
) -> dict[str, Any]:
    """Run IsoGraph module interpretation for a fitted artifact directory."""
    artifact_dir = Path(artifact_dir)
    bundle = load_dataset_bundle(Path(bundle_path))
    selected_genes = _module_genes(artifact_dir, module_ids)
    feature_table, feature_meta = build_transcript_usage_inputs(
        bundle,
        artifact_dir=artifact_dir,
        pseudocount=pseudocount,
        restrict_genes=selected_genes or None,
    )
    annotations = load_annotation_table(annotation_table) if annotation_table is not None else None
    output_dir = ensure_dir(Path(output_dir))
    if annotations is None and gtf_path is not None:
        annotations, switch_pairs = build_gtf_annotation_table(
            artifact_dir=artifact_dir,
            feature_table=feature_table,
            feature_meta=feature_meta,
            module_ids=module_ids,
            gtf_path=gtf_path,
            gtf_cache=gtf_cache,
            max_transcripts_per_gene=max_transcripts_per_gene,
            max_pairs_per_gene=max_pairs_per_gene,
        )
        switch_pairs.to_parquet(output_dir / "structure_switch_pairs.parquet", index=False, compression="zstd")
        annotations.to_parquet(output_dir / "structure_annotations.parquet", index=False, compression="zstd")
    config = ExplainConfig(
        split_percentile=split_percentile,
        min_complete_pairs=min_complete_pairs,
        fdr_method=fdr_method,
        plot=plot,
        output_format=output_format,
    )
    return explain_module(
        artifact_dir=artifact_dir,
        feature_table=feature_table,
        feature_meta=feature_meta,
        module_ids=module_ids,
        output_dir=output_dir,
        config=config,
        annotation_table=annotations,
    )


def summarize_explain_results(
    results: dict[str, Any],
    dataset: str,
    collection: str | None = None,
) -> pd.DataFrame:
    """Compact per-module summary for real-data interpretation outputs."""
    rows: list[dict[str, Any]] = []
    for module_id, result in sorted(results.items()):
        gene_drivers = result.gene_driver_table
        transcript_polarity = result.transcript_polarity_table
        top_gene = None
        top_gene_abs_r = np.nan
        if not gene_drivers.empty and "r" in gene_drivers.columns:
            top = gene_drivers.iloc[gene_drivers["r"].abs().fillna(0).to_numpy().argmax()]
            top_gene = top.get("gene_id")
            top_gene_abs_r = abs(float(top.get("r", np.nan)))
        max_switch_strength = np.nan
        if not transcript_polarity.empty and "switch_strength" in transcript_polarity.columns:
            max_switch_strength = float(pd.to_numeric(transcript_polarity["switch_strength"], errors="coerce").max())
        rows.append(
            {
                "collection": collection,
                "dataset": dataset,
                "module_id": module_id,
                "n_module_genes": int(result.n_module_genes),
                "n_gene_drivers": int(len(gene_drivers)),
                "n_transcript_features": int(len(transcript_polarity)),
                "top_gene_id": top_gene,
                "top_gene_abs_r": top_gene_abs_r,
                "max_switch_strength": max_switch_strength,
            }
        )
    return pd.DataFrame(rows)
