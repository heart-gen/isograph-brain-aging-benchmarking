"""Rebuild NOVA-family switch regulons with explicit motif-opportunity control."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from statsmodels.tools.sm_exceptions import PerfectSeparationError

from isograph_benchmark.config import load_yaml
from isograph_benchmark.paths import ensure_dir, region_store
from isograph_benchmark.real_data.neuronal_clip_fetch import sha256
from isograph_benchmark.real_data.neuronal_clip_motif_qc import (
    _sense_sequence,
    _unique_flanks,
)
from isograph_benchmark.real_data.neuronal_clip_validation import (
    DEFAULT_CONFIG,
    IndexedFasta,
    TranscriptModel,
    _configured_path,
    _load_transcript_models,
    _stable_id,
)
from isograph_benchmark.real_data.qtl_anchoring import _bare
from isograph_benchmark.real_data.rbp_regulon import _REGIONS, _gene_tags


CANDIDATE_IDENTITY_COLUMNS = [
    "tree",
    "region",
    "module_id",
    "gene",
    "transcript_id_1",
    "transcript_id_2",
    "cluster_positive_transcript",
    "cluster_negative_transcript",
]


def _safe_exp(value: float) -> float:
    """Exponentiate a finite log estimate without aborting report generation."""
    if value > math.log(np.finfo(float).max):
        return math.inf
    if value < math.log(np.finfo(float).tiny):
        return 0.0
    return math.exp(value)


def _candidate_identity_sha256(frame: pd.DataFrame) -> str:
    view = frame[CANDIDATE_IDENTITY_COLUMNS].astype("string").fillna("")
    view = view.sort_values(CANDIDATE_IDENTITY_COLUMNS, kind="mergesort")
    digest = hashlib.sha256()
    for row in view.itertuples(index=False, name=None):
        digest.update(
            (
                json.dumps(list(row), ensure_ascii=True, separators=(",", ":")) + "\n"
            ).encode()
        )
    return digest.hexdigest()


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(f".{path.name}.tmp")
    frame.to_parquet(temporary, index=False, compression="zstd")
    os.replace(temporary, path)


def _atomic_text(text: str, path: Path) -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text)
    os.replace(temporary, path)


def _load_pairs() -> tuple[pd.DataFrame, dict[str, str]]:
    frames: list[pd.DataFrame] = []
    hashes: dict[str, str] = {}
    for tree, region in _REGIONS:
        path = region_store(tree, region, "isograph_vae",
            "module_interpret",
            "structure_switch_pairs.parquet",
        )
        if not path.is_file():
            continue
        frame = pd.read_parquet(
            path, columns=["gene_id", "transcript_id_1", "transcript_id_2"]
        )
        frame["tree"] = tree
        frame["region"] = region
        frame["gene"] = _bare(frame["gene_id"])
        frame["transcript_a"] = frame[["transcript_id_1", "transcript_id_2"]].min(
            axis=1
        )
        frame["transcript_b"] = frame[["transcript_id_1", "transcript_id_2"]].max(
            axis=1
        )
        frame = frame.drop_duplicates(
            ["region", "gene", "transcript_a", "transcript_b"]
        )
        frames.append(frame)
        hashes[f"{tree}/{region}"] = sha256(path)
    if not frames:
        raise FileNotFoundError("No structure_switch_pairs.parquet inputs found")
    return pd.concat(frames, ignore_index=True), hashes


def _cluster_offsets(
    offsets: list[int], motif_length: int, min_instances: int, max_span: int
) -> list[tuple[int, int]]:
    clusters: list[tuple[int, int]] = []
    index = 0
    while index < len(offsets):
        final = index
        while final + 1 < len(offsets):
            span = offsets[final + 1] + motif_length - offsets[index]
            if span > max_span:
                break
            final += 1
        if final - index + 1 < min_instances:
            index += 1
            continue
        clusters.append((index, final))
        index = final + 1
    return clusters


def _motif_coordinates(
    model: TranscriptModel,
    flank_start: int,
    flank_end: int,
    offset: int,
    motif_length: int,
) -> tuple[int, int]:
    if model.strand == "+":
        start = flank_start + offset
        return start, start + motif_length
    return flank_end - offset - motif_length, flank_end - offset


def _score_transcript(
    model: TranscriptModel,
    fasta: IndexedFasta,
    width: int,
    motifs: set[str],
    min_instances: int,
    max_span: int,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    motif_lengths = {len(motif) for motif in motifs}
    if len(motif_lengths) != 1:
        raise ValueError("NOVA-family motifs must have one common length")
    motif_length = motif_lengths.pop()
    opportunity = 0
    total_instances = 0
    cluster_rows: list[dict[str, Any]] = []
    for flank in _unique_flanks(model, width):
        sequence = _sense_sequence(
            fasta.fetch(flank.chrom, flank.start, flank.end), model.strand
        )
        opportunity += len(sequence)
        offsets = [
            offset
            for offset in range(max(0, len(sequence) - motif_length + 1))
            if sequence[offset : offset + motif_length] in motifs
        ]
        total_instances += len(offsets)
        for first, final in _cluster_offsets(
            offsets, motif_length, min_instances, max_span
        ):
            cluster_offsets = offsets[first : final + 1]
            coordinates = [
                _motif_coordinates(model, flank.start, flank.end, offset, motif_length)
                for offset in cluster_offsets
            ]
            start = min(item[0] for item in coordinates)
            end = max(item[1] for item in coordinates)
            sense_span = cluster_offsets[-1] + motif_length - cluster_offsets[0]
            cluster_rows.append(
                {
                    "cluster_id": _stable_id(
                        "NRC",
                        model.transcript_id,
                        flank.flank_id,
                        start,
                        end,
                        len(cluster_offsets),
                    ),
                    "transcript_id": model.transcript_id,
                    "gene_id": model.gene_id,
                    "chrom": flank.chrom,
                    "start": start,
                    "end": end,
                    "strand": model.strand,
                    "coordinate_system": "0-based_half-open",
                    "flank_id": flank.flank_id,
                    "structural_class": flank.structural_class,
                    "intron_index": flank.intron_index,
                    "n_instances": len(cluster_offsets),
                    "span_nt": sense_span,
                    "cluster_rule": f">={min_instances}_YCAY_within_{max_span}nt",
                }
            )
    summary = {
        "transcript_id": model.transcript_id,
        "gene_id": model.gene_id,
        "exon_count": len(model.exons),
        "intronic_opportunity_nt": opportunity,
        "nova_family_instance_count": total_instances,
        "nova_family_cluster_count": len(cluster_rows),
        "max_cluster_instance_count": max(
            (int(row["n_instances"]) for row in cluster_rows), default=0
        ),
        "cluster_density_per_kb": (
            1000.0 * len(cluster_rows) / opportunity if opportunity else 0.0
        ),
        "model_status": "scored",
    }
    return summary, cluster_rows


def _score_transcripts(
    transcript_ids: set[str],
    models: dict[str, TranscriptModel],
    fasta_path: Path,
    stage: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    motifs = {str(motif).upper() for motif in stage["motifs_dna_sense"]}
    summaries: list[dict[str, Any]] = []
    clusters: list[dict[str, Any]] = []
    with IndexedFasta(fasta_path) as fasta:
        for index, transcript_id in enumerate(sorted(transcript_ids), start=1):
            model = models.get(transcript_id)
            if model is None:
                summaries.append(
                    {
                        "transcript_id": transcript_id,
                        "gene_id": "",
                        "exon_count": 0,
                        "intronic_opportunity_nt": 0,
                        "nova_family_instance_count": 0,
                        "nova_family_cluster_count": 0,
                        "max_cluster_instance_count": 0,
                        "cluster_density_per_kb": 0.0,
                        "model_status": "missing_gtf_model",
                    }
                )
                continue
            summary, rows = _score_transcript(
                model,
                fasta,
                int(stage["primary_width"]),
                motifs,
                int(stage["min_instances_per_cluster"]),
                int(stage["max_cluster_span_nt"]),
            )
            summaries.append(summary)
            clusters.extend(rows)
            if index % 5000 == 0:
                print(
                    f"scored {index:,}/{len(transcript_ids):,} transcripts", flush=True
                )
    return pd.DataFrame(summaries), pd.DataFrame(clusters)


def _pair_scores(pairs: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    score_columns = [
        "transcript_id",
        "exon_count",
        "intronic_opportunity_nt",
        "nova_family_instance_count",
        "nova_family_cluster_count",
        "max_cluster_instance_count",
        "cluster_density_per_kb",
        "model_status",
    ]
    first = scores[score_columns].rename(
        columns={
            column: f"{column}_1"
            for column in score_columns
            if column != "transcript_id"
        }
        | {"transcript_id": "transcript_id_1"}
    )
    second = scores[score_columns].rename(
        columns={
            column: f"{column}_2"
            for column in score_columns
            if column != "transcript_id"
        }
        | {"transcript_id": "transcript_id_2"}
    )
    scored = pairs.merge(
        first, on="transcript_id_1", how="left", validate="many_to_one"
    ).merge(second, on="transcript_id_2", how="left", validate="many_to_one")
    for column in [
        "intronic_opportunity_nt_1",
        "intronic_opportunity_nt_2",
        "nova_family_instance_count_1",
        "nova_family_instance_count_2",
        "nova_family_cluster_count_1",
        "nova_family_cluster_count_2",
    ]:
        scored[column] = scored[column].fillna(0).astype(int)
    return scored


def _unique_transcript_opportunity(group: pd.DataFrame) -> tuple[int, float]:
    first = group[["transcript_id_1", "intronic_opportunity_nt_1"]].rename(
        columns={
            "transcript_id_1": "transcript_id",
            "intronic_opportunity_nt_1": "intronic_opportunity_nt",
        }
    )
    second = group[["transcript_id_2", "intronic_opportunity_nt_2"]].rename(
        columns={
            "transcript_id_2": "transcript_id",
            "intronic_opportunity_nt_2": "intronic_opportunity_nt",
        }
    )
    opportunities = pd.concat([first, second], ignore_index=True)
    inconsistent = opportunities.groupby("transcript_id")[
        "intronic_opportunity_nt"
    ].nunique()
    if inconsistent.gt(1).any():
        raise ValueError("A transcript has inconsistent intronic-opportunity scores")
    opportunities = opportunities.drop_duplicates("transcript_id")
    return len(opportunities), float(opportunities["intronic_opportunity_nt"].sum())


def _tagged_pair_calls(
    all_pairs: pd.DataFrame, stage: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pair_frames: list[pd.DataFrame] = []
    gene_frames: list[pd.DataFrame] = []
    minimum_opportunity = int(stage["min_intronic_opportunity_nt_per_transcript"])
    for tree, region in _REGIONS:
        region_pairs = all_pairs[
            all_pairs["tree"].eq(tree) & all_pairs["region"].eq(region)
        ].copy()
        if region_pairs.empty:
            continue
        artifact = region_store(tree, region, "isograph_vae")
        tags, pool_source = _gene_tags(artifact, float(stage["discovery_fdr"]))
        if tags.empty:
            continue
        tags = tags.groupby("gene", as_index=False).agg(
            module_id=("module_id", "first"), go_invisible=("go_invisible", "max")
        )
        region_pairs = region_pairs.merge(
            tags, on="gene", how="inner", validate="many_to_one"
        )
        region_pairs["pool_source"] = pool_source
        region_pairs["eligible_pair"] = region_pairs["intronic_opportunity_nt_1"].ge(
            minimum_opportunity
        ) & region_pairs["intronic_opportunity_nt_2"].ge(minimum_opportunity)
        positive_1 = region_pairs["nova_family_cluster_count_1"].gt(0)
        positive_2 = region_pairs["nova_family_cluster_count_2"].gt(0)
        region_pairs["nova_family_switched_pair"] = region_pairs[
            "eligible_pair"
        ] & positive_1.ne(positive_2)
        region_pairs["cluster_positive_transcript"] = np.where(
            region_pairs["nova_family_switched_pair"] & positive_1,
            region_pairs["transcript_id_1"],
            np.where(
                region_pairs["nova_family_switched_pair"] & positive_2,
                region_pairs["transcript_id_2"],
                "",
            ),
        )
        region_pairs["cluster_negative_transcript"] = np.where(
            region_pairs["nova_family_switched_pair"] & positive_1,
            region_pairs["transcript_id_2"],
            np.where(
                region_pairs["nova_family_switched_pair"] & positive_2,
                region_pairs["transcript_id_1"],
                "",
            ),
        )
        pair_frames.append(region_pairs)
        eligible = region_pairs[region_pairs["eligible_pair"]].copy()
        if eligible.empty:
            continue
        gene_rows: list[dict[str, Any]] = []
        for gene, group in eligible.groupby("gene", sort=False):
            transcript_count, total_opportunity = _unique_transcript_opportunity(group)
            first = group.iloc[0]
            gene_rows.append(
                {
                    "region": region,
                    "tree": tree,
                    "gene": gene,
                    "module_id": first["module_id"],
                    "go_invisible": bool(first["go_invisible"]),
                    "pool_source": pool_source,
                    "nova_family_switched": bool(
                        group["nova_family_switched_pair"].any()
                    ),
                    "eligible_pair_count": int(len(group)),
                    "transcript_count": transcript_count,
                    "total_intronic_opportunity": total_opportunity,
                }
            )
        gene_frame = pd.DataFrame(gene_rows)
        gene_frame["log_total_intronic_opportunity"] = np.log1p(
            gene_frame["total_intronic_opportunity"]
        )
        gene_frame["log_eligible_pair_count"] = np.log1p(
            gene_frame["eligible_pair_count"]
        )
        gene_frame["log_transcript_count"] = np.log1p(gene_frame["transcript_count"])
        gene_frames.append(gene_frame)
    pair_calls = (
        pd.concat(pair_frames, ignore_index=True) if pair_frames else pd.DataFrame()
    )
    gene_calls = (
        pd.concat(gene_frames, ignore_index=True) if gene_frames else pd.DataFrame()
    )
    return pair_calls, gene_calls


def _fit_regulons(gene_calls: pd.DataFrame, stage: dict[str, Any]) -> pd.DataFrame:
    covariates = [str(value) for value in stage["model_covariates"]]
    rows: list[dict[str, Any]] = []
    for (tree, region), region_genes in gene_calls.groupby(
        ["tree", "region"], sort=True
    ):
        for module_id, module_genes in region_genes.groupby("module_id", sort=True):
            in_module = region_genes["module_id"].eq(module_id).astype(float)
            base = {
                "tree": tree,
                "region": region,
                "module_id": module_id,
                "rbp": "NOVA_FAMILY",
                "go_invisible": bool(module_genes["go_invisible"].mode().iloc[0]),
                "pool_source": str(module_genes["pool_source"].iloc[0]),
                "universe_genes": int(len(region_genes)),
                "module_genes": int(len(module_genes)),
                "universe_switched": int(region_genes["nova_family_switched"].sum()),
                "module_switched": int(module_genes["nova_family_switched"].sum()),
            }
            if len(module_genes) < int(stage["min_module_genes"]):
                rows.append({**base, "model_status": "module_too_small"})
                continue
            outcome = region_genes["nova_family_switched"].astype(int)
            if outcome.nunique() < 2 or in_module.nunique() < 2:
                rows.append({**base, "model_status": "outcome_or_predictor_constant"})
                continue
            design = pd.DataFrame({"intercept": 1.0, "in_module": in_module})
            for covariate in covariates:
                design[covariate] = pd.to_numeric(
                    region_genes[covariate], errors="raise"
                )
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fit = sm.GLM(
                        outcome,
                        design,
                        family=sm.families.Binomial(),
                    ).fit(maxiter=200)
                coefficient = float(fit.params["in_module"])
                standard_error = float(fit.bse["in_module"])
                pvalue = float(fit.pvalues["in_module"])
                ci = fit.conf_int().loc["in_module"]
                if not all(
                    math.isfinite(value)
                    for value in [coefficient, standard_error, pvalue, *ci]
                ):
                    raise ValueError("non-finite logistic estimate")
                rows.append(
                    {
                        **base,
                        "coefficient": coefficient,
                        "standard_error": standard_error,
                        "odds_ratio": _safe_exp(coefficient),
                        "ci_low": _safe_exp(float(ci.iloc[0])),
                        "ci_high": _safe_exp(float(ci.iloc[1])),
                        "pvalue": pvalue,
                        "model_status": "fit",
                    }
                )
            except (PerfectSeparationError, ValueError, np.linalg.LinAlgError) as error:
                rows.append(
                    {**base, "model_status": f"fit_failed:{type(error).__name__}"}
                )
    results = pd.DataFrame(rows)
    results["qvalue"] = np.nan
    fitted = results["model_status"].eq("fit")
    if fitted.any():
        results.loc[fitted, "qvalue"] = multipletests(
            results.loc[fitted, "pvalue"], method="fdr_bh"
        )[1]
    return results


def _freeze_candidates(
    pair_calls: pd.DataFrame,
    regulons: pd.DataFrame,
    stage: dict[str, Any],
    pair_hashes: dict[str, str],
    config_path: Path,
) -> pd.DataFrame:
    nominations = regulons[
        regulons["model_status"].eq("fit")
        & regulons["qvalue"].le(float(stage["regulon_q_max"]))
        & regulons["coefficient"].gt(0)
    ].copy()
    candidates = pair_calls[pair_calls["nova_family_switched_pair"]].merge(
        nominations[
            [
                "tree",
                "region",
                "module_id",
                "go_invisible",
                "odds_ratio",
                "pvalue",
                "qvalue",
                "module_genes",
                "module_switched",
                "pool_source",
            ]
        ],
        on=["tree", "region", "module_id"],
        how="inner",
        suffixes=("", "_regulon"),
        validate="many_to_one",
    )
    if candidates.empty:
        return candidates
    candidates["rbp"] = "NOVA_FAMILY"
    candidates["motif_scope"] = "intronic_ycay_cluster"
    candidates["motif_count_1"] = candidates["nova_family_cluster_count_1"]
    candidates["motif_count_2"] = candidates["nova_family_cluster_count_2"]
    candidates["switch_pairs_sha256"] = [
        pair_hashes[f"{tree}/{region}"]
        for tree, region in zip(candidates["tree"], candidates["region"])
    ]
    candidates["config_sha256"] = sha256(config_path)
    candidates["frozen_date"] = str(load_yaml(config_path)["frozen_date"])
    keep = [
        "tree",
        "region",
        "module_id",
        "go_invisible",
        "rbp",
        "gene",
        "gene_id",
        "transcript_id_1",
        "transcript_id_2",
        "cluster_positive_transcript",
        "cluster_negative_transcript",
        "motif_count_1",
        "motif_count_2",
        "intronic_opportunity_nt_1",
        "intronic_opportunity_nt_2",
        "nova_family_instance_count_1",
        "nova_family_instance_count_2",
        "odds_ratio",
        "pvalue",
        "qvalue",
        "module_genes",
        "module_switched",
        "pool_source",
        "motif_scope",
        "switch_pairs_sha256",
        "config_sha256",
        "frozen_date",
    ]
    return candidates[keep].sort_values(
        ["tree", "region", "module_id", "gene", "transcript_id_1", "transcript_id_2"]
    )


def _guard_expected(stage: dict[str, Any], observed: dict[str, int]) -> None:
    for key, value in observed.items():
        configured = stage.get(f"expected_{key}")
        if configured is not None and int(configured) != int(value):
            raise RuntimeError(
                f"NOVA-family freeze changed for {key}: {value} observed, {configured} expected"
            )


def _guard_candidate_identity(stage: dict[str, Any], observed: str) -> None:
    configured = stage.get("expected_candidate_identity_sha256")
    if configured is not None and str(configured) != observed:
        raise RuntimeError(
            "NOVA-family candidate identity changed: "
            f"{observed} observed, {configured} expected"
        )


def run(cfg: dict[str, Any], config_path: Path) -> dict[str, Path]:
    stage = cfg["nova_family_renomination"]
    output_dir = _configured_path(stage["output_dir"])
    candidate_path = _configured_path(stage["candidate_manifest"])
    pairs, pair_hashes = _load_pairs()
    transcript_ids = set(pairs["transcript_id_1"].astype(str)) | set(
        pairs["transcript_id_2"].astype(str)
    )
    models = _load_transcript_models(
        _configured_path(stage["gtf_cache"]), transcript_ids
    )
    scores, clusters = _score_transcripts(
        transcript_ids,
        models,
        _configured_path(stage["genome_fasta"]),
        stage,
    )
    scored_pairs = _pair_scores(pairs, scores)
    pair_calls, gene_calls = _tagged_pair_calls(scored_pairs, stage)
    regulons = _fit_regulons(gene_calls, stage)
    candidates = _freeze_candidates(
        pair_calls, regulons, stage, pair_hashes, config_path
    )
    observed = {
        "transcripts": len(scores),
        "eligible_pairs": int(pair_calls["eligible_pair"].sum()),
        "regulon_nominations": int(
            (
                regulons["model_status"].eq("fit")
                & regulons["qvalue"].le(float(stage["regulon_q_max"]))
                & regulons["coefficient"].gt(0)
            ).sum()
        ),
        "candidate_rows": len(candidates),
    }
    candidate_identity_sha256 = _candidate_identity_sha256(candidates)
    _guard_expected(stage, observed)
    _guard_candidate_identity(stage, candidate_identity_sha256)
    ensure_dir(output_dir)
    outputs = {
        "scores": output_dir / "nova_family_transcript_scores.parquet",
        "clusters": output_dir / "nova_family_clusters.parquet",
        "pairs": output_dir / "nova_family_pair_calls.parquet",
        "genes": output_dir / "nova_family_gene_calls.parquet",
        "regulons": output_dir / "nova_family_regulon.parquet",
        "candidates": candidate_path,
        "json": output_dir / "nova_family_renomination.json",
        "report": output_dir / "NOVA_FAMILY_RENOMINATION.md",
    }
    _atomic_parquet(scores, outputs["scores"])
    _atomic_parquet(clusters, outputs["clusters"])
    _atomic_parquet(pair_calls, outputs["pairs"])
    _atomic_parquet(gene_calls, outputs["genes"])
    _atomic_parquet(regulons, outputs["regulons"])
    _atomic_parquet(candidates, outputs["candidates"])
    summary = {
        "config": str(config_path),
        "config_sha256": sha256(config_path),
        "gtf_cache_sha256": sha256(_configured_path(stage["gtf_cache"])),
        "genome_fai_sha256": sha256(
            Path(f"{_configured_path(stage['genome_fasta'])}.fai")
        ),
        "pair_input_sha256": pair_hashes,
        "observed": observed,
        "candidate_identity_sha256": candidate_identity_sha256,
        "stage": stage,
        "output_sha256": {
            name: sha256(path)
            for name, path in outputs.items()
            if name not in {"json", "report"}
        },
    }
    _atomic_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", outputs["json"])
    fitted = regulons[regulons["model_status"].eq("fit")]
    nominations = fitted[
        fitted["qvalue"].le(float(stage["regulon_q_max"])) & fitted["coefficient"].gt(0)
    ].sort_values("qvalue")
    lines = [
        "# Opportunity-controlled NOVA-family re-nomination",
        "",
        "This analysis is independent of the invalidated transcript-wide NOVA1 presence nomination. Motif evidence is labeled NOVA-family.",
        "",
        f"- Switch transcripts scored: **{observed['transcripts']:,}**",
        f"- Opportunity-eligible region-specific pairs: **{observed['eligible_pairs']:,}**",
        f"- Adjusted module nominations at q <= {stage['regulon_q_max']}: **{observed['regulon_nominations']:,}**",
        f"- Frozen candidate rows: **{observed['candidate_rows']:,}**",
        "",
        "Module enrichment uses gene-level logistic regression adjusted for total intronic opportunity, eligible-pair count, and transcript count. A nomination additionally requires a positive module coefficient.",
        "",
        "## Nominations",
        "",
    ]
    if nominations.empty:
        lines.append(
            "No adjusted NOVA-family regulon passed the prespecified threshold."
        )
    else:
        display = nominations[
            [
                "region",
                "tree",
                "module_id",
                "go_invisible",
                "module_genes",
                "module_switched",
                "odds_ratio",
                "pvalue",
                "qvalue",
            ]
        ]
        headers = [str(column) for column in display.columns]
        lines.append("| " + " | ".join(headers) + " |")
        lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
        for values in display.itertuples(index=False, name=None):
            lines.append("| " + " | ".join(str(value) for value in values) + " |")
    lines.append("")
    _atomic_text("\n".join(lines), outputs["report"])
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = _configured_path(args.config)
    outputs = run(load_yaml(config_path), config_path)
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
