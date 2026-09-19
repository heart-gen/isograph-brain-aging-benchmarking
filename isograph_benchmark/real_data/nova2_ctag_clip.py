"""Test frozen NOVA-family candidates with public neuronal NOVA2 cTag-CLIP."""

from __future__ import annotations

import argparse
import bisect
import copy
import gzip
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

from isograph_benchmark.config import load_yaml
from isograph_benchmark.paths import ensure_dir
from isograph_benchmark.real_data.neuronal_clip_fetch import sha256
from isograph_benchmark.real_data.neuronal_clip_validation import (
    DEFAULT_CONFIG,
    _configured_path,
    prepare_windows,
)
from isograph_benchmark.real_data.nova_family_renomination import (
    _candidate_identity_sha256,
)


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(f".{path.name}.tmp")
    frame.to_parquet(temporary, index=False, compression="zstd")
    os.replace(temporary, path)


def _atomic_text(value: str, path: Path) -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(value)
    os.replace(temporary, path)


def _markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "No estimable rows."
    headers = [str(column) for column in frame.columns]
    rows = [headers, ["---"] * len(headers)]
    for values in frame.itertuples(index=False, name=None):
        rows.append([str(value) for value in values])
    return "\n".join("| " + " | ".join(row) + " |" for row in rows)


def _flip_strand(strand: str) -> str:
    return {"+": "-", "-": "+"}.get(strand, strand)


def _target_intervals(frame: pd.DataFrame) -> dict[str, tuple[list[int], list[Any]]]:
    grouped: dict[str, tuple[list[int], list[Any]]] = {}
    for chrom, chrom_frame in frame.groupby("chrom", sort=False):
        ordered = chrom_frame.sort_values(["start", "end", "interval_id"])
        rows = list(ordered.itertuples(index=False))
        grouped[str(chrom)] = ([int(row.start) for row in rows], rows)
    return grouped


def _map_chain(path: Path, intervals: pd.DataFrame) -> pd.DataFrame:
    required = {"interval_id", "chrom", "start", "end", "strand"}
    missing = sorted(required - set(intervals.columns))
    if missing:
        raise ValueError(f"Intervals missing columns: {missing}")
    if intervals["interval_id"].duplicated().any():
        raise ValueError("interval_id must be unique")
    targets = _target_intervals(intervals)
    mappings: dict[str, list[dict[str, Any]]] = {
        str(value): [] for value in intervals["interval_id"]
    }
    opener = gzip.open if path.suffix == ".gz" else open
    header: dict[str, Any] | None = None
    target_position = 0
    query_position = 0
    with opener(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                header = None
                continue
            fields = line.split()
            if fields[0] == "chain":
                if len(fields) != 13:
                    raise ValueError(f"Malformed chain header in {path}: {line}")
                header = {
                    "score": int(fields[1]),
                    "target_chrom": fields[2],
                    "query_chrom": fields[7],
                    "query_size": int(fields[8]),
                    "query_strand": fields[9],
                    "chain_id": fields[12],
                }
                target_position = int(fields[5])
                query_position = int(fields[10])
                continue
            if header is None:
                raise ValueError(f"Chain block without header in {path}")
            size = int(fields[0])
            target_end = target_position + size
            target = targets.get(str(header["target_chrom"]))
            if target is not None:
                starts, rows = target
                first = bisect.bisect_left(starts, target_position)
                final = bisect.bisect_left(starts, target_end)
                for row in rows[first:final]:
                    if int(row.end) > target_end:
                        continue
                    relative_start = int(row.start) - target_position
                    relative_end = int(row.end) - target_position
                    if header["query_strand"] == "+":
                        query_start = query_position + relative_start
                        query_end = query_position + relative_end
                        query_interval_strand = str(row.strand)
                    else:
                        query_start = int(header["query_size"]) - (
                            query_position + relative_end
                        )
                        query_end = int(header["query_size"]) - (
                            query_position + relative_start
                        )
                        query_interval_strand = _flip_strand(str(row.strand))
                    mappings[str(row.interval_id)].append(
                        {
                            "mapped_chrom": str(header["query_chrom"]),
                            "mapped_start": query_start,
                            "mapped_end": query_end,
                            "mapped_strand": query_interval_strand,
                            "chain_score": int(header["score"]),
                            "chain_id": str(header["chain_id"]),
                        }
                    )
            target_position = target_end
            query_position += size
            if len(fields) == 3:
                target_position += int(fields[1])
                query_position += int(fields[2])
            elif len(fields) != 1:
                raise ValueError(f"Malformed chain block in {path}: {line}")

    rows: list[dict[str, Any]] = []
    for interval in intervals.itertuples(index=False):
        options = mappings[str(interval.interval_id)]
        if not options:
            rows.append(
                {
                    "interval_id": str(interval.interval_id),
                    "mapping_status": "unmapped_or_crosses_alignment_gap",
                }
            )
            continue
        options.sort(
            key=lambda row: (
                -int(row["chain_score"]),
                str(row["mapped_chrom"]),
                int(row["mapped_start"]),
                int(row["mapped_end"]),
                str(row["chain_id"]),
            )
        )
        best_score = int(options[0]["chain_score"])
        best = [row for row in options if int(row["chain_score"]) == best_score]
        selected = best[0]
        rows.append(
            {
                "interval_id": str(interval.interval_id),
                **selected,
                "mapping_status": (
                    "mapped_unique_best" if len(best) == 1 else "ambiguous_best_chain"
                ),
                "mapping_count": len(options),
                "best_mapping_count": len(best),
            }
        )
    return pd.DataFrame(rows)


def _reciprocal_map(
    windows: pd.DataFrame,
    forward_chain: Path,
    reverse_chain: Path,
    minimum_overlap: float,
) -> pd.DataFrame:
    source = windows[["window_id", "chrom", "start", "end", "strand"]].rename(
        columns={"window_id": "interval_id"}
    )
    forward = (
        _map_chain(forward_chain, source)
        .add_prefix("forward_")
        .rename(columns={"forward_interval_id": "window_id"})
    )
    mapped = forward[forward["forward_mapping_status"].eq("mapped_unique_best")]
    reverse_input = mapped[
        [
            "window_id",
            "forward_mapped_chrom",
            "forward_mapped_start",
            "forward_mapped_end",
            "forward_mapped_strand",
        ]
    ].rename(
        columns={
            "window_id": "interval_id",
            "forward_mapped_chrom": "chrom",
            "forward_mapped_start": "start",
            "forward_mapped_end": "end",
            "forward_mapped_strand": "strand",
        }
    )
    reverse = (
        _map_chain(reverse_chain, reverse_input)
        .add_prefix("reverse_")
        .rename(columns={"reverse_interval_id": "window_id"})
    )
    result = windows.merge(forward, on="window_id", how="left", validate="one_to_one")
    result = result.merge(reverse, on="window_id", how="left", validate="one_to_one")
    result["reciprocal_overlap_fraction"] = 0.0
    reverse_mapped = result["reverse_mapping_status"].eq("mapped_unique_best")
    overlap = np.maximum(
        0,
        np.minimum(result["end"], result["reverse_mapped_end"])
        - np.maximum(result["start"], result["reverse_mapped_start"]),
    )
    result.loc[reverse_mapped, "reciprocal_overlap_fraction"] = overlap[
        reverse_mapped
    ] / (result.loc[reverse_mapped, "end"] - result.loc[reverse_mapped, "start"])
    result["reciprocal_mapped"] = (
        result["forward_mapping_status"].eq("mapped_unique_best")
        & reverse_mapped
        & result["reverse_mapped_chrom"].eq(result["chrom"])
        & result["reverse_mapped_strand"].eq(result["strand"])
        & result["reciprocal_overlap_fraction"].ge(minimum_overlap)
    )
    result["orthology_status"] = np.select(
        [
            ~result["forward_mapping_status"].eq("mapped_unique_best"),
            ~reverse_mapped,
            result["reverse_mapped_chrom"].ne(result["chrom"])
            | result["reverse_mapped_strand"].ne(result["strand"]),
            result["reciprocal_overlap_fraction"].lt(minimum_overlap),
        ],
        [
            "forward_not_uniquely_mapped",
            "reverse_not_uniquely_mapped",
            "reciprocal_chrom_or_strand_mismatch",
            "reciprocal_overlap_below_threshold",
        ],
        default="reciprocal_pass",
    )
    return result


def _bedgraph_calls(
    windows: pd.DataFrame, path: Path, context_id: str, minimum_coverage: float
) -> pd.DataFrame:
    mapped = windows[windows["reciprocal_mapped"]].copy()
    mapped["chrom"] = mapped["forward_mapped_chrom"]
    mapped["start"] = mapped["forward_mapped_start"].astype(int)
    mapped["end"] = mapped["forward_mapped_end"].astype(int)
    grouped: dict[str, tuple[list[int], list[Any], int]] = {}
    for chrom, frame in mapped.groupby("chrom", sort=False):
        ordered = frame.sort_values(["start", "end", "window_id"])
        rows = list(ordered.itertuples(index=False))
        grouped[str(chrom)] = (
            [int(row.start) for row in rows],
            rows,
            max(int(row.end) - int(row.start) for row in rows),
        )
    summaries = {
        str(row.window_id): {"covered_nt": 0, "coverage_area": 0.0, "max_coverage": 0.0}
        for row in mapped.itertuples(index=False)
    }
    with gzip.open(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.split("\t")
            if len(fields) != 4:
                raise ValueError(f"Malformed bedGraph row in {path}: {line}")
            chrom, start_value, end_value, coverage_value = fields
            target = grouped.get(chrom)
            if target is None:
                continue
            start = int(start_value)
            end = int(end_value)
            coverage = float(coverage_value)
            starts, rows, maximum_length = target
            first = bisect.bisect_left(starts, start - maximum_length + 1)
            final = bisect.bisect_left(starts, end)
            for row in rows[first:final]:
                overlap = max(0, min(int(row.end), end) - max(int(row.start), start))
                if not overlap:
                    continue
                summary = summaries[str(row.window_id)]
                summary["covered_nt"] += overlap
                summary["coverage_area"] += overlap * coverage
                summary["max_coverage"] = max(summary["max_coverage"], coverage)
    rows = []
    for window in windows.itertuples(index=False):
        summary = summaries.get(
            str(window.window_id),
            {"covered_nt": 0, "coverage_area": 0.0, "max_coverage": 0.0},
        )
        length = int(window.end) - int(window.start)
        covered_nt = min(length, int(summary["covered_nt"]))
        rows.append(
            {
                "context_id": context_id,
                "window_id": str(window.window_id),
                "callable": bool(window.reciprocal_mapped),
                "bound": bool(
                    window.reciprocal_mapped
                    and float(summary["max_coverage"]) >= minimum_coverage
                ),
                "covered_nt": covered_nt,
                "covered_fraction": covered_nt / length if length else math.nan,
                "mean_coverage": (
                    float(summary["coverage_area"]) / length if length else math.nan
                ),
                "max_coverage": float(summary["max_coverage"]),
                "callability_basis": "reciprocal_liftover_only",
                "signal_resolution": "pooled_only",
                "signal_type": "study_processed_unique_tag_coverage",
                "source_path": str(path),
                "source_sha256": sha256(path),
            }
        )
    return pd.DataFrame(rows)


def _effect_row(frame: pd.DataFrame, unit: str) -> dict[str, Any]:
    callable_frame = frame[frame["case_callable"] & frame["control_callable"]]
    case = callable_frame["case_bound"].astype(bool)
    control = callable_frame["control_bound"].astype(bool)
    case_only = int((case & ~control).sum())
    control_only = int((~case & control).sum())
    discordant = case_only + control_only
    odds_ratio = (case_only + 0.5) / (control_only + 0.5)
    standard_error = math.sqrt(1 / (case_only + 0.5) + 1 / (control_only + 0.5))
    return {
        "statistical_unit": unit,
        "callable_units": int(len(callable_frame)),
        "case_bound": int(case.sum()),
        "control_bound": int(control.sum()),
        "case_bound_fraction": float(case.mean()) if len(case) else math.nan,
        "control_bound_fraction": float(control.mean()) if len(control) else math.nan,
        "absolute_bound_fraction_difference": (
            float(case.mean() - control.mean()) if len(case) else math.nan
        ),
        "case_only_discordant": case_only,
        "control_only_discordant": control_only,
        "informative_discordant": discordant,
        "matched_odds_ratio_haldane": odds_ratio,
        "ci_low_haldane": math.exp(math.log(odds_ratio) - 1.96 * standard_error),
        "ci_high_haldane": math.exp(math.log(odds_ratio) + 1.96 * standard_error),
        "exact_pvalue": (
            float(binomtest(case_only, discordant, 0.5).pvalue)
            if discordant
            else math.nan
        ),
    }


def _paired_calls(
    calls: pd.DataFrame, windows: pd.DataFrame, membership: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    joined = calls.merge(
        windows[
            [
                "window_id",
                "candidate_id",
                "gene_id",
                "match_id",
                "window_role",
                "window_width",
                "match_status",
            ]
        ],
        on="window_id",
        how="inner",
        validate="many_to_one",
    )
    joined = joined[joined["match_status"].eq("matched")].copy()
    pairs = joined.pivot(
        index=["context_id", "candidate_id", "gene_id", "match_id", "window_width"],
        columns="window_role",
        values=["callable", "bound"],
    )
    required = [
        (field, role) for field in ["callable", "bound"] for role in ["case", "control"]
    ]
    if any(column not in pairs.columns for column in required):
        raise ValueError("No complete case/control cTag window pairs")
    pairs = pairs.reset_index()
    pairs.columns = [
        (
            "_".join(str(value) for value in column if value).rstrip("_")
            if isinstance(column, tuple)
            else str(column)
        )
        for column in pairs.columns
    ]
    group_columns = ["context_id", "candidate_id", "gene_id", "window_width"]
    candidates = pairs.groupby(group_columns, as_index=False).agg(
        case_callable=("callable_case", "any"),
        control_callable=("callable_control", "any"),
        case_bound=("bound_case", "any"),
        control_bound=("bound_control", "any"),
        callable_matched_sets=("match_id", "nunique"),
    )
    candidate_metadata = membership[
        ["candidate_id", "tree", "region", "module_id", "go_invisible", "gene"]
    ].drop_duplicates("candidate_id")
    candidates = candidates.merge(
        candidate_metadata, on="candidate_id", how="left", validate="many_to_one"
    )
    genes = candidates.groupby(
        [
            "context_id",
            "tree",
            "region",
            "module_id",
            "go_invisible",
            "gene",
            "window_width",
        ],
        as_index=False,
    ).agg(
        case_callable=("case_callable", "any"),
        control_callable=("control_callable", "any"),
        case_bound=("case_bound", "any"),
        control_bound=("control_bound", "any"),
        candidate_pairs=("candidate_id", "nunique"),
    )
    return pairs, candidates, genes


def _effect_tables(
    pairs: pd.DataFrame,
    genes: pd.DataFrame,
    primary_width: int,
    min_discordant: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    overall_rows = []
    for (context_id, width), frame in genes.groupby(["context_id", "window_width"]):
        row = _effect_row(frame, "gene")
        overall_rows.append(
            {"context_id": context_id, "window_width": int(width), **row}
        )
    overall = pd.DataFrame(overall_rows)
    overall["inference_status"] = np.where(
        overall["informative_discordant"].ge(min_discordant),
        "descriptive_pooled_signal_estimable",
        "descriptive_pooled_signal_sparse_discordance",
    )

    module_rows = []
    group = [
        "context_id",
        "tree",
        "region",
        "module_id",
        "go_invisible",
        "window_width",
    ]
    for keys, frame in genes.groupby(group):
        row = _effect_row(frame, "gene")
        module_rows.append({**dict(zip(group, keys)), **row})
    modules = pd.DataFrame(module_rows)
    modules["qvalue"] = np.nan
    primary = (
        modules["window_width"].eq(primary_width) & modules["exact_pvalue"].notna()
    )
    if primary.any():
        modules.loc[primary, "qvalue"] = multipletests(
            modules.loc[primary, "exact_pvalue"], method="fdr_bh"
        )[1]
    modules["inference_status"] = np.where(
        modules["informative_discordant"].ge(min_discordant),
        "descriptive_pooled_signal_estimable",
        "descriptive_pooled_signal_sparse_discordance",
    )

    window_rows = []
    for (context_id, width), frame in pairs.groupby(["context_id", "window_width"]):
        renamed = frame.rename(
            columns={
                "callable_case": "case_callable",
                "callable_control": "control_callable",
                "bound_case": "case_bound",
                "bound_control": "control_bound",
            }
        )
        row = _effect_row(renamed, "matched_window")
        window_rows.append(
            {"context_id": context_id, "window_width": int(width), **row}
        )
    return overall, modules, pd.DataFrame(window_rows)


def _cell_type_interaction(genes: pd.DataFrame, stage: dict[str, Any]) -> pd.DataFrame:
    primary = genes[
        genes["window_width"].eq(int(stage["primary_width"]))
        & genes["context_id"].isin(stage["primary_contexts"])
        & genes["case_callable"]
        & genes["control_callable"]
    ].copy()
    primary["delta"] = primary["case_bound"].astype(int) - primary[
        "control_bound"
    ].astype(int)
    index = ["tree", "region", "module_id", "gene"]
    pivot = primary.pivot(index=index, columns="context_id", values="delta")
    first, second = [str(value) for value in stage["primary_contexts"]]
    if first not in pivot or second not in pivot:
        return pd.DataFrame()
    paired = pivot[[first, second]].dropna()
    difference = paired[first] - paired[second]
    positive = int(difference.gt(0).sum())
    negative = int(difference.lt(0).sum())
    informative = positive + negative
    return pd.DataFrame(
        [
            {
                "context_1": first,
                "context_2": second,
                "jointly_callable_genes": int(len(paired)),
                "mean_case_control_delta_context_1": float(paired[first].mean()),
                "mean_case_control_delta_context_2": float(paired[second].mean()),
                "positive_difference": positive,
                "negative_difference": negative,
                "informative_differences": informative,
                "exact_sign_pvalue": (
                    float(binomtest(positive, informative, 0.5).pvalue)
                    if informative
                    else math.nan
                ),
                "interpretation_scope": "pooled_signal_cell_type_interaction",
            }
        ]
    )


def _configured_bedgraphs(cfg: dict[str, Any], dataset_id: str) -> dict[str, Path]:
    root = _configured_path(cfg["output_root"])
    dataset = next(
        (row for row in cfg["datasets"] if str(row["dataset_id"]) == dataset_id),
        None,
    )
    if dataset is None:
        raise ValueError(f"Dataset not configured: {dataset_id}")
    paths: dict[str, Path] = {}
    for sample in dataset.get("samples", []):
        files = [
            row for row in sample.get("files", []) if row.get("format") == "bedgraph_gz"
        ]
        if len(files) != 1:
            raise ValueError(f"{sample['sample_id']} must have one pooled bedGraph")
        paths[str(sample["context_id"])] = root / str(files[0]["path"])
    return paths


def _prepare_human_windows(
    cfg: dict[str, Any], config_path: Path, stage: dict[str, Any], output_dir: Path
) -> dict[str, Path]:
    derived = copy.deepcopy(cfg)
    derived["candidate_manifest"] = str(stage["candidate_manifest"])
    derived["window_stage"]["expected_candidate_rows"] = int(
        stage["expected_candidate_rows"]
    )
    # Regions can nominate the same RBP/gene/transcript pair, so unique canonical
    # candidates may be fewer than manifest rows; legacy configs pinned only rows.
    derived["window_stage"]["expected_unique_candidates"] = int(
        stage.get("expected_unique_candidates", stage["expected_candidate_rows"])
    )
    derived["window_stage"]["output_dir"] = str(output_dir)
    return prepare_windows(derived, config_path, output_dir, set(), None)


def run(cfg: dict[str, Any], config_path: Path) -> dict[str, Path]:
    stage = cfg["nova2_ctag_stage"]
    stage["primary_width"] = int(cfg["window_stage"]["primary_width"])
    stage["min_informative_discordant_sets"] = int(
        cfg["overlap_stage"]["min_informative_discordant_sets"]
    )
    output_dir = _configured_path(stage["output_dir"])
    ensure_dir(output_dir)
    candidate_path = _configured_path(stage["candidate_manifest"])
    candidates = pd.read_parquet(candidate_path)
    if len(candidates) != int(stage["expected_candidate_rows"]):
        raise RuntimeError("NOVA-family candidate row count changed")
    identity = _candidate_identity_sha256(candidates)
    if identity != str(stage["expected_candidate_identity_sha256"]):
        raise RuntimeError(f"NOVA-family candidate identity changed: {identity}")

    human_dir = output_dir / "human_windows"
    human_outputs = _prepare_human_windows(cfg, config_path, stage, human_dir)
    human_windows = pd.read_parquet(human_outputs["windows"])
    matched = human_windows[human_windows["match_status"].eq("matched")].copy()
    reciprocal = _reciprocal_map(
        matched,
        _configured_path(stage["forward_chain"]),
        _configured_path(stage["reverse_chain"]),
        float(stage["reciprocal_overlap_min"]),
    )
    membership = pd.read_parquet(human_outputs["membership"])

    bedgraphs = _configured_bedgraphs(cfg, str(stage["dataset_id"]))
    expected_contexts = set(stage["primary_contexts"]) | set(
        stage["secondary_contexts"]
    )
    if set(bedgraphs) != expected_contexts:
        raise RuntimeError(f"Configured cTag contexts changed: {sorted(bedgraphs)}")
    call_frames = []
    context_call_outputs: dict[str, Path] = {}
    context_call_dir = output_dir / "context_calls"
    ensure_dir(context_call_dir)
    for context_id in sorted(bedgraphs):
        path = bedgraphs[context_id]
        if not path.is_file():
            raise FileNotFoundError(path)
        print(f"scoring {context_id}: {path}", flush=True)
        context_calls = _bedgraph_calls(
            reciprocal,
            path,
            context_id,
            float(stage["minimum_coverage"]),
        )
        context_path = context_call_dir / f"{context_id}.parquet"
        _atomic_parquet(context_calls, context_path)
        context_call_outputs[context_id] = context_path
        call_frames.append(context_calls)
    calls = pd.concat(call_frames, ignore_index=True)
    pairs, candidate_calls, gene_calls = _paired_calls(calls, reciprocal, membership)
    overall, modules, window_effects = _effect_tables(
        pairs,
        gene_calls,
        int(stage["primary_width"]),
        int(stage["min_informative_discordant_sets"]),
    )
    interaction = _cell_type_interaction(gene_calls, stage)

    outputs = {
        "reciprocal_windows": output_dir / "reciprocal_windows.parquet",
        "overlap_calls": output_dir / "overlap_calls.parquet",
        "pair_calls": output_dir / "matched_pair_calls.parquet",
        "candidate_support": output_dir / "candidate_support.parquet",
        "gene_support": output_dir / "gene_support.parquet",
        "overall_effects": output_dir / "overall_effects.parquet",
        "module_effects": output_dir / "module_effects.parquet",
        "window_effects": output_dir / "window_effects.parquet",
        "interaction": output_dir / "cell_type_interaction.parquet",
        "json": output_dir / "nova2_ctag_validation.json",
        "report": output_dir / "NOVA2_CTAG_VALIDATION.md",
    }
    for name, frame in [
        ("reciprocal_windows", reciprocal),
        ("overlap_calls", calls),
        ("pair_calls", pairs),
        ("candidate_support", candidate_calls),
        ("gene_support", gene_calls),
        ("overall_effects", overall),
        ("module_effects", modules),
        ("window_effects", window_effects),
        ("interaction", interaction),
    ]:
        _atomic_parquet(frame, outputs[name])

    provenance = {
        "config": str(config_path),
        "config_sha256": sha256(config_path),
        "candidate_manifest": str(candidate_path),
        "candidate_manifest_sha256": sha256(candidate_path),
        "candidate_identity_sha256": identity,
        "source_sha256": sha256(Path(__file__)),
        "forward_chain_sha256": sha256(_configured_path(stage["forward_chain"])),
        "reverse_chain_sha256": sha256(_configured_path(stage["reverse_chain"])),
        "bedgraph_sha256": {key: sha256(value) for key, value in bedgraphs.items()},
        "context_call_sha256": {
            key: sha256(value) for key, value in context_call_outputs.items()
        },
        "stage": stage,
        "observed": {
            "human_matched_window_rows": int(len(matched)),
            "reciprocal_window_rows": int(reciprocal["reciprocal_mapped"].sum()),
            "reciprocal_window_fraction": float(reciprocal["reciprocal_mapped"].mean()),
            "candidate_support_rows": int(len(candidate_calls)),
            "gene_support_rows": int(len(gene_calls)),
        },
        "output_sha256": {
            name: sha256(path)
            for name, path in outputs.items()
            if name not in {"json", "report"}
        },
    }
    _atomic_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", outputs["json"]
    )

    primary = overall[overall["window_width"].eq(int(stage["primary_width"]))].copy()
    primary_display = primary[
        [
            "context_id",
            "callable_units",
            "case_bound_fraction",
            "control_bound_fraction",
            "absolute_bound_fraction_difference",
            "informative_discordant",
            "matched_odds_ratio_haldane",
            "ci_low_haldane",
            "ci_high_haldane",
            "exact_pvalue",
            "inference_status",
        ]
    ]
    module_display = modules[modules["window_width"].eq(int(stage["primary_width"]))][
        [
            "context_id",
            "module_id",
            "go_invisible",
            "callable_units",
            "informative_discordant",
            "matched_odds_ratio_haldane",
            "exact_pvalue",
            "qvalue",
        ]
    ]
    lines = [
        "# NOVA2 neuronal cTag-CLIP validation",
        "",
        "This orthogonal Tier-2 analysis tests the frozen NOVA-family candidates against exact-RBP mouse NOVA2 cTag-CLIP. It does not convert mouse occupancy into evidence of adult human cortical occupancy.",
        "",
        f"- Frozen candidate rows: **{len(candidates):,}**",
        f"- Matched human window rows: **{len(matched):,}**",
        f"- Reciprocally mapped window rows: **{int(reciprocal['reciprocal_mapped'].sum()):,}** ({reciprocal['reciprocal_mapped'].mean():.1%})",
        "- Signal source: study-processed pooled unique-tag coverage from three biological replicates per context.",
        "- Replicate limitation: processed replicate-level peaks are unavailable; all p-values are labeled descriptive pooled-signal tests.",
        "- Callability limitation: no matched input is provided; absence of coverage is not interpreted as proof of no binding.",
        "",
        "## Gene-level primary-width effects",
        "",
        _markdown_table(primary_display),
        "",
        "## Module-level primary-width effects",
        "",
        _markdown_table(module_display),
        "",
        "## Cortical cell-type interaction",
        "",
        _markdown_table(interaction),
        "",
        "Candidates remain labeled NOVA_FAMILY in discovery. Exact NOVA2 is used only for this orthogonal occupancy evidence tier.",
        "",
    ]
    _atomic_text("\n".join(lines), outputs["report"])
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = Path(args.config)
    outputs = run(load_yaml(config_path), config_path)
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
