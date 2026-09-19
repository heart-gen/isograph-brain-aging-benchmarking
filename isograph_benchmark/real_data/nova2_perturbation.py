"""Localize Nova2-cKO splicing responses to frozen NOVA-family switch windows."""

from __future__ import annotations

import argparse
import bisect
import json
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
    _stable_id,
)
from isograph_benchmark.real_data.nova2_ctag_clip import (
    _atomic_parquet,
    _atomic_text,
    _effect_row,
    _markdown_table,
    _reciprocal_map,
)
from isograph_benchmark.real_data.nova_family_renomination import (
    _candidate_identity_sha256,
)


BED12_COLUMNS = [
    "chrom",
    "chrom_start",
    "chrom_end",
    "event_id",
    "score",
    "strand",
    "thick_start",
    "thick_end",
    "item_rgb",
    "block_count",
    "block_sizes",
    "block_starts",
]


def _read_event_statistics(
    path: Path, context_id: str, fdr_max: float, absolute_delta_min: float
) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="gzip")
    frame.columns = [str(column).strip() for column in frame.columns]
    required = {
        "Ensemble_gene//Gene_symbol",
        "chr",
        "chr.Start",
        "chr.End",
        "AS.name",
        "strand",
        "AS.type",
        "delta.exon.inclusion.rate(dI)",
        "pvalue",
        "FDR",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Quantas table {path} is missing columns: {missing}")
    frame = frame.reset_index(names="source_stat_row")
    frame = frame.rename(
        columns={
            "Ensemble_gene//Gene_symbol": "source_gene_annotation",
            "chr": "source_chrom",
            "chr.Start": "source_start",
            "chr.End": "source_end",
            "AS.name": "event_id",
            "strand": "source_strand",
            "AS.type": "event_type",
            "delta.exon.inclusion.rate(dI)": "delta_inclusion",
            "FDR": "fdr",
        }
    )
    for column in ["source_start", "source_end", "delta_inclusion", "pvalue", "fdr"]:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    if frame["event_id"].isna().any() or frame["event_id"].astype(str).eq("").any():
        raise ValueError(f"Quantas table {path} contains empty event IDs")
    if ~frame["source_strand"].isin(["+", "-"]).all():
        raise ValueError(f"Quantas table {path} contains invalid strands")
    frame["event_id"] = frame["event_id"].astype(str)
    frame["absolute_delta_inclusion"] = frame["delta_inclusion"].abs()
    frame["source_stat_rows"] = frame.groupby("event_id")["event_id"].transform("size")
    frame["passes_event_thresholds"] = frame["fdr"].le(fdr_max) & frame[
        "absolute_delta_inclusion"
    ].ge(absolute_delta_min)
    frame = frame.sort_values(
        ["event_id", "fdr", "absolute_delta_inclusion", "source_stat_row"],
        ascending=[True, True, False, True],
        kind="mergesort",
    ).drop_duplicates("event_id", keep="first")
    frame.insert(0, "context_id", context_id)
    frame["statistics_selection_rule"] = (
        "lowest_fdr_then_largest_absolute_delta_then_source_row"
    )
    columns = [
        "context_id",
        "event_id",
        "source_gene_annotation",
        "source_chrom",
        "source_start",
        "source_end",
        "source_strand",
        "event_type",
        "delta_inclusion",
        "absolute_delta_inclusion",
        "pvalue",
        "fdr",
        "passes_event_thresholds",
        "source_stat_rows",
        "source_stat_row",
        "statistics_selection_rule",
    ]
    return (
        frame[columns]
        .sort_values("event_id", kind="mergesort")
        .reset_index(drop=True)
    )


def _parse_integer_list(value: Any, expected: int, label: str) -> list[int]:
    values = [item for item in str(value).rstrip(",").split(",") if item]
    try:
        parsed = [int(item) for item in values]
    except ValueError as error:
        raise ValueError(f"Malformed BED12 {label}: {value!r}") from error
    if len(parsed) != expected:
        raise ValueError(
            f"BED12 {label} has {len(parsed)} values but block_count is {expected}"
        )
    return parsed


def _read_bed12(path: Path, context_id: str) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        names=BED12_COLUMNS,
        header=None,
        dtype={"chrom": str, "event_id": str, "strand": str},
    ).reset_index(names="source_bed_row")
    for column in ["chrom_start", "chrom_end", "block_count"]:
        frame[column] = pd.to_numeric(frame[column], errors="raise").astype(int)
    if ~frame["strand"].isin(["+", "-"]).all():
        raise ValueError(f"BED12 table {path} contains invalid strands")
    if (frame["chrom_end"] <= frame["chrom_start"]).any():
        raise ValueError(f"BED12 table {path} contains non-positive intervals")
    if (frame["block_count"] < 2).any():
        raise ValueError(f"BED12 table {path} contains fewer than two blocks")
    frame["event_coordinate_rows"] = frame.groupby("event_id")["event_id"].transform(
        "size"
    )
    frame.insert(0, "context_id", context_id)
    return frame


def _event_flanks(
    statistics: pd.DataFrame, bed12: pd.DataFrame, widths: list[int]
) -> pd.DataFrame:
    statistics_ids = set(statistics["event_id"])
    coordinate_ids = set(bed12["event_id"])
    if statistics_ids != coordinate_ids:
        raise ValueError(
            "Quantas/BED12 unique event IDs do not reconcile: "
            f"statistics_only={len(statistics_ids - coordinate_ids)}, "
            f"bed12_only={len(coordinate_ids - statistics_ids)}"
        )
    eligible = statistics[statistics["passes_event_thresholds"]].copy()
    bed = bed12[bed12["event_id"].isin(eligible["event_id"])].copy()
    missing_coordinates = sorted(set(eligible["event_id"]) - set(bed["event_id"]))
    if missing_coordinates:
        raise ValueError(
            f"{len(missing_coordinates)} eligible Quantas events lack BED12 coordinates"
        )
    metadata = eligible.set_index("event_id")
    rows: list[dict[str, Any]] = []
    for bed_row in bed.itertuples(index=False):
        sizes = _parse_integer_list(
            bed_row.block_sizes, int(bed_row.block_count), "block_sizes"
        )
        starts = _parse_integer_list(
            bed_row.block_starts, int(bed_row.block_count), "block_starts"
        )
        if any(size <= 0 for size in sizes) or any(start < 0 for start in starts):
            raise ValueError(f"BED12 has invalid blocks for {bed_row.event_id}")
        blocks = sorted(
            (
                int(bed_row.chrom_start) + start,
                int(bed_row.chrom_start) + start + size,
            )
            for start, size in zip(starts, sizes)
        )
        if any(
            start < int(bed_row.chrom_start) or end > int(bed_row.chrom_end)
            for start, end in blocks
        ):
            raise ValueError(
                f"BED12 blocks exceed parent interval for {bed_row.event_id}"
            )
        stat = metadata.loc[str(bed_row.event_id)]
        coordinate_id = _stable_id(
            "E",
            str(bed_row.context_id),
            str(bed_row.event_id),
            str(bed_row.chrom),
            int(bed_row.chrom_start),
            int(bed_row.chrom_end),
            str(bed_row.strand),
            str(bed_row.block_sizes),
            str(bed_row.block_starts),
            int(bed_row.source_bed_row),
        )
        introns = [
            (left[1], right[0])
            for left, right in zip(blocks[:-1], blocks[1:])
            if left[1] < right[0]
        ]
        for genomic_index, (intron_start, intron_end) in enumerate(introns, start=1):
            transcript_index = (
                genomic_index
                if str(bed_row.strand) == "+"
                else len(introns) - genomic_index + 1
            )
            for width in widths:
                length = intron_end - intron_start
                if length <= 2 * width:
                    intervals = [(intron_start, intron_end, "short_intron")]
                else:
                    low_class = (
                        "donor_flank"
                        if str(bed_row.strand) == "+"
                        else "acceptor_flank"
                    )
                    high_class = (
                        "acceptor_flank"
                        if str(bed_row.strand) == "+"
                        else "donor_flank"
                    )
                    intervals = [
                        (intron_start, intron_start + width, low_class),
                        (intron_end - width, intron_end, high_class),
                    ]
                for start, end, structural_class in intervals:
                    window_id = _stable_id(
                        "F",
                        coordinate_id,
                        transcript_index,
                        width,
                        structural_class,
                        start,
                        end,
                    )
                    rows.append(
                        {
                            "context_id": str(bed_row.context_id),
                            "event_id": str(bed_row.event_id),
                            "event_coordinate_id": coordinate_id,
                            "window_id": window_id,
                            "source_bed_row": int(bed_row.source_bed_row),
                            "event_coordinate_rows": int(bed_row.event_coordinate_rows),
                            "source_gene_annotation": str(stat.source_gene_annotation),
                            "event_type": str(stat.event_type),
                            "delta_inclusion": float(stat.delta_inclusion),
                            "absolute_delta_inclusion": float(
                                stat.absolute_delta_inclusion
                            ),
                            "fdr": float(stat.fdr),
                            "chrom": str(bed_row.chrom),
                            "start": int(start),
                            "end": int(end),
                            "strand": str(bed_row.strand),
                            "window_width": int(width),
                            "structural_class": structural_class,
                            "intron_index": int(transcript_index),
                            "coordinate_system": "0-based_half-open",
                        }
                    )
    columns = [
        "context_id",
        "event_id",
        "event_coordinate_id",
        "window_id",
        "source_bed_row",
        "event_coordinate_rows",
        "source_gene_annotation",
        "event_type",
        "delta_inclusion",
        "absolute_delta_inclusion",
        "fdr",
        "chrom",
        "start",
        "end",
        "strand",
        "window_width",
        "structural_class",
        "intron_index",
        "coordinate_system",
    ]
    frame = pd.DataFrame(rows, columns=columns)
    if frame.empty:
        raise ValueError("No eligible intronic event flanks were derived")
    if frame["window_id"].duplicated().any():
        raise ValueError("Event flank window IDs are not unique")
    return frame.sort_values(
        ["context_id", "window_width", "chrom", "start", "end", "window_id"],
        kind="mergesort",
    ).reset_index(drop=True)


def _overlap_event_flanks(
    human_windows: pd.DataFrame, reciprocal_flanks: pd.DataFrame, context_id: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    windows = human_windows[
        human_windows["match_status"].eq("matched")
    ].copy()
    required = [
        "window_id",
        "chrom",
        "start",
        "end",
        "strand",
        "window_width",
        "structural_class",
    ]
    if windows["window_id"].duplicated().any():
        raise ValueError("Frozen human window IDs are not unique")
    groups: dict[tuple[Any, ...], tuple[list[int], list[Any], int]] = {}
    for key, frame in windows.groupby(
        ["chrom", "strand", "window_width", "structural_class"], sort=False
    ):
        ordered = frame.sort_values(["start", "end", "window_id"], kind="mergesort")
        records = list(ordered.itertuples(index=False))
        groups[key] = (
            [int(row.start) for row in records],
            records,
            max(int(row.end) - int(row.start) for row in records),
        )
    overlap_rows: list[dict[str, Any]] = []
    mapped = reciprocal_flanks[reciprocal_flanks["reciprocal_mapped"]]
    for event in mapped.itertuples(index=False):
        key = (
            str(event.forward_mapped_chrom),
            str(event.forward_mapped_strand),
            int(event.window_width),
            str(event.structural_class),
        )
        target = groups.get(key)
        if target is None:
            continue
        starts, records, maximum_length = target
        event_start = int(event.forward_mapped_start)
        event_end = int(event.forward_mapped_end)
        first = bisect.bisect_left(starts, event_start - maximum_length + 1)
        final = bisect.bisect_left(starts, event_end)
        for window in records[first:final]:
            overlap = max(
                0,
                min(int(window.end), event_end) - max(int(window.start), event_start),
            )
            if not overlap:
                continue
            overlap_rows.append(
                {
                    "context_id": context_id,
                    "window_id": str(window.window_id),
                    "event_flank_id": str(event.window_id),
                    "event_id": str(event.event_id),
                    "event_coordinate_id": str(event.event_coordinate_id),
                    "overlap_nt": int(overlap),
                    "event_fdr": float(event.fdr),
                    "event_delta_inclusion": float(event.delta_inclusion),
                    "event_absolute_delta_inclusion": float(
                        event.absolute_delta_inclusion
                    ),
                }
            )
    overlaps = pd.DataFrame(
        overlap_rows,
        columns=[
            "context_id",
            "window_id",
            "event_flank_id",
            "event_id",
            "event_coordinate_id",
            "overlap_nt",
            "event_fdr",
            "event_delta_inclusion",
            "event_absolute_delta_inclusion",
        ],
    )
    if overlaps.empty:
        summary = pd.DataFrame(
            {
                "window_id": pd.Series(dtype="string"),
                "localized_event_flanks": pd.Series(dtype="int64"),
                "localized_events": pd.Series(dtype="int64"),
                "maximum_overlap_nt": pd.Series(dtype="int64"),
                "minimum_event_fdr": pd.Series(dtype="float64"),
                "maximum_absolute_delta_inclusion": pd.Series(dtype="float64"),
            }
        )
    else:
        summary = overlaps.groupby("window_id", as_index=False).agg(
            localized_event_flanks=("event_flank_id", "nunique"),
            localized_events=("event_id", "nunique"),
            maximum_overlap_nt=("overlap_nt", "max"),
            minimum_event_fdr=("event_fdr", "min"),
            maximum_absolute_delta_inclusion=(
                "event_absolute_delta_inclusion",
                "max",
            ),
        )
    calls = windows[
        [
            "window_id",
            "candidate_id",
            "gene_id",
            "match_id",
            "window_role",
            "window_width",
            "structural_class",
        ]
    ].merge(summary, on="window_id", how="left", validate="one_to_one")
    for column in ["localized_event_flanks", "localized_events", "maximum_overlap_nt"]:
        calls[column] = (
            pd.to_numeric(calls[column], errors="coerce").fillna(0).astype(int)
        )
    for column in ["minimum_event_fdr", "maximum_absolute_delta_inclusion"]:
        calls[column] = pd.to_numeric(calls[column], errors="coerce")
    calls.insert(0, "context_id", context_id)
    calls["callable"] = True
    calls["localized"] = calls["localized_events"].gt(0)
    calls["localization_basis"] = (
        "significant_event_flank_reciprocal_liftover_same_strand_and_class"
    )
    return calls, overlaps


def _collapse_support(
    calls: pd.DataFrame, membership: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pairs = calls.pivot(
        index=["context_id", "candidate_id", "gene_id", "match_id", "window_width"],
        columns="window_role",
        values=["callable", "localized"],
    )
    required = [
        (field, role)
        for field in ["callable", "localized"]
        for role in ["case", "control"]
    ]
    if any(column not in pairs.columns for column in required):
        raise ValueError("No complete case/control perturbation window pairs")
    pairs = pairs.reset_index()
    pairs.columns = [
        "_".join(str(value) for value in column if value).rstrip("_")
        if isinstance(column, tuple)
        else str(column)
        for column in pairs.columns
    ]
    group_columns = ["context_id", "candidate_id", "gene_id", "window_width"]
    candidates = pairs.groupby(group_columns, as_index=False).agg(
        case_callable=("callable_case", "any"),
        control_callable=("callable_control", "any"),
        case_localized=("localized_case", "any"),
        control_localized=("localized_control", "any"),
        matched_sets=("match_id", "nunique"),
    )
    metadata = membership[
        ["candidate_id", "tree", "region", "module_id", "go_invisible", "gene"]
    ].drop_duplicates("candidate_id")
    candidates = candidates.merge(
        metadata, on="candidate_id", how="left", validate="many_to_one"
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
        dropna=False,
    ).agg(
        case_callable=("case_callable", "any"),
        control_callable=("control_callable", "any"),
        case_localized=("case_localized", "any"),
        control_localized=("control_localized", "any"),
        candidate_pairs=("candidate_id", "nunique"),
    )
    return pairs, candidates, genes


def _localization_effect(frame: pd.DataFrame, unit: str) -> dict[str, Any]:
    renamed = frame.rename(
        columns={
            "case_localized": "case_bound",
            "control_localized": "control_bound",
        }
    )
    row = _effect_row(renamed, unit)
    return {
        key.replace("bound", "localized"): value for key, value in row.items()
    }


def _effect_tables(
    pairs: pd.DataFrame,
    genes: pd.DataFrame,
    primary_width: int,
    minimum_discordant: int,
    primary_contexts: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    primary_context_set = set(primary_contexts or [])
    overall_rows = []
    for (context_id, width), frame in genes.groupby(["context_id", "window_width"]):
        overall_rows.append(
            {
                "context_id": context_id,
                "window_width": int(width),
                **_localization_effect(frame, "gene"),
            }
        )
    overall = pd.DataFrame(overall_rows)
    overall["inference_status"] = np.where(
        overall["informative_discordant"].ge(minimum_discordant),
        "formally_estimable",
        "not_estimable_sparse_discordance",
    )
    overall["analysis_scope"] = np.where(
        overall["context_id"].isin(primary_context_set), "primary", "secondary"
    )
    overall["positive_result"] = (
        overall["window_width"].eq(primary_width)
        & overall["analysis_scope"].eq("primary")
        & overall["inference_status"].eq("formally_estimable")
        & overall["matched_odds_ratio_haldane"].gt(1)
        & overall["exact_pvalue"].le(0.05)
    )

    group = [
        "context_id",
        "tree",
        "region",
        "module_id",
        "go_invisible",
        "window_width",
    ]
    module_rows = []
    for keys, frame in genes.groupby(group, dropna=False):
        module_rows.append(
            {**dict(zip(group, keys)), **_localization_effect(frame, "gene")}
        )
    modules = pd.DataFrame(module_rows)
    modules["qvalue"] = np.nan
    for _, index in modules[
        modules["window_width"].eq(primary_width) & modules["exact_pvalue"].notna()
    ].groupby("context_id").groups.items():
        modules.loc[index, "qvalue"] = multipletests(
            modules.loc[index, "exact_pvalue"], method="fdr_bh"
        )[1]
    modules["inference_status"] = np.where(
        modules["informative_discordant"].ge(minimum_discordant),
        "formally_estimable",
        "not_estimable_sparse_discordance",
    )
    modules["analysis_scope"] = np.where(
        modules["context_id"].isin(primary_context_set), "primary", "secondary"
    )
    modules["positive_result"] = (
        modules["window_width"].eq(primary_width)
        & modules["analysis_scope"].eq("primary")
        & modules["inference_status"].eq("formally_estimable")
        & modules["matched_odds_ratio_haldane"].gt(1)
        & modules["qvalue"].le(0.05)
    )

    window_rows = []
    for (context_id, width), frame in pairs.groupby(["context_id", "window_width"]):
        renamed = frame.rename(
            columns={
                "callable_case": "case_callable",
                "callable_control": "control_callable",
                "localized_case": "case_localized",
                "localized_control": "control_localized",
            }
        )
        window_rows.append(
            {
                "context_id": context_id,
                "window_width": int(width),
                **_localization_effect(renamed, "matched_window"),
            }
        )
    return overall, modules, pd.DataFrame(window_rows)


def _cross_context_replication(
    genes: pd.DataFrame, stage: dict[str, Any]
) -> pd.DataFrame:
    primary_contexts = [str(value) for value in stage["primary_contexts"]]
    frame = genes[
        genes["window_width"].eq(int(stage["primary_width"]))
        & genes["context_id"].isin(primary_contexts)
    ].copy()
    frame["case_specific_localization"] = frame["case_localized"] & ~frame[
        "control_localized"
    ]
    index = ["tree", "region", "module_id", "go_invisible", "gene"]
    pivot = frame.pivot(
        index=index, columns="context_id", values="case_specific_localization"
    )
    for context in primary_contexts:
        if context not in pivot:
            pivot[context] = False
    pivot = pivot[primary_contexts].fillna(False).reset_index()
    pivot["primary_contexts_supported"] = (
        pivot[primary_contexts].sum(axis=1).astype(int)
    )
    pivot["replicates_in_both_primary_contexts"] = pivot[
        "primary_contexts_supported"
    ].eq(len(primary_contexts))
    return pivot.sort_values(index, kind="mergesort").reset_index(drop=True)


def _direct_plus_responsive(
    window_calls: pd.DataFrame, ctag_calls: pd.DataFrame, stage: dict[str, Any]
) -> pd.DataFrame:
    mapping = {str(key): str(value) for key, value in stage["ctag_context_map"].items()}
    frames = []
    for perturbation_context, ctag_context in mapping.items():
        perturbation = window_calls[
            window_calls["context_id"].eq(perturbation_context)
        ].copy()
        occupancy = ctag_calls[ctag_calls["context_id"].eq(ctag_context)][
            ["window_id", "callable", "bound", "covered_fraction", "max_coverage"]
        ].rename(
            columns={
                "callable": "ctag_callable",
                "bound": "ctag_bound",
                "covered_fraction": "ctag_covered_fraction",
                "max_coverage": "ctag_max_coverage",
            }
        )
        joined = perturbation.merge(
            occupancy, on="window_id", how="left", validate="one_to_one"
        )
        joined["ctag_context_id"] = ctag_context
        joined["ctag_callable"] = joined["ctag_callable"].fillna(False).astype(bool)
        joined["ctag_bound"] = joined["ctag_bound"].fillna(False).astype(bool)
        joined["direct_plus_responsive"] = (
            joined["localized"] & joined["ctag_callable"] & joined["ctag_bound"]
        )
        frames.append(joined)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _configured_event_files(
    cfg: dict[str, Any], dataset_id: str
) -> dict[str, dict[str, Path]]:
    dataset = next(
        (row for row in cfg["datasets"] if str(row["dataset_id"]) == dataset_id),
        None,
    )
    if dataset is None:
        raise ValueError(f"Dataset not configured: {dataset_id}")
    root = _configured_path(cfg["output_root"])
    paths: dict[str, dict[str, Path]] = {}
    role_map = {
        "significant_event_statistics": "statistics",
        "significant_event_coordinates": "bed12",
    }
    for sample in dataset.get("samples", []):
        context_id = str(sample["context_id"])
        context_paths: dict[str, Path] = {}
        for file_cfg in sample.get("files", []):
            role = role_map.get(str(file_cfg["role"]))
            if role is not None:
                context_paths[role] = root / str(file_cfg["path"])
        if set(context_paths) != {"statistics", "bed12"}:
            raise ValueError(f"{context_id} must configure statistics and BED12 files")
        paths[context_id] = context_paths
    return paths


def _mapping_qc(
    statistics: pd.DataFrame,
    flanks: pd.DataFrame,
    reciprocal: pd.DataFrame,
    window_calls: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    widths = sorted(int(value) for value in flanks["window_width"].unique())
    for context_id in sorted(statistics["context_id"].unique()):
        stat = statistics[statistics["context_id"].eq(context_id)]
        for width in widths:
            event_flanks = flanks[
                flanks["context_id"].eq(context_id)
                & flanks["window_width"].eq(width)
            ]
            mapped = reciprocal[
                reciprocal["context_id"].eq(context_id)
                & reciprocal["window_width"].eq(width)
            ]
            calls = window_calls[
                window_calls["context_id"].eq(context_id)
                & window_calls["window_width"].eq(width)
            ]
            rows.append(
                {
                    "context_id": context_id,
                    "window_width": width,
                    "unique_events": int(len(stat)),
                    "eligible_events": int(stat["passes_event_thresholds"].sum()),
                    "event_flanks": int(len(event_flanks)),
                    "reciprocal_flanks": int(mapped["reciprocal_mapped"].sum()),
                    "reciprocal_fraction": float(mapped["reciprocal_mapped"].mean()),
                    "localized_human_windows": int(calls["localized"].sum()),
                }
            )
    return pd.DataFrame(rows)


def _validate_expected_counts(
    mapping_qc: pd.DataFrame, expected_counts: dict[str, Any]
) -> None:
    for context_id, expected in expected_counts.items():
        observed = mapping_qc[mapping_qc["context_id"].eq(str(context_id))]
        if observed.empty:
            raise RuntimeError(f"Missing expected perturbation context: {context_id}")
        unique_events = int(observed["unique_events"].iloc[0])
        if unique_events != int(expected["unique_events"]):
            raise RuntimeError(
                f"{context_id} unique event count changed: {unique_events}"
            )
        for field in ["event_flanks", "reciprocal_flanks", "localized_windows"]:
            observed_field = (
                "localized_human_windows" if field == "localized_windows" else field
            )
            for width_value, expected_value in expected[field].items():
                width = int(width_value)
                row = observed[observed["window_width"].eq(width)]
                if len(row) != 1 or int(row[observed_field].iloc[0]) != int(
                    expected_value
                ):
                    actual = None if row.empty else int(row[observed_field].iloc[0])
                    raise RuntimeError(
                        f"{context_id} {field} at {width} nt changed: {actual}"
                    )


def _write_report(
    outputs: dict[str, Path],
    candidates: pd.DataFrame,
    mapping_qc: pd.DataFrame,
    overall: pd.DataFrame,
    modules: pd.DataFrame,
    replication: pd.DataFrame,
    direct: pd.DataFrame,
    stage: dict[str, Any],
) -> None:
    primary_width = int(stage["primary_width"])
    display = overall[
        [
            "context_id",
            "window_width",
            "callable_units",
            "case_localized_fraction",
            "control_localized_fraction",
            "informative_discordant",
            "matched_odds_ratio_haldane",
            "ci_low_haldane",
            "ci_high_haldane",
            "exact_pvalue",
            "inference_status",
            "analysis_scope",
            "positive_result",
        ]
    ]
    module_display = modules[modules["window_width"].eq(primary_width)][
        [
            "context_id",
            "module_id",
            "callable_units",
            "informative_discordant",
            "matched_odds_ratio_haldane",
            "exact_pvalue",
            "qvalue",
            "inference_status",
            "analysis_scope",
            "positive_result",
        ]
    ]
    replicated = int(replication["replicates_in_both_primary_contexts"].sum())
    direct_count = (
        int(direct["direct_plus_responsive"].sum()) if not direct.empty else 0
    )
    primary_direct_count = int(
        direct[
            direct["window_width"].eq(primary_width)
            & direct["context_id"].isin(stage["primary_contexts"])
        ]["direct_plus_responsive"].sum()
    )
    lines = [
        "# NOVA2 perturbation validation",
        "",
        "This perturbation tier asks whether significant Nova2-cKO splice-event flanks localize preferentially to the case member of frozen human NOVA_FAMILY matched windows. The public Quantas catalogs contain significant events only, so this is a matched spatial-localization test, not genome-wide enrichment or proof of non-response outside cataloged events.",
        "",
        f"- Frozen candidate transcript pairs: **{len(candidates):,}**",
        f"- Primary splice-flank width: **{primary_width} nt**",
        f"- Genes with case-specific localization in both cortical contexts: **{replicated:,}**",
        f"- Primary-width exact window-level direct-plus-responsive calls: **{primary_direct_count:,}**",
        f"- Direct-plus-responsive calls across all prespecified widths: **{direct_count:,}**",
        "- Discovery candidates remain labeled **NOVA_FAMILY**; NOVA2 is the exact perturbation used only in this orthogonal tier.",
        "- An unlocalized window means absent from the significant-only catalog after reciprocal mapping, not proof of no response.",
        "",
        "## Context and mapping QC",
        "",
        _markdown_table(mapping_qc),
        "",
        "## Gene-level primary and sensitivity effects",
        "",
        _markdown_table(display),
        "",
        "## Module-level primary effects",
        "",
        _markdown_table(module_display),
        "",
    ]
    _atomic_text("\n".join(lines), outputs["report"])


def run(cfg: dict[str, Any], config_path: Path) -> dict[str, Path]:
    stage = cfg["nova2_perturbation_stage"]
    output_dir = _configured_path(stage["output_dir"])
    ensure_dir(output_dir)
    candidate_path = _configured_path(stage["candidate_manifest"])
    candidates = pd.read_parquet(candidate_path)
    if len(candidates) != int(stage["expected_candidate_rows"]):
        raise RuntimeError("NOVA-family candidate row count changed")
    identity = _candidate_identity_sha256(candidates)
    if identity != str(stage["expected_candidate_identity_sha256"]):
        raise RuntimeError(f"NOVA-family candidate identity changed: {identity}")

    human_window_path = _configured_path(stage["human_window_manifest"])
    membership_path = _configured_path(stage["human_window_membership"])
    ctag_path = _configured_path(stage["ctag_overlap_calls"])
    human_windows = pd.read_parquet(human_window_path)
    membership = pd.read_parquet(membership_path)
    ctag_calls = pd.read_parquet(ctag_path)
    membership_candidate_ids = set(membership["candidate_id"].dropna())
    # Membership holds canonical candidate ids, which regions can share.
    expected_unique = stage.get("expected_unique_candidates", stage["expected_candidate_rows"])
    if len(membership_candidate_ids) != int(expected_unique):
        raise RuntimeError("Frozen human-window candidate membership changed")
    if set(human_windows["candidate_id"].dropna()) - membership_candidate_ids:
        raise RuntimeError("Frozen human windows contain unexpected candidates")

    event_files = _configured_event_files(cfg, str(stage["dataset_id"]))
    expected_contexts = set(stage["primary_contexts"]) | set(
        stage["secondary_contexts"]
    )
    if set(event_files) != expected_contexts:
        raise RuntimeError(f"Configured perturbation contexts changed: {sorted(event_files)}")
    widths = sorted(
        {
            int(stage["primary_width"]),
            *(int(value) for value in stage["sensitivity_widths"]),
        }
    )
    statistics_frames = []
    bed_frames = []
    flank_frames = []
    for context_id in sorted(event_files):
        paths = event_files[context_id]
        for path in paths.values():
            if not path.is_file():
                raise FileNotFoundError(path)
        print(f"parsing {context_id}", flush=True)
        statistics = _read_event_statistics(
            paths["statistics"],
            context_id,
            float(stage["event_fdr_max"]),
            float(stage["absolute_delta_inclusion_min"]),
        )
        bed12 = _read_bed12(paths["bed12"], context_id)
        flanks = _event_flanks(statistics, bed12, widths)
        statistics_frames.append(statistics)
        bed_frames.append(bed12)
        flank_frames.append(flanks)
    statistics = pd.concat(statistics_frames, ignore_index=True)
    bed12 = pd.concat(bed_frames, ignore_index=True)
    flanks = pd.concat(flank_frames, ignore_index=True)
    reciprocal = _reciprocal_map(
        flanks,
        _configured_path(stage["forward_chain"]),
        _configured_path(stage["reverse_chain"]),
        float(stage["reciprocal_overlap_min"]),
    )

    call_frames = []
    overlap_frames = []
    for context_id in sorted(expected_contexts):
        print(f"localizing {context_id}", flush=True)
        context_flanks = reciprocal[reciprocal["context_id"].eq(context_id)]
        calls, overlaps = _overlap_event_flanks(
            human_windows, context_flanks, context_id
        )
        call_frames.append(calls)
        overlap_frames.append(overlaps)
    window_calls = pd.concat(call_frames, ignore_index=True)
    nonempty_overlaps = [frame for frame in overlap_frames if not frame.empty]
    event_window_overlaps = (
        pd.concat(nonempty_overlaps, ignore_index=True)
        if nonempty_overlaps
        else overlap_frames[0].iloc[0:0].copy()
    )
    pairs, candidate_support, gene_support = _collapse_support(window_calls, membership)
    overall, modules, window_effects = _effect_tables(
        pairs,
        gene_support,
        int(stage["primary_width"]),
        int(stage["min_informative_discordant_genes"]),
        [str(value) for value in stage["primary_contexts"]],
    )
    replication = _cross_context_replication(gene_support, stage)
    direct = _direct_plus_responsive(window_calls, ctag_calls, stage)
    mapping_qc = _mapping_qc(statistics, flanks, reciprocal, window_calls)
    _validate_expected_counts(mapping_qc, stage["expected_context_counts"])

    outputs = {
        "events": output_dir / "parsed_events.parquet",
        "bed12": output_dir / "parsed_event_coordinates.parquet",
        "event_flanks": output_dir / "mouse_event_flanks.parquet",
        "reciprocal_flanks": output_dir / "reciprocal_event_flanks.parquet",
        "event_window_overlaps": output_dir / "event_window_overlaps.parquet",
        "window_calls": output_dir / "window_localization_calls.parquet",
        "pair_calls": output_dir / "matched_pair_calls.parquet",
        "candidate_support": output_dir / "candidate_support.parquet",
        "gene_support": output_dir / "gene_support.parquet",
        "overall_effects": output_dir / "overall_effects.parquet",
        "module_effects": output_dir / "module_effects.parquet",
        "window_effects": output_dir / "window_effects.parquet",
        "mapping_qc": output_dir / "mapping_qc.parquet",
        "replication": output_dir / "cross_context_replication.parquet",
        "direct": output_dir / "direct_plus_responsive.parquet",
        "json": output_dir / "nova2_perturbation_validation.json",
        "report": output_dir / "NOVA2_PERTURBATION_VALIDATION.md",
    }
    frames = {
        "events": statistics,
        "bed12": bed12,
        "event_flanks": flanks,
        "reciprocal_flanks": reciprocal,
        "event_window_overlaps": event_window_overlaps,
        "window_calls": window_calls,
        "pair_calls": pairs,
        "candidate_support": candidate_support,
        "gene_support": gene_support,
        "overall_effects": overall,
        "module_effects": modules,
        "window_effects": window_effects,
        "mapping_qc": mapping_qc,
        "replication": replication,
        "direct": direct,
    }
    for name, frame in frames.items():
        _atomic_parquet(frame, outputs[name])

    provenance = {
        "config": str(config_path),
        "config_sha256": sha256(config_path),
        "source_sha256": sha256(Path(__file__)),
        "candidate_manifest": str(candidate_path),
        "candidate_manifest_sha256": sha256(candidate_path),
        "candidate_identity_sha256": identity,
        "human_window_manifest_sha256": sha256(human_window_path),
        "human_window_membership_sha256": sha256(membership_path),
        "ctag_overlap_calls_sha256": sha256(ctag_path),
        "forward_chain_sha256": sha256(_configured_path(stage["forward_chain"])),
        "reverse_chain_sha256": sha256(_configured_path(stage["reverse_chain"])),
        "event_source_sha256": {
            context: {role: sha256(path) for role, path in paths.items()}
            for context, paths in event_files.items()
        },
        "stage": stage,
        "observed": {
            "candidate_rows": int(len(candidates)),
            "matched_human_window_rows": int(
                human_windows["match_status"].eq("matched").sum()
            ),
            "unique_event_rows": int(len(statistics)),
            "eligible_event_rows": int(statistics["passes_event_thresholds"].sum()),
            "mouse_event_flank_rows": int(len(flanks)),
            "reciprocal_event_flank_rows": int(reciprocal["reciprocal_mapped"].sum()),
            "localized_human_window_rows": int(window_calls["localized"].sum()),
            "gene_support_rows": int(len(gene_support)),
            "direct_plus_responsive_rows": int(
                direct["direct_plus_responsive"].sum()
            ),
        },
        "output_sha256": {name: sha256(outputs[name]) for name in frames},
    }
    _atomic_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", outputs["json"]
    )
    _write_report(
        outputs,
        candidates,
        mapping_qc,
        overall,
        modules,
        replication,
        direct,
        stage,
    )
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
