"""Prepare deterministic switch-relevant windows for neuronal CLIP validation."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

import numpy as np
import pandas as pd

from isograph_benchmark.config import load_yaml
from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.neuronal_clip_fetch import sha256


DEFAULT_CONFIG = "configs/neuronal_clip.yaml"


@dataclass(frozen=True)
class Flank:
    transcript_id: str
    chrom: str
    start: int
    end: int
    strand: str
    structural_class: str
    intron_index: int
    flank_id: str

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def coordinate_key(self) -> tuple[str, int, int, str]:
        return self.chrom, self.start, self.end, self.strand


@dataclass(frozen=True)
class TranscriptModel:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: tuple[tuple[int, int], ...]


class IndexedFasta:
    def __init__(self, path: Path):
        self.path = path
        self.index_path = Path(f"{path}.fai")
        if not self.path.is_file():
            raise FileNotFoundError(self.path)
        if not self.index_path.is_file():
            raise FileNotFoundError(self.index_path)
        self.index: dict[str, tuple[int, int, int, int]] = {}
        with self.index_path.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                self.index[fields[0]] = tuple(map(int, fields[1:5]))
        self.handle: BinaryIO = self.path.open("rb")

    def close(self) -> None:
        self.handle.close()

    def __enter__(self) -> IndexedFasta:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def fetch(self, chrom: str, start: int, end: int) -> str:
        if chrom not in self.index:
            raise KeyError(f"Chromosome {chrom!r} is absent from {self.index_path}")
        chrom_length, offset, line_bases, line_width = self.index[chrom]
        if start < 0 or end < start or end > chrom_length:
            raise ValueError(
                f"Invalid interval {chrom}:{start}-{end} for length {chrom_length}"
            )
        if start == end:
            return ""
        byte_start = offset + (start // line_bases) * line_width + start % line_bases
        newline_bytes = line_width - line_bases
        lines_spanned = (end - 1) // line_bases - start // line_bases
        bytes_needed = end - start + lines_spanned * newline_bytes
        self.handle.seek(byte_start)
        sequence = self.handle.read(bytes_needed).replace(b"\n", b"").replace(b"\r", b"")
        return sequence[: end - start].decode("ascii").upper()


def _configured_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else rel(value)


def _stable_id(prefix: str, *values: object) -> str:
    payload = "\x1f".join(str(value) for value in values)
    return f"{prefix}{hashlib.sha256(payload.encode()).hexdigest()}"


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


def _markdown_table(frame: pd.DataFrame) -> str:
    headers = [str(column) for column in frame.columns]
    rows = [headers, ["---"] * len(headers)]
    for values in frame.itertuples(index=False, name=None):
        rows.append([str(value) for value in values])
    return "\n".join("| " + " | ".join(row) + " |" for row in rows)


def _canonical_candidates(candidate: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {
        "region",
        "module_id",
        "gene_id",
        "rbp",
        "transcript_id_1",
        "transcript_id_2",
        "motif_count_1",
        "motif_count_2",
    }
    missing = sorted(required - set(candidate.columns))
    if missing:
        raise ValueError(f"Candidate manifest is missing columns: {missing}")

    membership = candidate.reset_index(drop=True).copy()
    membership["candidate_manifest_row"] = np.arange(len(membership), dtype=np.int64)
    positive_1 = membership["motif_count_1"] > 0
    positive_2 = membership["motif_count_2"] > 0
    invalid = positive_1 == positive_2
    if invalid.any():
        raise ValueError(
            f"Candidate manifest contains {int(invalid.sum())} rows without an exact "
            "one-transcript motif-presence switch"
        )
    membership["motif_positive_transcript"] = np.where(
        positive_1, membership["transcript_id_1"], membership["transcript_id_2"]
    )
    membership["motif_negative_transcript"] = np.where(
        positive_1, membership["transcript_id_2"], membership["transcript_id_1"]
    )
    membership["motif_positive_count"] = np.where(
        positive_1, membership["motif_count_1"], membership["motif_count_2"]
    ).astype(int)
    membership["transcript_a"] = membership[
        ["transcript_id_1", "transcript_id_2"]
    ].min(axis=1)
    membership["transcript_b"] = membership[
        ["transcript_id_1", "transcript_id_2"]
    ].max(axis=1)
    membership["candidate_id"] = [
        _stable_id("C", rbp, gene_id, tx_a, tx_b)
        for rbp, gene_id, tx_a, tx_b in zip(
            membership["rbp"],
            membership["gene_id"],
            membership["transcript_a"],
            membership["transcript_b"],
        )
    ]

    invariant = [
        "rbp",
        "gene_id",
        "transcript_a",
        "transcript_b",
        "motif_positive_transcript",
        "motif_negative_transcript",
        "motif_positive_count",
    ]
    inconsistent = (
        membership.groupby("candidate_id")[invariant]
        .nunique(dropna=False)
        .gt(1)
        .any(axis=1)
    )
    if inconsistent.any():
        raise ValueError(
            f"{int(inconsistent.sum())} canonical candidates have inconsistent metadata"
        )

    canonical = (
        membership.sort_values(
            ["rbp", "gene_id", "transcript_a", "transcript_b", "candidate_manifest_row"]
        )
        .drop_duplicates("candidate_id", keep="first")
        [
            [
                "candidate_id",
                "rbp",
                "gene_id",
                "transcript_a",
                "transcript_b",
                "motif_positive_transcript",
                "motif_negative_transcript",
                "motif_positive_count",
            ]
        ]
        .reset_index(drop=True)
    )
    return membership, canonical


def _load_transcript_models(
    gtf_cache: Path, transcript_ids: set[str]
) -> dict[str, TranscriptModel]:
    gtf = pd.read_parquet(
        gtf_cache,
        columns=["transcript_id", "gene_id", "chrom", "strand", "feature", "start", "end"],
    )
    exons = gtf[
        (gtf["feature"] == "exon") & gtf["transcript_id"].isin(transcript_ids)
    ].copy()
    models: dict[str, TranscriptModel] = {}
    for transcript_id, group in exons.groupby("transcript_id", sort=False):
        chroms = group["chrom"].astype(str).unique()
        strands = group["strand"].astype(str).unique()
        genes = group["gene_id"].unique()
        if len(chroms) != 1 or len(strands) != 1 or len(genes) != 1:
            continue
        exon_intervals = tuple(
            sorted(
                (int(row.start) - 1, int(row.end))
                for row in group[["start", "end"]].itertuples(index=False)
            )
        )
        models[str(transcript_id)] = TranscriptModel(
            transcript_id=str(transcript_id),
            gene_id=str(genes[0]),
            chrom=str(chroms[0]),
            strand=str(strands[0]),
            exons=exon_intervals,
        )
    return models


def _intron_flanks(model: TranscriptModel, width: int) -> list[Flank]:
    introns = [
        (left[1], right[0])
        for left, right in zip(model.exons[:-1], model.exons[1:])
        if left[1] < right[0]
    ]
    flanks: list[Flank] = []
    n_introns = len(introns)
    for genomic_index, (start, end) in enumerate(introns, start=1):
        transcript_index = (
            genomic_index if model.strand == "+" else n_introns - genomic_index + 1
        )
        length = end - start
        if length <= 2 * width:
            flank_id = f"{model.transcript_id}:I{transcript_index}:short:{width}"
            flanks.append(
                Flank(
                    model.transcript_id,
                    model.chrom,
                    start,
                    end,
                    model.strand,
                    "short_intron",
                    transcript_index,
                    flank_id,
                )
            )
            continue
        low_class = "donor_flank" if model.strand == "+" else "acceptor_flank"
        high_class = "acceptor_flank" if model.strand == "+" else "donor_flank"
        flanks.extend(
            [
                Flank(
                    model.transcript_id,
                    model.chrom,
                    start,
                    start + width,
                    model.strand,
                    low_class,
                    transcript_index,
                    f"{model.transcript_id}:I{transcript_index}:{low_class}:{width}",
                ),
                Flank(
                    model.transcript_id,
                    model.chrom,
                    end - width,
                    end,
                    model.strand,
                    high_class,
                    transcript_index,
                    f"{model.transcript_id}:I{transcript_index}:{high_class}:{width}",
                ),
            ]
        )
    return flanks


def _overlaps(left: Flank, right: Flank) -> bool:
    return (
        left.chrom == right.chrom
        and left.strand == right.strand
        and left.start < right.end
        and right.start < left.end
    )


def _sequence_metrics(
    fasta: IndexedFasta,
    flank: Flank,
    cache: dict[tuple[str, int, int], tuple[float, float]],
) -> tuple[float, float]:
    key = flank.chrom, flank.start, flank.end
    if key in cache:
        return cache[key]
    sequence = fasta.fetch(*key)
    if not sequence:
        metrics = math.nan, 1.0
    else:
        canonical = sum(sequence.count(base) for base in "ACGT")
        gc = (sequence.count("G") + sequence.count("C")) / canonical if canonical else math.nan
        ambiguous = 1.0 - canonical / len(sequence)
        metrics = gc, ambiguous
    cache[key] = metrics
    return metrics


def _match_flanks(
    cases: list[Flank],
    controls: list[Flank],
    fasta: IndexedFasta,
    stage_cfg: dict[str, Any],
    metric_cache: dict[tuple[str, int, int], tuple[float, float]],
) -> tuple[dict[int, tuple[int, float]], dict[int, str]]:
    length_min = float(stage_cfg["length_ratio_min"])
    length_max = float(stage_cfg["length_ratio_max"])
    gc_tolerance = float(stage_cfg["gc_tolerance"])
    max_ambiguous = float(stage_cfg["max_ambiguous_fraction"])
    options: dict[int, list[tuple[int, float]]] = {}
    exclusions: dict[int, str] = {}

    for case_index, case in enumerate(cases):
        case_gc, case_ambiguous = _sequence_metrics(fasta, case, metric_cache)
        if not math.isfinite(case_gc) or case_ambiguous > max_ambiguous:
            options[case_index] = []
            exclusions[case_index] = "case_sequence_not_callable"
            continue
        viable: list[tuple[int, float]] = []
        for control_index, control in enumerate(controls):
            if case.structural_class != control.structural_class or _overlaps(case, control):
                continue
            ratio = control.length / case.length
            if ratio < length_min or ratio > length_max:
                continue
            control_gc, control_ambiguous = _sequence_metrics(fasta, control, metric_cache)
            if not math.isfinite(control_gc) or control_ambiguous > max_ambiguous:
                continue
            gc_difference = abs(control_gc - case_gc)
            if gc_difference > gc_tolerance:
                continue
            score = abs(math.log(ratio)) + gc_difference
            viable.append((control_index, score))
        options[case_index] = sorted(
            viable,
            key=lambda item: (
                item[1],
                controls[item[0]].chrom,
                controls[item[0]].start,
                controls[item[0]].end,
                controls[item[0]].flank_id,
            ),
        )
        if not viable:
            exclusions[case_index] = "no_strict_within_gene_control"

    used_controls: set[int] = set()
    matches: dict[int, tuple[int, float]] = {}
    case_order = sorted(
        range(len(cases)),
        key=lambda index: (
            len(options[index]),
            cases[index].chrom,
            cases[index].start,
            cases[index].end,
            cases[index].flank_id,
        ),
    )
    for case_index in case_order:
        available = [
            item for item in options[case_index] if item[0] not in used_controls
        ]
        if not available:
            if options[case_index]:
                exclusions[case_index] = "control_exhausted_without_replacement"
            continue
        control_index, score = available[0]
        matches[case_index] = control_index, score
        used_controls.add(control_index)
        exclusions.pop(case_index, None)
    return matches, exclusions


def _window_row(
    candidate: Any,
    flank: Flank,
    role: str,
    width: int,
    match_id: str | None,
    match_status: str,
    exclusion_reason: str,
    match_score: float | None,
    fasta: IndexedFasta,
    stage_cfg: dict[str, Any],
    metric_cache: dict[tuple[str, int, int], tuple[float, float]],
) -> dict[str, Any]:
    gc, ambiguous = _sequence_metrics(fasta, flank, metric_cache)
    return {
        "window_id": _stable_id(
            "W",
            candidate.candidate_id,
            width,
            role,
            flank.chrom,
            flank.start,
            flank.end,
            flank.strand,
            flank.flank_id,
        ),
        "candidate_id": candidate.candidate_id,
        "rbp": candidate.rbp,
        "gene_id": candidate.gene_id,
        "transcript_a": candidate.transcript_a,
        "transcript_b": candidate.transcript_b,
        "motif_positive_transcript": candidate.motif_positive_transcript,
        "motif_negative_transcript": candidate.motif_negative_transcript,
        "motif_positive_count": int(candidate.motif_positive_count),
        "window_width": width,
        "is_primary_width": width == int(stage_cfg["primary_width"]),
        "match_id": match_id or "",
        "window_role": role,
        "match_status": match_status,
        "exclusion_reason": exclusion_reason,
        "chrom": flank.chrom,
        "start": flank.start,
        "end": flank.end,
        "strand": flank.strand,
        "coordinate_system": "0-based_half-open",
        "source_transcript": flank.transcript_id,
        "flank_id": flank.flank_id,
        "structural_class": flank.structural_class,
        "event_class": "internal_splice_flank",
        "intron_index": flank.intron_index,
        "length": flank.length,
        "callable_length": (
            flank.length
            if ambiguous <= float(stage_cfg["max_ambiguous_fraction"])
            else 0
        ),
        "gc_fraction": gc,
        "ambiguous_fraction": ambiguous,
        "mappability_mean": np.nan,
        "mappability_status": stage_cfg["mappability_policy"],
        "input_coverage": np.nan,
        "input_coverage_status": stage_cfg["input_coverage_policy"],
        "orthology_status": "not_applicable_human",
        "match_score": match_score,
    }


def _candidate_windows(
    candidate: Any,
    models: dict[str, TranscriptModel],
    widths: list[int],
    fasta: IndexedFasta,
    stage_cfg: dict[str, Any],
    metric_cache: dict[tuple[str, int, int], tuple[float, float]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    positive = models.get(candidate.motif_positive_transcript)
    negative = models.get(candidate.motif_negative_transcript)
    rows: list[dict[str, Any]] = []
    statuses: list[dict[str, Any]] = []

    base_status = {
        "candidate_id": candidate.candidate_id,
        "rbp": candidate.rbp,
        "gene_id": candidate.gene_id,
        "transcript_a": candidate.transcript_a,
        "transcript_b": candidate.transcript_b,
        "motif_positive_transcript": candidate.motif_positive_transcript,
        "motif_negative_transcript": candidate.motif_negative_transcript,
    }
    if positive is None or negative is None:
        reason = (
            "missing_positive_transcript"
            if positive is None
            else "missing_negative_transcript"
        )
        for width in widths:
            statuses.append(
                {
                    **base_status,
                    "window_width": width,
                    "n_positive_flanks": 0,
                    "n_control_flanks": 0,
                    "n_case_windows": 0,
                    "n_matched_sets": 0,
                    "n_unmatched_cases": 0,
                    "window_status": reason,
                }
            )
        return rows, statuses
    if positive.gene_id != candidate.gene_id or negative.gene_id != candidate.gene_id:
        for width in widths:
            statuses.append(
                {
                    **base_status,
                    "window_width": width,
                    "n_positive_flanks": 0,
                    "n_control_flanks": 0,
                    "n_case_windows": 0,
                    "n_matched_sets": 0,
                    "n_unmatched_cases": 0,
                    "window_status": "gtf_gene_mismatch",
                }
            )
        return rows, statuses
    if positive.chrom != negative.chrom or positive.strand != negative.strand:
        for width in widths:
            statuses.append(
                {
                    **base_status,
                    "window_width": width,
                    "n_positive_flanks": 0,
                    "n_control_flanks": 0,
                    "n_case_windows": 0,
                    "n_matched_sets": 0,
                    "n_unmatched_cases": 0,
                    "window_status": "transcript_coordinate_mismatch",
                }
            )
        return rows, statuses

    for width in widths:
        positive_flanks = _intron_flanks(positive, width)
        negative_flanks = _intron_flanks(negative, width)
        negative_keys = {flank.coordinate_key for flank in negative_flanks}
        cases = [
            flank for flank in positive_flanks if flank.coordinate_key not in negative_keys
        ]
        controls = [
            flank
            for flank in negative_flanks
            if all(not _overlaps(flank, case) for case in cases)
        ]
        if not cases:
            statuses.append(
                {
                    **base_status,
                    "window_width": width,
                    "n_positive_flanks": len(positive_flanks),
                    "n_control_flanks": len(controls),
                    "n_case_windows": 0,
                    "n_matched_sets": 0,
                    "n_unmatched_cases": 0,
                    "window_status": "no_pair_specific_positive_flank",
                }
            )
            continue

        matches, exclusions = _match_flanks(
            cases, controls, fasta, stage_cfg, metric_cache
        )
        for case_index, case in enumerate(cases):
            matched = case_index in matches
            if matched:
                control_index, score = matches[case_index]
                match_id = _stable_id("M", candidate.candidate_id, width, case.flank_id)
                rows.append(
                    _window_row(
                        candidate,
                        case,
                        "case",
                        width,
                        match_id,
                        "matched",
                        "",
                        score,
                        fasta,
                        stage_cfg,
                        metric_cache,
                    )
                )
                rows.append(
                    _window_row(
                        candidate,
                        controls[control_index],
                        "control",
                        width,
                        match_id,
                        "matched",
                        "",
                        score,
                        fasta,
                        stage_cfg,
                        metric_cache,
                    )
                )
            else:
                rows.append(
                    _window_row(
                        candidate,
                        case,
                        "case",
                        width,
                        None,
                        "unmatched",
                        exclusions[case_index],
                        None,
                        fasta,
                        stage_cfg,
                        metric_cache,
                    )
                )
        statuses.append(
            {
                **base_status,
                "window_width": width,
                "n_positive_flanks": len(positive_flanks),
                "n_control_flanks": len(controls),
                "n_case_windows": len(cases),
                "n_matched_sets": len(matches),
                "n_unmatched_cases": len(cases) - len(matches),
                "window_status": (
                    "preclip_matched" if matches else "no_strict_within_gene_match"
                ),
            }
        )
    return rows, statuses


def _write_report(
    output_dir: Path,
    membership: pd.DataFrame,
    canonical: pd.DataFrame,
    windows: pd.DataFrame,
    statuses: pd.DataFrame,
    stage_cfg: dict[str, Any],
    provenance: dict[str, Any],
) -> None:
    primary_width = int(stage_cfg["primary_width"])
    primary = statuses[statuses["window_width"] == primary_width]
    status_counts = (
        statuses.groupby(["window_width", "window_status"], as_index=False)
        .size()
        .sort_values(["window_width", "window_status"])
    )
    rbp_summary = (
        primary.groupby("rbp", as_index=False)
        .agg(
            candidates=("candidate_id", "size"),
            preclip_matched=(
                "window_status",
                lambda values: int((values == "preclip_matched").sum()),
            ),
            matched_sets=("n_matched_sets", "sum"),
        )
        .sort_values("rbp")
    )
    matched = windows[windows["match_status"] == "matched"]
    summary = {
        "candidate_manifest_rows": int(len(membership)),
        "unique_candidates": int(len(canonical)),
        "window_rows": int(len(windows)),
        "primary_width": primary_width,
        "primary_preclip_matched_candidates": int(
            (primary["window_status"] == "preclip_matched").sum()
        ),
        "primary_matched_sets": int(primary["n_matched_sets"].sum()),
        "matched_rows": int(len(matched)),
        "mappability_policy": stage_cfg["mappability_policy"],
        "input_coverage_policy": stage_cfg["input_coverage_policy"],
        "status_counts": status_counts.to_dict(orient="records"),
        "rbp_summary": rbp_summary.to_dict(orient="records"),
        "provenance": provenance,
    }
    _atomic_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        output_dir / "window_qc.json",
    )

    lines = [
        "# Neuronal CLIP switch-window QC",
        "",
        f"- Candidate-manifest rows accounted for: **{len(membership):,}**",
        f"- Canonical RBP/gene/transcript-pair candidates: **{len(canonical):,}**",
        f"- Window rows written: **{len(windows):,}**",
        f"- Primary width: **{primary_width} nt**",
        f"- Primary candidates with at least one strict pre-CLIP match: "
        f"**{summary['primary_preclip_matched_candidates']:,}**",
        f"- Primary matched sets: **{summary['primary_matched_sets']:,}**",
        "",
        "Mappability and context-specific input coverage are deliberately deferred to the "
        "callable-window stage. IP peaks were not read during coordinate construction or matching.",
        "",
        "## Primary-width summary by RBP",
        "",
        _markdown_table(rbp_summary),
        "",
        "## Status accounting",
        "",
        _markdown_table(status_counts),
        "",
    ]
    _atomic_text("\n".join(lines), output_dir / "WINDOW_QC.md")


def prepare_windows(
    cfg: dict[str, Any],
    config_path: Path,
    output_dir: Path,
    selected_rbps: set[str],
    max_candidates: int | None,
) -> dict[str, Path]:
    stage_cfg = cfg["window_stage"]
    if not bool(stage_cfg["match_without_replacement"]):
        raise ValueError("The prespecified window stage requires matching without replacement")
    candidate_path = rel(cfg["candidate_manifest"])
    gtf_cache = _configured_path(stage_cfg["gtf_cache"])
    genome_fasta = _configured_path(stage_cfg["genome_fasta"])
    candidate = pd.read_parquet(candidate_path)
    full_run = not selected_rbps and max_candidates is None

    if full_run and len(candidate) != int(stage_cfg["expected_candidate_rows"]):
        raise RuntimeError(
            f"Candidate row count changed: {len(candidate)} observed, "
            f"{stage_cfg['expected_candidate_rows']} expected"
        )
    if selected_rbps:
        candidate = candidate[candidate["rbp"].isin(selected_rbps)].copy()
    membership, canonical = _canonical_candidates(candidate)
    if max_candidates is not None:
        canonical = canonical.head(max_candidates).copy()
        selected_ids = set(canonical["candidate_id"])
        membership = membership[membership["candidate_id"].isin(selected_ids)].copy()
    if full_run and len(canonical) != int(stage_cfg["expected_unique_candidates"]):
        raise RuntimeError(
            f"Canonical unordered RBP/gene/transcript-pair count changed: "
            f"{len(canonical)} observed, "
            f"{stage_cfg['expected_unique_candidates']} expected"
        )
    if canonical.empty:
        raise ValueError("No candidates selected")

    ensure_dir(output_dir)
    print(f"output directory: {output_dir}", flush=True)
    transcript_ids = set(canonical["motif_positive_transcript"]) | set(
        canonical["motif_negative_transcript"]
    )
    models = _load_transcript_models(gtf_cache, transcript_ids)
    widths = sorted(
        {
            int(stage_cfg["primary_width"]),
            *(int(width) for width in stage_cfg["sensitivity_widths"]),
        }
    )
    window_rows: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []
    metric_cache: dict[tuple[str, int, int], tuple[float, float]] = {}
    with IndexedFasta(genome_fasta) as fasta:
        for index, row in enumerate(canonical.itertuples(index=False), start=1):
            rows, statuses = _candidate_windows(
                row, models, widths, fasta, stage_cfg, metric_cache
            )
            window_rows.extend(rows)
            status_rows.extend(statuses)
            if index % 1000 == 0:
                print(f"prepared {index:,}/{len(canonical):,} candidates", flush=True)

    windows = pd.DataFrame(window_rows)
    statuses = pd.DataFrame(status_rows)
    primary = statuses[statuses["window_width"] == int(stage_cfg["primary_width"])][
        ["candidate_id", "window_status", "n_case_windows", "n_matched_sets", "n_unmatched_cases"]
    ].rename(
        columns={
            "window_status": "primary_window_status",
            "n_case_windows": "primary_case_windows",
            "n_matched_sets": "primary_matched_sets",
            "n_unmatched_cases": "primary_unmatched_cases",
        }
    )
    membership = membership.merge(primary, on="candidate_id", how="left", validate="many_to_one")
    membership["candidate_manifest_sha256"] = sha256(candidate_path)
    membership["window_config_sha256"] = sha256(config_path)

    membership = membership.sort_values("candidate_manifest_row").reset_index(drop=True)
    canonical = canonical.sort_values(["rbp", "gene_id", "transcript_a", "transcript_b"])
    statuses = statuses.sort_values(["window_width", "rbp", "gene_id", "candidate_id"])
    if not windows.empty:
        windows = windows.sort_values(
            [
                "window_width",
                "rbp",
                "gene_id",
                "candidate_id",
                "match_id",
                "window_role",
                "window_id",
            ]
        ).reset_index(drop=True)

    outputs = {
        "membership": output_dir / "candidate_window_membership.parquet",
        "canonical": output_dir / "canonical_candidates.parquet",
        "status": output_dir / "candidate_window_status.parquet",
        "windows": output_dir / "window_manifest.parquet",
    }
    _atomic_parquet(membership, outputs["membership"])
    _atomic_parquet(canonical, outputs["canonical"])
    _atomic_parquet(statuses, outputs["status"])
    _atomic_parquet(windows, outputs["windows"])

    provenance = {
        "candidate_manifest": str(candidate_path),
        "candidate_manifest_sha256": sha256(candidate_path),
        "config": str(config_path),
        "config_sha256": sha256(config_path),
        "gtf_cache": str(gtf_cache),
        "gtf_cache_sha256": sha256(gtf_cache),
        "genome_fasta": str(genome_fasta),
        "genome_fasta_bytes": genome_fasta.stat().st_size,
        "genome_fai_sha256": sha256(Path(f"{genome_fasta}.fai")),
        "seed": int(cfg["seed"]),
        "window_stage": stage_cfg,
        "partial_run": not full_run,
        "selected_rbps": sorted(selected_rbps),
        "max_candidates": max_candidates,
        "output_sha256": {name: sha256(path) for name, path in outputs.items()},
    }
    _write_report(output_dir, membership, canonical, windows, statuses, stage_cfg, provenance)
    outputs["json"] = output_dir / "window_qc.json"
    outputs["report"] = output_dir / "WINDOW_QC.md"
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser(
        "prepare-windows",
        help="build pair-specific intronic windows and deterministic within-gene matches",
    )
    prepare.add_argument(
        "--output-dir",
        help="override configured output directory; required for partial runs",
    )
    prepare.add_argument(
        "--rbp",
        action="append",
        default=[],
        help="restrict to an RBP for a partial run; repeat as needed",
    )
    prepare.add_argument(
        "--max-candidates",
        type=int,
        help="restrict canonical candidates for a smoke test",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = _configured_path(args.config)
    cfg = load_yaml(str(config_path))
    stage_cfg = cfg["window_stage"]
    partial = bool(args.rbp) or args.max_candidates is not None
    if partial and not args.output_dir:
        raise ValueError("Partial runs require --output-dir to protect official outputs")
    output_dir = (
        _configured_path(args.output_dir)
        if args.output_dir
        else _configured_path(stage_cfg["output_dir"])
    )
    outputs = prepare_windows(
        cfg,
        config_path,
        output_dir,
        set(args.rbp),
        args.max_candidates,
    )
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
