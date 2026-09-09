"""Audit motif opportunity and quarantine invalid neuronal CLIP nominations."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from isograph_benchmark.config import load_yaml
from isograph_benchmark.paths import cohort_dir, ensure_dir
from isograph_benchmark.real_data.neuronal_clip_fetch import sha256
from isograph_benchmark.real_data.neuronal_clip_validation import (
    DEFAULT_CONFIG,
    Flank,
    IndexedFasta,
    TranscriptModel,
    _configured_path,
    _intron_flanks,
    _load_transcript_models,
    _stable_id,
)
from isograph_benchmark.real_data.rbp_regulon import _REGIONS


_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


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


def _unique_flanks(model: TranscriptModel, width: int) -> list[Flank]:
    unique: dict[tuple[str, int, int, str], Flank] = {}
    for flank in _intron_flanks(model, width):
        unique.setdefault(flank.coordinate_key, flank)
    return sorted(
        unique.values(),
        key=lambda item: (item.chrom, item.start, item.end, item.structural_class),
    )


def _sense_sequence(sequence: str, strand: str) -> str:
    return sequence if strand == "+" else sequence.translate(_COMPLEMENT)[::-1]


def _scan_transcript(
    model: TranscriptModel,
    width: int,
    motifs: set[str],
    fasta: IndexedFasta,
) -> tuple[list[dict[str, Any]], int]:
    rows: list[dict[str, Any]] = []
    opportunity = 0
    motif_length = {len(motif) for motif in motifs}
    if len(motif_length) != 1:
        raise ValueError("NOVA-family motifs must have a common length")
    k = motif_length.pop()
    for flank in _unique_flanks(model, width):
        genomic = fasta.fetch(flank.chrom, flank.start, flank.end)
        sense = _sense_sequence(genomic, flank.strand)
        opportunity += len(sense)
        for offset in range(max(0, len(sense) - k + 1)):
            motif = sense[offset : offset + k]
            if motif not in motifs:
                continue
            if flank.strand == "+":
                start = flank.start + offset
                end = start + k
            else:
                start = flank.end - offset - k
                end = flank.end - offset
            rows.append(
                {
                    "instance_id": _stable_id(
                        "NI", model.transcript_id, flank.flank_id, start, end, motif
                    ),
                    "transcript_id": model.transcript_id,
                    "gene_id": model.gene_id,
                    "chrom": flank.chrom,
                    "start": start,
                    "end": end,
                    "strand": flank.strand,
                    "coordinate_system": "0-based_half-open",
                    "motif": motif,
                    "motif_family": "NOVA_FAMILY_YCAY",
                    "flank_id": flank.flank_id,
                    "structural_class": flank.structural_class,
                    "intron_index": flank.intron_index,
                    "sense_offset": offset,
                    "flank_length": flank.length,
                }
            )
    return rows, opportunity


def _cluster_instances(
    instances: pd.DataFrame,
    min_instances: int,
    max_span: int,
) -> pd.DataFrame:
    columns = [
        "cluster_id",
        "transcript_id",
        "gene_id",
        "chrom",
        "start",
        "end",
        "strand",
        "coordinate_system",
        "flank_id",
        "structural_class",
        "intron_index",
        "n_instances",
        "span_nt",
        "motifs",
        "cluster_rule",
        "evidence_scope",
    ]
    if instances.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, Any]] = []
    for (_, flank_id), group in instances.groupby(
        ["transcript_id", "flank_id"], sort=True
    ):
        group = group.sort_values(["sense_offset", "start", "motif"]).reset_index(drop=True)
        index = 0
        while index < len(group):
            final = index
            while final + 1 < len(group):
                span = (
                    int(group.loc[final + 1, "sense_offset"])
                    + len(str(group.loc[final + 1, "motif"]))
                    - int(group.loc[index, "sense_offset"])
                )
                if span > max_span:
                    break
                final += 1
            count = final - index + 1
            if count < min_instances:
                index += 1
                continue
            cluster = group.iloc[index : final + 1]
            start = int(cluster["start"].min())
            end = int(cluster["end"].max())
            first = cluster.iloc[0]
            sense_span = (
                int(cluster.iloc[-1]["sense_offset"])
                + len(str(cluster.iloc[-1]["motif"]))
                - int(first["sense_offset"])
            )
            rows.append(
                {
                    "cluster_id": _stable_id(
                        "NC", first.transcript_id, flank_id, start, end, count
                    ),
                    "transcript_id": first.transcript_id,
                    "gene_id": first.gene_id,
                    "chrom": first.chrom,
                    "start": start,
                    "end": end,
                    "strand": first.strand,
                    "coordinate_system": "0-based_half-open",
                    "flank_id": flank_id,
                    "structural_class": first.structural_class,
                    "intron_index": int(first.intron_index),
                    "n_instances": count,
                    "span_nt": sense_span,
                    "motifs": ":".join(cluster["motif"].astype(str)),
                    "cluster_rule": f">={min_instances}_YCAY_within_{max_span}nt",
                    "evidence_scope": "exploratory_qc_only",
                }
            )
            index = final + 1
    return pd.DataFrame(rows, columns=columns)


def _global_structural_audit(counts_path: Path) -> dict[str, Any]:
    counts = pd.read_parquet(counts_path, columns=["transcript_id", "rbp"])
    any_intronic = set(counts["transcript_id"].astype(str))
    nova = set(counts.loc[counts["rbp"].eq("NOVA1"), "transcript_id"].astype(str))
    pair_total = 0
    pair_agree = 0
    nova_switches = 0
    opportunity_switches = 0
    for tree, region in _REGIONS:
        path = cohort_dir(
            tree,
            region,
            "_m",
            "isograph_vae",
            "module_interpret",
            "structure_switch_pairs.parquet",
        )
        if not path.is_file():
            continue
        pairs = pd.read_parquet(
            path, columns=["transcript_id_1", "transcript_id_2"]
        )
        for transcript_1, transcript_2 in pairs.itertuples(index=False, name=None):
            nova_switch = (transcript_1 in nova) != (transcript_2 in nova)
            opportunity_switch = (transcript_1 in any_intronic) != (
                transcript_2 in any_intronic
            )
            pair_total += 1
            nova_switches += int(nova_switch)
            opportunity_switches += int(opportunity_switch)
            pair_agree += int(nova_switch == opportunity_switch)
    return {
        "transcripts_with_any_intronic_motif": len(any_intronic),
        "nova1_positive_transcripts": len(nova),
        "nova1_positive_fraction": len(nova) / len(any_intronic),
        "switch_pairs": pair_total,
        "nova1_switch_pairs": nova_switches,
        "intronic_opportunity_switch_pairs": opportunity_switches,
        "classification_agreement_pairs": pair_agree,
        "classification_agreement_fraction": pair_agree / pair_total,
    }


def run_audit(cfg: dict[str, Any], config_path: Path) -> dict[str, Path]:
    stage = cfg["motif_opportunity_stage"]
    window_stage = cfg["window_stage"]
    output_dir = _configured_path(stage["output_dir"])
    window_dir = _configured_path(stage["window_dir"])
    canonical_path = window_dir / "canonical_candidates.parquet"
    status_path = window_dir / "candidate_window_status.parquet"
    canonical = pd.read_parquet(canonical_path)
    statuses = pd.read_parquet(status_path)
    primary_width = int(stage["primary_width"])
    primary = statuses.loc[
        statuses["window_width"].eq(primary_width),
        ["candidate_id", "window_status", "n_matched_sets"],
    ].rename(
        columns={
            "window_status": "primary_window_status",
            "n_matched_sets": "primary_matched_sets",
        }
    )
    transcript_ids = set(canonical["motif_positive_transcript"].astype(str)) | set(
        canonical["motif_negative_transcript"].astype(str)
    )
    models = _load_transcript_models(
        _configured_path(window_stage["gtf_cache"]), transcript_ids
    )
    nova_cfg = stage["nova_family"]
    motifs = {str(motif).upper() for motif in nova_cfg["motifs_dna_sense"]}
    nova_transcripts = set(
        canonical.loc[canonical["rbp"].eq("NOVA1"), "motif_positive_transcript"]
    ) | set(canonical.loc[canonical["rbp"].eq("NOVA1"), "motif_negative_transcript"])
    instance_rows: list[dict[str, Any]] = []
    fasta_path = _configured_path(window_stage["genome_fasta"])
    with IndexedFasta(fasta_path) as fasta:
        for transcript_id in sorted(nova_transcripts):
            model = models.get(str(transcript_id))
            if model is None:
                continue
            rows, _ = _scan_transcript(model, primary_width, motifs, fasta)
            instance_rows.extend(rows)
    instances = pd.DataFrame(instance_rows)
    if instances.empty:
        instances = pd.DataFrame(
            columns=[
                "instance_id",
                "transcript_id",
                "gene_id",
                "chrom",
                "start",
                "end",
                "strand",
                "coordinate_system",
                "motif",
                "motif_family",
                "flank_id",
                "structural_class",
                "intron_index",
                "sense_offset",
                "flank_length",
            ]
        )
    else:
        instances = instances.drop_duplicates("instance_id").sort_values(
            ["transcript_id", "flank_id", "sense_offset", "motif"]
        )
    clusters = _cluster_instances(
        instances,
        int(nova_cfg["min_instances_per_cluster"]),
        int(nova_cfg["max_cluster_span_nt"]),
    )
    cluster_counts = (
        clusters.groupby("transcript_id").size().to_dict()
        if not clusters.empty
        else {}
    )
    quarantined = stage["quarantined_rbps"]
    rows = []
    for candidate in canonical.itertuples(index=False):
        positive = models.get(candidate.motif_positive_transcript)
        negative = models.get(candidate.motif_negative_transcript)
        positive_opportunity = (
            sum(flank.length for flank in _unique_flanks(positive, primary_width))
            if positive is not None
            else 0
        )
        negative_opportunity = (
            sum(flank.length for flank in _unique_flanks(negative, primary_width))
            if negative is not None
            else 0
        )
        quarantine = quarantined.get(candidate.rbp)
        if quarantine:
            eligibility = str(quarantine["status"])
            reason = str(quarantine["reason"])
            evidence_label = str(quarantine["replacement_label"])
        else:
            eligibility = "pending_window_qc"
            reason = ""
            evidence_label = candidate.rbp
        rows.append(
            {
                "candidate_id": candidate.candidate_id,
                "rbp": candidate.rbp,
                "evidence_label": evidence_label,
                "gene_id": candidate.gene_id,
                "motif_positive_transcript": candidate.motif_positive_transcript,
                "motif_negative_transcript": candidate.motif_negative_transcript,
                "positive_exons": len(positive.exons) if positive else 0,
                "negative_exons": len(negative.exons) if negative else 0,
                "positive_intronic_opportunity_nt": positive_opportunity,
                "negative_intronic_opportunity_nt": negative_opportunity,
                "both_transcripts_scannable": bool(
                    positive_opportunity >= int(nova_cfg["min_intronic_opportunity_nt"])
                    and negative_opportunity
                    >= int(nova_cfg["min_intronic_opportunity_nt"])
                ),
                "positive_nova_family_clusters": (
                    int(cluster_counts.get(candidate.motif_positive_transcript, 0))
                    if candidate.rbp == "NOVA1"
                    else np.nan
                ),
                "negative_nova_family_clusters": (
                    int(cluster_counts.get(candidate.motif_negative_transcript, 0))
                    if candidate.rbp == "NOVA1"
                    else np.nan
                ),
                "upstream_eligibility": eligibility,
                "eligibility_reason": reason,
                "confirmatory_eligible": False,
            }
        )
    eligibility = pd.DataFrame(rows).merge(
        primary, on="candidate_id", how="left", validate="one_to_one"
    )
    eligible_upstream = eligibility["upstream_eligibility"].eq("pending_window_qc")
    eligible_windows = eligibility["primary_window_status"].eq("preclip_matched")
    eligibility.loc[eligible_upstream & eligible_windows, "upstream_eligibility"] = (
        "eligible_exact_rbp_nomination"
    )
    eligibility.loc[eligible_upstream & ~eligible_windows, "upstream_eligibility"] = (
        "not_testable_no_strict_window_match"
    )
    eligibility.loc[
        eligible_upstream & ~eligible_windows, "eligibility_reason"
    ] = eligibility.loc[eligible_upstream & ~eligible_windows, "primary_window_status"]
    eligibility["confirmatory_eligible"] = (
        eligibility["upstream_eligibility"].eq("eligible_exact_rbp_nomination")
        & eligibility["primary_matched_sets"].fillna(0).gt(0)
    )
    global_audit = _global_structural_audit(
        _configured_path(cfg["candidate_source"]["motif_counts"])
    )
    nova_eligibility = eligibility[eligibility["rbp"].eq("NOVA1")]
    summary = {
        "config": str(config_path),
        "config_sha256": sha256(config_path),
        "canonical_candidates_sha256": sha256(canonical_path),
        "primary_width": primary_width,
        "quarantined_rbps": quarantined,
        "global_structural_audit": global_audit,
        "candidates": int(len(eligibility)),
        "confirmatory_eligible_candidates": int(eligibility["confirmatory_eligible"].sum()),
        "nova1_candidates": int(len(nova_eligibility)),
        "nova1_negative_monoexonic": int(nova_eligibility["negative_exons"].eq(1).sum()),
        "nova1_both_transcripts_scannable": int(
            nova_eligibility["both_transcripts_scannable"].sum()
        ),
        "nova_family_instances": int(len(instances)),
        "nova_family_clusters": int(len(clusters)),
        "nova_family_evidence_scope": str(nova_cfg["evidence_scope"]),
    }
    ensure_dir(output_dir)
    outputs = {
        "eligibility": output_dir / "candidate_eligibility.parquet",
        "instances": output_dir / "nova_family_motif_instances.parquet",
        "clusters": output_dir / "nova_family_motif_clusters.parquet",
        "json": output_dir / "motif_opportunity_qc.json",
        "report": output_dir / "MOTIF_OPPORTUNITY_QC.md",
    }
    _atomic_parquet(
        eligibility.sort_values(["rbp", "gene_id", "candidate_id"]),
        outputs["eligibility"],
    )
    _atomic_parquet(instances, outputs["instances"])
    _atomic_parquet(clusters, outputs["clusters"])
    _atomic_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", outputs["json"])
    audit = global_audit
    lines = [
        "# Neuronal CLIP motif-opportunity QC",
        "",
        "## Gate decision",
        "",
        "NOVA1 is `not_testable_current_candidate_definition` and is excluded from confirmatory overlap. Frozen nominations remain traceable.",
        "",
        f"- NOVA1-positive transcripts among transcripts with any intronic motif: **{audit['nova1_positive_transcripts']:,}/{audit['transcripts_with_any_intronic_motif']:,} ({audit['nova1_positive_fraction']:.3%})**",
        f"- Agreement between NOVA1 switching and intronic-opportunity switching: **{audit['classification_agreement_pairs']:,}/{audit['switch_pairs']:,} ({audit['classification_agreement_fraction']:.3%})**",
        f"- Frozen NOVA1 candidates: **{len(nova_eligibility):,}**",
        f"- NOVA1 candidates with a monoexonic motif-negative transcript: **{summary['nova1_negative_monoexonic']:,}**",
        f"- NOVA1 candidates with both transcripts scannable: **{summary['nova1_both_transcripts_scannable']:,}**",
        "",
        "## Coordinate audit",
        "",
        f"- Strand-aware YCAY instances: **{len(instances):,}**",
        f"- Exploratory YCAY clusters: **{len(clusters):,}**",
        f"- Rule: at least {nova_cfg['min_instances_per_cluster']} YCAY instances within {nova_cfg['max_cluster_span_nt']} nt",
        "",
        "The YCAY coordinates are exploratory QC only. They do not rescue or replace the invalidated NOVA1 module nominations, and motif evidence alone is labeled NOVA-family rather than NOVA1-specific.",
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
    config_path = _configured_path(args.config)
    outputs = run_audit(load_yaml(config_path), config_path)
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
