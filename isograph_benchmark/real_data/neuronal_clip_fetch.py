"""Reproducible acquisition and freezing of neuronal CLIP validation inputs."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import tarfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from isograph_benchmark.config import load_yaml
from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.qtl_anchoring import _bare
from isograph_benchmark.real_data.rbp_regulon import _REGIONS

DEFAULT_CONFIG = "configs/neuronal_clip.yaml"
CHUNK_SIZE = 1024 * 1024
CANDIDATE_IDENTITY_COLUMNS = [
    "region",
    "module_id",
    "rbp",
    "gene",
    "transcript_id_1",
    "transcript_id_2",
]
RETRIEVAL_FIELDS = (
    "retrieved_at_utc",
    "last_modified",
    "retrieval_metadata_source",
)


def sha256(path: Path, chunk_size: int = CHUNK_SIZE) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_hashes(cfg: dict[str, Any]) -> dict[str, str]:
    source = cfg["candidate_source"]
    return {
        "regulon_sha256": sha256(rel(source["regulon"])),
        "switch_calls_sha256": sha256(rel(source["switch_calls"])),
        "motif_counts_sha256": sha256(rel(source["motif_counts"])),
    }


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


def _guard_candidate_identity(source: dict[str, Any], observed: str) -> None:
    configured = source.get("expected_candidate_identity_sha256")
    if configured is None:
        print(
            "candidate identity unpinned; observed sha256 "
            f"{observed} -- pin it as "
            "candidate_source.expected_candidate_identity_sha256"
        )
        return
    if str(configured) != observed:
        raise RuntimeError(
            "Neuronal-CLIP candidate identity changed: "
            f"{observed} observed, {configured} expected"
        )


def freeze_candidates(cfg: dict[str, Any]) -> Path:
    source = cfg["candidate_source"]
    rbps = set(source["rbps"])
    q_max = float(source["q_max"])
    regulon_path = rel(source["regulon"])
    calls_path = rel(source["switch_calls"])
    counts_path = rel(source["motif_counts"])

    regulon = pd.read_parquet(regulon_path)
    regulon = regulon[regulon["rbp"].isin(rbps) & (regulon["q"] <= q_max)].copy()
    calls = pd.read_parquet(calls_path)
    calls = calls[calls["rbp"].isin(rbps) & calls["switched"].fillna(False)].copy()
    keys = ["region", "module_id", "rbp"]
    candidate = calls.merge(
        regulon[
            keys
            + [
                "go_invisible",
                "q",
                "p",
                "enrichment",
                "module_size",
                "pool_source",
            ]
        ],
        on=keys,
        how="inner",
        suffixes=("", "_regulon"),
        validate="many_to_one",
    )
    candidate = candidate.drop_duplicates(["region", "module_id", "rbp", "gene"])

    region_tree = {region: tree for tree, region in _REGIONS}
    pair_frames: list[pd.DataFrame] = []
    pair_hashes: dict[str, str] = {}
    for region in sorted(candidate["region"].unique()):
        tree = region_tree.get(region)
        if tree is None:
            raise ValueError(f"No canonical analysis tree for region {region!r}")
        pair_path = rel(
            "real_data",
            tree,
            region,
            "_m",
            "isograph_vae",
            "module_interpret",
            "structure_switch_pairs.parquet",
        )
        pairs = pd.read_parquet(
            pair_path, columns=["gene_id", "transcript_id_1", "transcript_id_2"]
        )
        pairs["gene"] = _bare(pairs["gene_id"])
        pairs["region"] = region
        pair_frames.append(pairs)
        pair_hashes[region] = sha256(pair_path)

    pair_table = pd.concat(pair_frames, ignore_index=True)
    candidate = candidate.merge(
        pair_table,
        on=["region", "gene"],
        how="inner",
        validate="many_to_many",
    )
    counts = pd.read_parquet(counts_path, columns=["transcript_id", "rbp", "count"])
    counts = counts[counts["rbp"].isin(rbps)].copy()
    counts = (
        counts.groupby(["transcript_id", "rbp"], as_index=False)["count"]
        .sum()
        .rename(columns={"count": "motif_count"})
    )
    count_1 = counts.rename(
        columns={
            "transcript_id": "transcript_id_1",
            "motif_count": "motif_count_1",
        }
    )
    count_2 = counts.rename(
        columns={
            "transcript_id": "transcript_id_2",
            "motif_count": "motif_count_2",
        }
    )
    candidate = candidate.merge(
        count_1, on=["transcript_id_1", "rbp"], how="left", validate="many_to_one"
    )
    candidate = candidate.merge(
        count_2, on=["transcript_id_2", "rbp"], how="left", validate="many_to_one"
    )
    candidate[["motif_count_1", "motif_count_2"]] = (
        candidate[["motif_count_1", "motif_count_2"]].fillna(0).astype(int)
    )
    candidate = candidate[
        (candidate["motif_count_1"] > 0) != (candidate["motif_count_2"] > 0)
    ].copy()
    hashes = _source_hashes(cfg)
    candidate["regulon_sha256"] = hashes["regulon_sha256"]
    candidate["switch_calls_sha256"] = hashes["switch_calls_sha256"]
    candidate["motif_counts_sha256"] = hashes["motif_counts_sha256"]
    candidate["switch_pairs_sha256"] = candidate["region"].map(pair_hashes)
    candidate["frozen_date"] = str(cfg["frozen_date"])
    candidate["motif_scope"] = "intronic"
    candidate["q_max"] = q_max
    candidate = candidate.sort_values(
        ["rbp", "region", "module_id", "gene", "transcript_id_1", "transcript_id_2"]
    )
    observed_nominations = len(
        candidate[["region", "module_id", "rbp"]].drop_duplicates()
    )
    expected_nominations = int(source["expected_module_rbp_nominations"])
    if observed_nominations != expected_nominations:
        raise RuntimeError(
            "Candidate freeze changed: "
            f"{observed_nominations} module-RBP nominations observed, "
            f"{expected_nominations} expected"
        )
    _guard_candidate_identity(source, _candidate_identity_sha256(candidate))

    out = rel(cfg["candidate_manifest"])
    ensure_dir(out.parent)
    candidate.to_parquet(out, index=False, compression="zstd")
    return out


def _configured_rows(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    root = rel(cfg["output_root"])
    rows: list[dict[str, Any]] = []
    for dataset in cfg["datasets"]:
        common = {
            key: value
            for key, value in dataset.items()
            if key not in {"files", "samples"}
        }
        files = dataset.get("files", [])
        if not files:
            rows.append(
                {
                    **common,
                    "record_type": "dataset_access",
                    "url": "",
                    "local_path": "",
                    "role": "",
                    "expected_bytes": pd.NA,
                    "extract": False,
                    "status": (
                        "controlled_access"
                        if dataset["access"] == "controlled"
                        else "no_files_configured"
                    ),
                }
            )
            continue
        for file_cfg in files:
            path = root / file_cfg["path"]
            file_metadata = {
                key: value
                for key, value in file_cfg.items()
                if key
                not in {
                    "url",
                    "path",
                    "role",
                    "expected_bytes",
                    "sha256",
                    "extract",
                }
            }
            rows.append(
                {
                    **common,
                    **file_metadata,
                    "record_type": "dataset_file",
                    "url": file_cfg["url"],
                    "local_path": str(path),
                    "role": file_cfg["role"],
                    "expected_bytes": int(file_cfg["expected_bytes"]),
                    "expected_sha256": file_cfg.get("sha256", ""),
                    "extract": bool(file_cfg.get("extract", False)),
                    "status": "configured",
                }
            )
        archive_paths = [
            root / file_cfg["path"]
            for file_cfg in files
            if file_cfg.get("extract", False)
        ]
        source_archive = str(archive_paths[0]) if len(archive_paths) == 1 else ""
        for sample in dataset.get("samples", []):
            sample_common = {
                key: value for key, value in sample.items() if key not in {"files"}
            }
            for member in sample.get("files", []):
                member_path = root / member["path"]
                member_metadata = {
                    key: value
                    for key, value in member.items()
                    if key not in {"path", "role"}
                }
                rows.append(
                    {
                        **common,
                        **sample_common,
                        **member_metadata,
                        "record_type": "sample_file",
                        "url": "",
                        "local_path": str(member_path),
                        "source_archive": source_archive,
                        "role": member["role"],
                        "expected_bytes": pd.NA,
                        "expected_sha256": "",
                        "extract": False,
                        "status": "configured_member",
                    }
                )
    return rows


def _download(
    url: str,
    destination: Path,
    expected_bytes: int,
    expected_sha256: str = "",
) -> dict[str, Any]:
    ensure_dir(destination.parent)
    partial = destination.with_name(destination.name + ".part")

    if destination.exists():
        size = destination.stat().st_size
        if size != expected_bytes:
            raise RuntimeError(
                f"Existing file has {size} bytes, expected {expected_bytes}: {destination}"
            )
        digest = sha256(destination)
        if expected_sha256 and digest != expected_sha256:
            raise RuntimeError(
                f"SHA-256 mismatch for {destination}: {digest} != {expected_sha256}"
            )
        return {
            "status": "existing_verified",
            "bytes": size,
            "sha256": digest,
            "retrieved_at_utc": datetime.fromtimestamp(
                destination.stat().st_mtime, timezone.utc
            ).isoformat(),
            "last_modified": "",
            "retrieval_metadata_source": "file_mtime_recovered",
            "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        }

    offset = partial.stat().st_size if partial.exists() else 0
    if offset == expected_bytes:
        digest = sha256(partial)
        if expected_sha256 and digest != expected_sha256:
            raise RuntimeError(
                f"SHA-256 mismatch for complete partial {partial}: "
                f"{digest} != {expected_sha256}"
            )
        os.replace(partial, destination)
        return {
            "status": "resumed_verified",
            "bytes": expected_bytes,
            "sha256": digest,
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            "last_modified": "",
            "retrieval_metadata_source": "transfer_completed_at",
            "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        }
    if offset > expected_bytes:
        raise RuntimeError(
            f"Partial file exceeds expected size ({offset} > {expected_bytes}): {partial}"
        )
    headers = {"User-Agent": "IsoGraph-neuronal-CLIP/1.0"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=300) as response:
        status_code = response.getcode()
        append = offset > 0 and status_code == 206
        mode = "ab" if append else "wb"
        if offset and not append:
            offset = 0
        with partial.open(mode) as handle:
            shutil.copyfileobj(response, handle, length=CHUNK_SIZE)
        last_modified = response.headers.get("Last-Modified", "")

    size = partial.stat().st_size
    if size != expected_bytes:
        raise RuntimeError(
            f"Incomplete download has {size} bytes, expected {expected_bytes}: {partial}"
        )
    digest = sha256(partial)
    if expected_sha256 and digest != expected_sha256:
        raise RuntimeError(
            f"SHA-256 mismatch for {partial}: {digest} != {expected_sha256}"
        )
    os.replace(partial, destination)
    return {
        "status": "downloaded",
        "bytes": size,
        "sha256": digest,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "last_modified": last_modified,
        "retrieval_metadata_source": "http_response",
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def _tar_members_verified(archive: Path, destination: Path) -> tuple[bool, int]:
    with tarfile.open(archive) as handle:
        members = handle.getmembers()
    for member in members:
        target = destination / member.name
        if member.isfile() and (
            not target.is_file() or target.stat().st_size != member.size
        ):
            return False, len(members)
        if member.isdir() and not target.is_dir():
            return False, len(members)
    return True, len(members)


def _safe_extract_tar(archive: Path, destination: Path) -> dict[str, Any]:
    ensure_dir(destination)
    verified, member_count = _tar_members_verified(archive, destination)
    if verified:
        return {
            "extracted_files": member_count,
            "extract_path": str(destination),
            "extraction_status": "existing_verified",
        }
    root = destination.resolve()
    with tarfile.open(archive) as handle:
        members = handle.getmembers()
        for member in members:
            target = (destination / member.name).resolve()
            if root not in {target, *target.parents}:
                raise RuntimeError(f"Unsafe archive path {member.name!r} in {archive}")
        handle.extractall(destination, members=members, filter="data")
    verified, member_count = _tar_members_verified(archive, destination)
    if not verified:
        raise RuntimeError(f"Extracted members failed verification: {archive}")
    return {
        "extracted_files": member_count,
        "extract_path": str(destination),
        "extraction_status": "extracted_verified",
    }


def write_dataset_manifest(
    cfg: dict[str, Any], results: dict[str, dict[str, Any]] | None = None
) -> Path:
    results = results or {}
    rows = _configured_rows(cfg)
    for row in rows:
        outcome = results.get(row["local_path"])
        if outcome:
            row.update(outcome)
        else:
            row.setdefault("bytes", pd.NA)
            row.setdefault("sha256", "")
            row.setdefault("retrieved_at_utc", "")
            row.setdefault("last_modified", "")
            row.setdefault("retrieval_metadata_source", "")
            row.setdefault("verified_at_utc", "")
    manifest = pd.DataFrame(rows)
    out = rel(cfg["dataset_manifest"])
    ensure_dir(out.parent)
    if out.exists():
        previous = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
        if "record_type" not in previous.columns:
            previous["record_type"] = "dataset_file"
        previous_by_path = previous.drop_duplicates(
            ["record_type", "local_path"], keep="last"
        ).set_index(["record_type", "local_path"])
        for idx, row in manifest.iterrows():
            key = (str(row.get("record_type", "")), str(row.get("local_path", "")))
            if key not in previous_by_path.index:
                continue
            old = previous_by_path.loc[key]
            for field in RETRIEVAL_FIELDS:
                old_value = str(old.get(field, ""))
                current_source = str(row.get("retrieval_metadata_source", ""))
                if old_value and (
                    not str(row.get(field, ""))
                    or current_source == "file_mtime_recovered"
                ):
                    manifest.at[idx, field] = old_value

    archive_metadata = (
        manifest[manifest["record_type"] == "dataset_file"]
        .drop_duplicates("local_path", keep="last")
        .set_index("local_path")
    )
    for idx, row in manifest[manifest["record_type"] == "sample_file"].iterrows():
        source_archive = str(row.get("source_archive", ""))
        if not source_archive or source_archive not in archive_metadata.index:
            continue
        archive = archive_metadata.loc[source_archive]
        manifest.at[idx, "source_archive_sha256"] = archive.get("sha256", "")
        for field in ("retrieved_at_utc", "last_modified"):
            manifest.at[idx, field] = archive.get(field, "")
        manifest.at[idx, "retrieval_metadata_source"] = "source_archive"

    manifest.to_csv(out, sep="\t", index=False)
    return out


def _enriched_window_stats(path: Path, q_max: float) -> dict[str, int]:
    data_rows = 0
    significant_rows = 0
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if "qvalue" not in header:
            raise ValueError(f"Missing qvalue column in {path}")
        q_index = header.index("qvalue")
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            data_rows += 1
            if float(fields[q_index]) <= q_max:
                significant_rows += 1
    return {
        "data_rows": data_rows,
        "significant_rows": significant_rows,
    }


def _quantas_event_stats(path: Path) -> dict[str, int]:
    with gzip.open(path, "rt", newline=None) as handle:
        header = [
            value.rstrip("\r")
            for value in handle.readline().rstrip("\n").split("\t")
        ]
        required = {"AS.name", "FDR", "delta.exon.inclusion.rate(dI)"}
        if not required.issubset(header):
            missing = sorted(required - set(header))
            raise ValueError(f"Missing Quantas columns in {path}: {missing}")
        event_index = header.index("AS.name")
        event_ids: set[str] = set()
        rows = 0
        for line in handle:
            if not line.strip():
                continue
            rows += 1
            event_ids.add(line.rstrip("\n\r").split("\t")[event_index])
    return {"event_rows": rows, "unique_event_ids": len(event_ids)}


def _bed12_event_stats(path: Path) -> dict[str, int]:
    rows = 0
    event_ids: set[str] = set()
    with gzip.open(path, "rt", newline=None) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 12:
                raise ValueError(f"Expected BED12 row in {path}")
            rows += 1
            event_ids.add(fields[3])
    return {"coordinate_rows": rows, "coordinate_event_ids": len(event_ids)}


def _sample_results(cfg: dict[str, Any]) -> dict[str, dict[str, Any]]:
    root = rel(cfg["output_root"])
    now = datetime.now(timezone.utc).isoformat()
    results: dict[str, dict[str, Any]] = {}
    for dataset in cfg["datasets"]:
        q_max = float(dataset.get("enriched_window_q_max", 0.05))
        for sample in dataset.get("samples", []):
            for member in sample.get("files", []):
                path = root / member["path"]
                if not path.exists():
                    continue
                outcome: dict[str, Any] = {
                    "status": "member_verified",
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "verified_at_utc": now,
                    "retrieved_at_utc": "",
                    "last_modified": "",
                    "retrieval_metadata_source": "source_archive",
                }
                if member.get("format") == "enriched_windows_tsv_gz":
                    outcome.update(_enriched_window_stats(path, q_max))
                elif member.get("format") == "quantas_as_tsv_gz":
                    outcome.update(_quantas_event_stats(path))
                elif member.get("format") == "quantas_bed12_gz":
                    outcome.update(_bed12_event_stats(path))
                results[str(path)] = outcome
    return results


def write_context_qc(cfg: dict[str, Any], manifest_path: Path) -> Path:
    manifest = pd.read_csv(manifest_path, sep="\t", low_memory=False)
    samples = manifest[manifest["record_type"] == "sample_file"].copy()
    rows: list[dict[str, Any]] = []
    for context in cfg.get("contexts", []):
        context_id = context["context_id"]
        context_samples = samples[samples["context_id"] == context_id].copy()
        qc_type = context["qc_type"]
        row: dict[str, Any] = {
            **context,
            "configured_samples": int(context_samples["sample_id"].nunique()),
            "verified_sample_files": int(
                (context_samples["status"] == "member_verified").sum()
            ),
            "missing_sample_files": int(
                (context_samples["status"] != "member_verified").sum()
            ),
        }
        if qc_type == "paired_strand_bigwig":
            consensus_files = manifest[
                (manifest["record_type"] == "dataset_file")
                & (manifest["context_id"] == context_id)
            ]
            observed_consensus = (
                str(consensus_files.iloc[0]["sha256"])
                if len(consensus_files) == 1
                else ""
            )
            required_strands = set(context["required_strands"])
            grouped = context_samples.groupby("sample_id", dropna=False)
            complete_samples = grouped.filter(
                lambda frame: set(frame["strand"].dropna()) == required_strands
                and (frame["status"] == "member_verified").all()
            )["sample_id"].nunique()
            ip_replicates = context_samples.loc[
                context_samples["assay_role"] == "IP", "sample_id"
            ].nunique()
            input_replicates = context_samples.loc[
                context_samples["assay_role"] == "input", "sample_id"
            ].nunique()
            passed = (
                complete_samples == context_samples["sample_id"].nunique()
                and ip_replicates >= int(context["min_ip_replicates"])
                and input_replicates >= int(context["min_input_replicates"])
            )
            row.update(
                {
                    "complete_strand_samples": int(complete_samples),
                    "ip_replicates": int(ip_replicates),
                    "input_replicates": int(input_replicates),
                    "consensus_sha256_observed": observed_consensus,
                    "consensus_hash_matches": (
                        observed_consensus == context["consensus_sha256"]
                    ),
                    "qc_pass": passed,
                    "qc_status": (
                        "ready_for_context_specific_reconstruction"
                        if passed
                        else "fail_missing_replicate_signal"
                    ),
                    "primary_eligible": False,
                }
            )
        elif qc_type == "enriched_windows":
            per_sample = context_samples.drop_duplicates("sample_id")
            q_min = int(context["min_significant_windows_per_replicate"])
            passing = per_sample[
                (per_sample["status"] == "member_verified")
                & (per_sample["significant_rows"].fillna(0) >= q_min)
            ]
            passed = len(passing) >= int(context["min_replicates"])
            scope = context["analysis_scope"]
            row.update(
                {
                    "data_rows_total": int(per_sample["data_rows"].fillna(0).sum()),
                    "significant_rows_total": int(
                        per_sample["significant_rows"].fillna(0).sum()
                    ),
                    "replicates_passing": int(len(passing)),
                    "qc_pass": passed,
                    "qc_status": (
                        f"pass_{scope}" if passed else "not_estimable_low_yield"
                    ),
                    "primary_eligible": bool(passed and scope == "primary"),
                }
            )
        elif qc_type == "ctag_pooled_bedgraph":
            per_sample = context_samples.drop_duplicates("sample_id")
            verified = per_sample[per_sample["status"] == "member_verified"]
            pooled_replicates = pd.to_numeric(
                verified.get("biological_replicates_pooled", pd.Series(dtype=float)),
                errors="coerce",
            )
            passed = len(verified) == 1 and bool(
                pooled_replicates.ge(
                    int(context["min_biological_replicates_pooled"])
                ).all()
            )
            scope = context["analysis_scope"]
            row.update(
                {
                    "pooled_files_passing": int(len(verified)),
                    "biological_replicates_pooled": int(
                        pooled_replicates.max() if len(pooled_replicates) else 0
                    ),
                    "replicate_resolution": "pooled_only",
                    "qc_pass": passed,
                    "qc_status": (
                        f"pass_{scope}_pooled_only"
                        if passed
                        else "not_estimable_missing_pooled_bedgraph"
                    ),
                    "primary_eligible": False,
                }
            )
        elif qc_type == "perturbation_tables":
            statistics = context_samples[
                context_samples["role"] == "significant_event_statistics"
            ]
            coordinates = context_samples[
                context_samples["role"] == "significant_event_coordinates"
            ]
            event_rows = int(statistics["event_rows"].fillna(0).sum())
            event_ids = int(statistics["unique_event_ids"].fillna(0).sum())
            coordinate_rows = int(coordinates["coordinate_rows"].fillna(0).sum())
            coordinate_event_ids = int(
                coordinates["coordinate_event_ids"].fillna(0).sum()
            )
            passed = (
                len(statistics) == 1
                and len(coordinates) == 1
                and bool((statistics["status"] == "member_verified").all())
                and bool((coordinates["status"] == "member_verified").all())
                and event_rows >= int(context["min_events"])
                and coordinate_event_ids == event_ids
            )
            scope = context["analysis_scope"]
            row.update(
                {
                    "event_rows": event_rows,
                    "unique_event_ids": event_ids,
                    "coordinate_rows": coordinate_rows,
                    "coordinate_event_ids": coordinate_event_ids,
                    "event_coordinate_ids_match": event_ids == coordinate_event_ids,
                    "qc_pass": passed,
                    "qc_status": (
                        f"pass_{scope}_significant_event_catalog"
                        if passed
                        else "fail_perturbation_table_qc"
                    ),
                    "primary_eligible": bool(passed and scope == "orthogonal_primary"),
                }
            )
        elif qc_type == "controlled_access":
            scope = context["analysis_scope"]
            row.update(
                {
                    "qc_pass": False,
                    "qc_status": f"{scope}_controlled_access",
                    "primary_eligible": False,
                }
            )
        else:
            raise ValueError(f"Unknown context QC type {qc_type!r}")
        rows.append(row)

    qc = pd.DataFrame(rows)
    if "consensus_sha256_observed" in qc.columns:
        qc["consensus_duplicate_contexts"] = pd.NA
        consensus = qc["consensus_sha256_observed"].fillna("").astype(str)
        has_consensus = consensus != ""
        qc.loc[has_consensus, "consensus_duplicate_contexts"] = (
            qc.loc[has_consensus]
            .groupby("consensus_sha256_observed")["context_id"]
            .transform("nunique")
        )
    out = rel(cfg["context_qc_manifest"])
    ensure_dir(out.parent)
    qc.to_csv(out, sep="\t", index=False)
    return out


def existing_results(cfg: dict[str, Any]) -> dict[str, dict[str, Any]]:
    root = rel(cfg["output_root"])
    results: dict[str, dict[str, Any]] = {}
    for dataset in cfg["datasets"]:
        for file_cfg in dataset.get("files", []):
            destination = root / file_cfg["path"]
            if not destination.exists():
                continue
            size = destination.stat().st_size
            expected_bytes = int(file_cfg["expected_bytes"])
            if size != expected_bytes:
                raise RuntimeError(
                    f"Existing file has {size} bytes, expected {expected_bytes}: "
                    f"{destination}"
                )
            digest = sha256(destination)
            expected_sha256 = file_cfg.get("sha256", "")
            if expected_sha256 and digest != expected_sha256:
                raise RuntimeError(
                    f"SHA-256 mismatch for {destination}: "
                    f"{digest} != {expected_sha256}"
                )
            results[str(destination)] = {
                "status": "existing_verified",
                "bytes": size,
                "sha256": digest,
                "retrieved_at_utc": datetime.fromtimestamp(
                    destination.stat().st_mtime, timezone.utc
                ).isoformat(),
                "last_modified": "",
                "retrieval_metadata_source": "file_mtime_recovered",
                "verified_at_utc": datetime.now(timezone.utc).isoformat(),
            }
            if file_cfg.get("extract", False):
                extract_dir = destination.parent / f"{destination.stem}_files"
                verified, member_count = _tar_members_verified(destination, extract_dir)
                results[str(destination)].update(
                    {
                        "extracted_files": member_count if verified else pd.NA,
                        "extract_path": str(extract_dir) if verified else "",
                        "extraction_status": (
                            "existing_verified" if verified else "not_extracted"
                        ),
                    }
                )
    results.update(_sample_results(cfg))
    return results


def fetch(cfg: dict[str, Any], selected: set[str], extract: bool) -> Path:
    freeze_candidates(cfg)
    results = existing_results(cfg)
    write_dataset_manifest(cfg, results)
    root = rel(cfg["output_root"])

    for dataset in cfg["datasets"]:
        dataset_id = dataset["dataset_id"]
        if selected and dataset_id not in selected:
            continue
        if dataset["access"] != "public":
            continue
        for file_cfg in dataset.get("files", []):
            destination = root / file_cfg["path"]
            print(f"[{dataset_id}] {file_cfg['url']} -> {destination}", flush=True)
            outcome = _download(
                file_cfg["url"],
                destination,
                int(file_cfg["expected_bytes"]),
                file_cfg.get("sha256", ""),
            )
            if extract and file_cfg.get("extract", False):
                extract_dir = destination.parent / f"{destination.stem}_files"
                outcome.update(_safe_extract_tar(destination, extract_dir))
            results[str(destination)] = outcome
            print(
                f"[{dataset_id}] {outcome['status']}: "
                f"{outcome['bytes']} bytes sha256={outcome['sha256']}",
                flush=True,
            )
            write_dataset_manifest(cfg, results)
    results.update(_sample_results(cfg))
    manifest = write_dataset_manifest(cfg, results)
    write_context_qc(cfg, manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Freeze manifests and reproducibly acquire neuronal CLIP inputs."
    )
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser(
        "manifest", help="freeze candidate and dataset manifests only"
    )
    fetch_parser = subparsers.add_parser(
        "fetch", help="download public configured files with atomic resume"
    )
    fetch_parser.add_argument(
        "--dataset",
        action="append",
        default=[],
        help="dataset_id to fetch; repeat as needed (default: all public datasets)",
    )
    fetch_parser.add_argument(
        "--extract",
        action="store_true",
        help="safely extract archives marked extract=true after checksum calculation",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_yaml(args.config)
    if args.command == "manifest":
        candidate = freeze_candidates(cfg)
        dataset = write_dataset_manifest(cfg, existing_results(cfg))
        context_qc = write_context_qc(cfg, dataset)
        print(candidate)
        print(dataset)
        print(context_qc)
        return
    manifest = fetch(cfg, set(args.dataset), args.extract)
    print(manifest)


if __name__ == "__main__":
    main()
