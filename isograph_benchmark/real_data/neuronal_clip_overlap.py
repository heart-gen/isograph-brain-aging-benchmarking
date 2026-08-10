"""Call neuronal CLIP support in callable, matched switch windows."""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from isograph_benchmark.config import load_yaml
from isograph_benchmark.paths import ensure_dir
from isograph_benchmark.real_data.neuronal_clip_fetch import sha256
from isograph_benchmark.real_data.neuronal_clip_validation import (
    DEFAULT_CONFIG,
    _configured_path,
)


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


def _atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(f".{path.name}.tmp")
    frame.to_csv(temporary, sep="\t", index=False)
    os.replace(temporary, path)


def _markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return ""
    headers = [str(column) for column in frame.columns]
    rows = [headers, ["---"] * len(headers)]
    for values in frame.itertuples(index=False, name=None):
        rows.append([str(value) for value in values])
    return "\n".join("| " + " | ".join(row) + " |" for row in rows)


def _as_bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes"}
    return bool(value)


def _context_map(cfg: dict[str, Any]) -> dict[str, dict[str, Any]]:
    contexts = {str(row["context_id"]): dict(row) for row in cfg["contexts"]}
    if len(contexts) != len(cfg["contexts"]):
        raise ValueError("context_id values must be unique")
    return contexts


def _sample_files(
    cfg: dict[str, Any], context_id: str, formats: set[str]
) -> list[dict[str, Any]]:
    root = _configured_path(cfg["output_root"])
    rows: list[dict[str, Any]] = []
    for dataset in cfg["datasets"]:
        for sample in dataset.get("samples", []):
            if str(sample.get("context_id", "")) != context_id:
                continue
            for file_cfg in sample.get("files", []):
                if str(file_cfg.get("format", "")) not in formats:
                    continue
                path = root / str(file_cfg["path"])
                rows.append(
                    {
                        "dataset_id": str(dataset["dataset_id"]),
                        "sample_id": str(sample["sample_id"]),
                        "context_id": context_id,
                        "assay_role": str(sample["assay_role"]),
                        "replicate": int(sample["replicate"]),
                        "strand": str(file_cfg.get("strand", "")),
                        "format": str(file_cfg["format"]),
                        "path": str(path),
                        "sha256": sha256(path) if path.is_file() else "",
                    }
                )
    return rows


def _load_windows(stage: dict[str, Any]) -> pd.DataFrame:
    window_dir = _configured_path(stage["window_dir"])
    windows = pd.read_parquet(window_dir / "window_manifest.parquet")
    eligibility = pd.read_parquet(_configured_path(stage["eligibility_manifest"]))
    eligible_ids = set(
        eligibility.loc[eligibility["confirmatory_eligible"], "candidate_id"].astype(str)
    )
    windows = windows[
        windows["window_width"].eq(int(stage["primary_width"]))
        & windows["match_status"].eq("matched")
        & windows["candidate_id"].astype(str).isin(eligible_ids)
    ].copy()
    if windows.empty:
        raise ValueError("No confirmatory-eligible matched primary-width windows")
    duplicated = windows.duplicated(["match_id", "window_role"], keep=False)
    if duplicated.any():
        raise ValueError(
            f"{int(duplicated.sum())} rows violate one case and one control per match_id"
        )
    return windows


def _interval_calls(
    windows: pd.DataFrame,
    tested: pd.DataFrame,
    q_max: float,
    require_same_strand: bool,
    callable_gene_ids: set[str] | None = None,
) -> pd.DataFrame:
    tested = tested.copy()
    tested["start"] = pd.to_numeric(tested["start"], errors="raise").astype(int)
    tested["end"] = pd.to_numeric(tested["end"], errors="raise").astype(int)
    tested["qvalue"] = pd.to_numeric(tested["qvalue"], errors="coerce")
    tested["input"] = pd.to_numeric(tested.get("input", np.nan), errors="coerce")
    tested["clip"] = pd.to_numeric(tested.get("clip", np.nan), errors="coerce")
    chrom_column = "chr" if "chr" in tested.columns else "chrom"
    tested["_chrom"] = tested[chrom_column].astype(str)
    tested["_strand"] = (
        tested["strand"].astype(str) if "strand" in tested else "."
    )
    rows: list[dict[str, Any]] = []
    for window in windows.itertuples(index=False):
        subset = tested[tested["_chrom"].eq(str(window.chrom))]
        if require_same_strand:
            subset = subset[subset["_strand"].eq(str(window.strand))]
        overlap = subset[
            subset["start"].lt(int(window.end))
            & subset["end"].gt(int(window.start))
        ]
        significant = overlap[overlap["qvalue"].le(q_max)]
        if callable_gene_ids is None:
            callable_window = bool(len(overlap))
        else:
            callable_window = str(window.gene_id).split(".", 1)[0] in callable_gene_ids
        rows.append(
            {
                "window_id": window.window_id,
                "n_tested_overlap": int(len(overlap)),
                "n_significant_overlap": int(len(significant)),
                "callable": callable_window,
                "bound": callable_window and bool(len(significant)),
                "minimum_qvalue": (
                    float(overlap["qvalue"].min()) if len(overlap) else math.nan
                ),
                "maximum_input": (
                    float(overlap["input"].max()) if len(overlap) else math.nan
                ),
                "maximum_clip": (
                    float(overlap["clip"].max()) if len(overlap) else math.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def _tdp43_calls(
    cfg: dict[str, Any],
    context: dict[str, Any],
    windows: pd.DataFrame,
) -> pd.DataFrame:
    stage = cfg["overlap_stage"]["enriched_windows"]
    files = _sample_files(cfg, str(context["context_id"]), {"enriched_windows_tsv_gz"})
    if len(files) < int(context["min_replicates"]):
        raise FileNotFoundError(
            f"{context['context_id']} has {len(files)} configured enriched-window files"
        )
    selected = windows[windows["rbp"].eq(str(context["rbp"]))].copy()
    rows: list[pd.DataFrame] = []
    for file_row in files:
        path = Path(file_row["path"])
        if not path.is_file():
            raise FileNotFoundError(path)
        tested = pd.read_csv(path, sep="\t", compression="gzip")
        if "gene_id" not in tested:
            raise ValueError(f"{path} lacks gene_id for processed-gene callability")
        callable_genes = {
            str(gene_id).split(".", 1)[0]
            for gene_id in tested["gene_id"].dropna().astype(str)
        }
        call = _interval_calls(
            selected,
            tested,
            float(stage["q_max"]),
            bool(stage["require_same_strand"]),
            callable_genes,
        )
        call["context_id"] = str(context["context_id"])
        call["sample_id"] = file_row["sample_id"]
        call["replicate"] = int(file_row["replicate"])
        call["qc_type"] = "enriched_windows"
        call["callability_basis"] = str(stage["callability_basis"])
        call["mappability_status"] = str(stage["mappability_status"])
        call["source_path"] = str(path)
        call["source_sha256"] = file_row["sha256"]
        rows.append(call)
    return pd.concat(rows, ignore_index=True)


def _run_bigwig_helper(
    cfg: dict[str, Any],
    context: dict[str, Any],
    windows: pd.DataFrame,
    output_dir: Path,
) -> pd.DataFrame:
    stage = cfg["overlap_stage"]["paired_strand_bigwig"]
    files = _sample_files(cfg, str(context["context_id"]), {"bigwig"})
    expected = int(context["min_ip_replicates"]) + int(context["min_input_replicates"])
    complete_samples = len({row["sample_id"] for row in files})
    if complete_samples < expected:
        raise FileNotFoundError(
            f"{context['context_id']} has {complete_samples}/{expected} complete samples"
        )
    missing = [row["path"] for row in files if not Path(row["path"]).is_file()]
    if missing:
        raise FileNotFoundError(missing[0])
    selected = windows[windows["rbp"].eq(str(context["rbp"]))].copy()
    rscript = Path(str(stage["rscript"]))
    helper = _configured_path(stage["helper"])
    if not rscript.is_file():
        raise FileNotFoundError(rscript)
    if not helper.is_file():
        raise FileNotFoundError(helper)
    with tempfile.TemporaryDirectory(prefix="ptbp2_signal_", dir=output_dir) as temporary:
        temporary_dir = Path(temporary)
        window_path = temporary_dir / "windows.tsv"
        file_path = temporary_dir / "files.tsv"
        result_path = temporary_dir / "signals.tsv"
        selected[["window_id", "chrom", "start", "end", "strand"]].to_csv(
            window_path, sep="\t", index=False
        )
        pd.DataFrame(files)[
            ["sample_id", "assay_role", "replicate", "strand", "path"]
        ].to_csv(file_path, sep="\t", index=False)
        subprocess.run(
            [str(rscript), str(helper), str(window_path), str(file_path), str(result_path)],
            check=True,
        )
        signal = pd.read_csv(result_path, sep="\t")
    scale = float(stage["normalization_scale"])
    signal["signal_cpm"] = signal["mean_signal"] / (
        signal["library_total"] / scale
    )
    key = ["window_id", "replicate"]
    duplicate = signal.duplicated(key + ["assay_role"], keep=False)
    if duplicate.any():
        raise ValueError("PTBP2 signal rows are not unique by window/replicate/assay role")
    input_signal = signal[signal["assay_role"].eq("input")].rename(
        columns={
            "sample_id": "input_sample_id",
            "signal_cpm": "input_cpm",
            "bigwig_path": "input_path",
        }
    )
    ip_signal = signal[signal["assay_role"].eq("IP")].rename(
        columns={
            "sample_id": "ip_sample_id",
            "signal_cpm": "ip_cpm",
            "bigwig_path": "ip_path",
        }
    )
    paired = ip_signal[key + ["ip_sample_id", "ip_cpm", "ip_path"]].merge(
        input_signal[key + ["input_sample_id", "input_cpm", "input_path"]],
        on=key,
        how="inner",
        validate="one_to_one",
    )
    pseudocount = float(stage["pseudocount_cpm"])
    paired["log2_ip_input"] = np.log2(
        (paired["ip_cpm"] + pseudocount)
        / (paired["input_cpm"] + pseudocount)
    )
    paired["callable"] = paired["input_cpm"].ge(float(stage["min_input_cpm"]))
    if bool(stage["require_positive_input"]):
        paired["callable"] &= paired["input_cpm"].gt(0)
    paired["bound"] = paired["callable"] & paired["log2_ip_input"].ge(
        float(stage["min_log2_ip_input"])
    )
    paired["context_id"] = str(context["context_id"])
    paired["sample_id"] = paired["ip_sample_id"] + "|" + paired["input_sample_id"]
    paired["qc_type"] = "paired_strand_bigwig"
    paired["callability_basis"] = str(stage["callability_basis"])
    paired["mappability_status"] = str(stage["mappability_status"])
    paired["source_path"] = paired["ip_path"] + "|" + paired["input_path"]
    paired["source_sha256"] = "see_dataset_manifest"
    return paired


def _consensus(
    calls: pd.DataFrame, context: dict[str, Any], cfg: dict[str, Any]
) -> pd.DataFrame:
    stage_key = (
        "enriched_windows"
        if context["qc_type"] == "enriched_windows"
        else "paired_strand_bigwig"
    )
    stage = cfg["overlap_stage"][stage_key]
    summary = (
        calls.groupby(["context_id", "window_id"], as_index=False)
        .agg(
            configured_replicates=("replicate", "nunique"),
            callable_replicates=("callable", "sum"),
            bound_replicates=("bound", "sum"),
            callability_basis=("callability_basis", "first"),
            mappability_status=("mappability_status", "first"),
        )
    )
    summary["context_callable"] = summary["callable_replicates"].ge(
        int(stage["min_callable_replicates"])
    )
    summary["context_bound"] = summary["context_callable"] & summary[
        "bound_replicates"
    ].ge(int(stage["min_bound_replicates"]))
    summary["consensus_rule"] = (
        f"callable>={stage['min_callable_replicates']};bound>={stage['min_bound_replicates']}"
    )
    return summary


def _matched_effects(
    consensus: pd.DataFrame,
    windows: pd.DataFrame,
    context: dict[str, Any],
    min_discordant: int,
) -> dict[str, Any]:
    joined = windows[
        ["window_id", "candidate_id", "rbp", "match_id", "window_role"]
    ].merge(consensus, on="window_id", how="inner", validate="one_to_one")
    pairs = joined.pivot(
        index=["candidate_id", "rbp", "match_id"],
        columns="window_role",
        values=["context_callable", "context_bound"],
    )
    required = [
        ("context_callable", "case"),
        ("context_callable", "control"),
        ("context_bound", "case"),
        ("context_bound", "control"),
    ]
    if any(column not in pairs.columns for column in required):
        callable_pairs = pairs.iloc[0:0]
    else:
        callable_pairs = pairs[
            pairs[("context_callable", "case")]
            & pairs[("context_callable", "control")]
        ]
    case_bound = callable_pairs.get(("context_bound", "case"), pd.Series(dtype=bool))
    control_bound = callable_pairs.get(
        ("context_bound", "control"), pd.Series(dtype=bool)
    )
    case_only = int((case_bound & ~control_bound).sum())
    control_only = int((~case_bound & control_bound).sum())
    discordant = case_only + control_only
    exact = (
        float(binomtest(case_only, discordant, 0.5, alternative="greater").pvalue)
        if discordant
        else math.nan
    )
    odds_ratio = (
        case_only / control_only
        if control_only
        else (case_only + 0.5) / (control_only + 0.5)
    )
    analysis_scope = str(context["analysis_scope"])
    if discordant < min_discordant:
        status = "descriptive_underpowered"
        confirmatory_pvalue = math.nan
    elif analysis_scope == "primary":
        status = "confirmatory"
        confirmatory_pvalue = exact
    else:
        status = "sensitivity"
        confirmatory_pvalue = math.nan
    return {
        "context_id": str(context["context_id"]),
        "rbp": str(context["rbp"]),
        "analysis_scope": analysis_scope,
        "evidence_tier": int(context["evidence_tier"]),
        "qc_type": str(context["qc_type"]),
        "callable_matched_sets": int(len(callable_pairs)),
        "case_bound_sets": int(case_bound.sum()),
        "control_bound_sets": int(control_bound.sum()),
        "case_only_discordant": case_only,
        "control_only_discordant": control_only,
        "informative_discordant_sets": discordant,
        "matched_odds_ratio": float(odds_ratio),
        "exact_one_sided_pvalue": exact,
        "confirmatory_pvalue": confirmatory_pvalue,
        "effect_status": status,
    }


def run_overlap(
    cfg: dict[str, Any],
    config_path: Path,
    output_dir: Path,
    selected_contexts: set[str],
) -> dict[str, Path]:
    stage = cfg["overlap_stage"]
    ensure_dir(output_dir)
    windows = _load_windows(stage)
    context_map = _context_map(cfg)
    context_qc_path = _configured_path(stage["context_qc_manifest"])
    context_qc = pd.read_csv(context_qc_path, sep="\t").set_index("context_id")
    configured = selected_contexts or set(stage["human_contexts"])
    unknown = configured - set(context_map)
    if unknown:
        raise ValueError(f"Unknown context IDs: {sorted(unknown)}")
    call_frames: list[pd.DataFrame] = []
    consensus_frames: list[pd.DataFrame] = []
    effects: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []
    for context_id in sorted(configured):
        context = context_map[context_id]
        qc = context_qc.loc[context_id]
        qc_pass = _as_bool(qc["qc_pass"])
        if not qc_pass:
            status_rows.append(
                {
                    "context_id": context_id,
                    "analysis_status": str(qc["qc_status"]),
                    "included": False,
                    "reason": "context_qc_failed",
                }
            )
            continue
        rbp = str(context["rbp"])
        if ";" in rbp or not windows["rbp"].eq(rbp).any():
            status_rows.append(
                {
                    "context_id": context_id,
                    "analysis_status": "not_testable_no_eligible_candidates",
                    "included": False,
                    "reason": "no_confirmatory_eligible_windows_for_exact_rbp",
                }
            )
            continue
        print(f"processing {context_id}", flush=True)
        if context["qc_type"] == "enriched_windows":
            calls = _tdp43_calls(cfg, context, windows)
        elif context["qc_type"] == "paired_strand_bigwig":
            calls = _run_bigwig_helper(cfg, context, windows, output_dir)
        else:
            status_rows.append(
                {
                    "context_id": context_id,
                    "analysis_status": "not_implemented_for_qc_type",
                    "included": False,
                    "reason": str(context["qc_type"]),
                }
            )
            continue
        context_consensus = _consensus(calls, context, cfg)
        effect = _matched_effects(
            context_consensus,
            windows[windows["rbp"].eq(rbp)],
            context,
            int(stage["min_informative_discordant_sets"]),
        )
        call_frames.append(calls)
        consensus_frames.append(context_consensus)
        effects.append(effect)
        status_rows.append(
            {
                "context_id": context_id,
                "analysis_status": effect["effect_status"],
                "included": True,
                "reason": "",
            }
        )
    calls = pd.concat(call_frames, ignore_index=True) if call_frames else pd.DataFrame()
    consensus = (
        pd.concat(consensus_frames, ignore_index=True)
        if consensus_frames
        else pd.DataFrame()
    )
    effects_frame = pd.DataFrame(effects)
    statuses = pd.DataFrame(status_rows)
    outputs = {
        "calls": output_dir / "overlap_calls.parquet",
        "consensus": output_dir / "context_consensus.parquet",
        "effects": output_dir / "primary_effects.parquet",
        "contexts": output_dir / "context_analysis_status.tsv",
        "json": output_dir / "neuronal_clip_overlap.json",
        "report": output_dir / "NEURONAL_CLIP_OVERLAP.md",
    }
    _atomic_parquet(calls, outputs["calls"])
    _atomic_parquet(consensus, outputs["consensus"])
    _atomic_parquet(effects_frame, outputs["effects"])
    _atomic_tsv(statuses, outputs["contexts"])
    provenance = {
        "config": str(config_path),
        "config_sha256": sha256(config_path),
        "window_manifest_sha256": sha256(
            _configured_path(stage["window_dir"]) / "window_manifest.parquet"
        ),
        "eligibility_manifest_sha256": sha256(
            _configured_path(stage["eligibility_manifest"])
        ),
        "context_qc_sha256": sha256(context_qc_path),
        "selected_contexts": sorted(configured),
        "analyzed_contexts": sorted(effects_frame["context_id"].tolist())
        if not effects_frame.empty
        else [],
        "seed": int(cfg["seed"]),
    }
    _atomic_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", outputs["json"])
    lines = [
        "# Neuronal CLIP callable-window overlap",
        "",
        "Current NOVA1 nominations are quarantined upstream and cannot enter these results.",
        "",
        "## Context status",
        "",
        _markdown_table(statuses) if not statuses.empty else "No contexts configured.",
        "",
        "## Matched effects",
        "",
        _markdown_table(effects_frame)
        if not effects_frame.empty
        else "No context produced an estimable effect.",
        "",
        "TDP-43 callability requires the matched pair's gene in each replicate's study-processed CLIP universe; binding still requires strand-correct overlap with a q<=0.05 enriched window. PTBP2 uses replicate-paired, library-normalized IP/input signal and is reported as signal support rather than a reconstructed CLAM peak.",
        "",
    ]
    _atomic_text("\n".join(lines), outputs["report"])
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir")
    parser.add_argument("--context", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = _configured_path(args.config)
    cfg = load_yaml(config_path)
    configured_output = _configured_path(cfg["overlap_stage"]["output_dir"])
    if args.context and not args.output_dir:
        raise ValueError("Partial context runs require --output-dir")
    output_dir = _configured_path(args.output_dir) if args.output_dir else configured_output
    outputs = run_overlap(cfg, config_path, output_dir, set(args.context))
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
