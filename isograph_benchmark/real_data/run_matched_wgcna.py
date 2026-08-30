"""Matched-feature WGCNA real-data baselines.

This is the real-data counterpart to the synthetic benchmark's same-input WGCNA
comparison. The existing real-data ``wgcna_gene`` scripts are classical gene-
abundance WGCNA baselines. This module instead builds the exact IsoGraph feature
matrix from each dataset bundle and runs WGCNA on selected feature rows:

* ``switch``: switch-coordinate rows only (multi-isoform genes).
* ``multiplex``: abundance + switch rows, matching the benchmark WGCNA input.

Outputs are side-by-side directories under each region's ``_m/`` folder:
``wgcna_switch_only`` and ``wgcna_multiplex``. The existing ``wgcna_gene`` outputs
are never overwritten.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from isograph.features.channels import feature_sample_columns, gene_feature_channels, make_feature_scores
from isograph.features.residualize import build_design_matrix, residualize_rows
from isograph.io.artifacts import load_dataset_bundle
from isograph.models.base import FitArtifacts, compute_module_gene_roles
from isograph.models.wgcna import _RUNNER_R
from isograph.workflow.config import WgcnaModelConfig

from isograph_benchmark.paths import ensure_dir, region_store, rel
from isograph_benchmark.real_data.run_models import (
    GTEX_REGIONS,
    _filter_expressed_transcripts,
    _gtex_qc_covariate_table,
    _rnaseqc_covariate_table,
    _save_age_artifacts,
    _save_diagnosis_artifacts,
)

BRAINSEQ_AGING_REGIONS = ["caudate", "hippocampus", "dlpfc"]
VARIANT_DIRS = {
    "switch": "wgcna_switch_only",
    "multiplex": "wgcna_multiplex",
}

BRAINSEQ_COVARIATES = [
    "Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
    "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5",
]
GTEX_COVARIATES = ["SEX", "SMRIN", "SMTSISCH", "SMMAPRT"]
SCZD_COVARIATES = ["Age", *BRAINSEQ_COVARIATES]


def _check_rscript() -> None:
    if shutil.which("Rscript") is None:
        raise RuntimeError("Rscript not found in PATH; activate the R/WGCNA-capable environment.")


def _feature_subset(matrix: np.ndarray, info: pd.DataFrame, variant: str) -> tuple[np.ndarray, pd.DataFrame]:
    if variant == "multiplex":
        keep = np.ones(len(info), dtype=bool)
    elif variant == "switch":
        keep = (info["feature_type"].astype(str) == "switch").to_numpy()
    else:
        raise ValueError(f"unknown WGCNA feature variant: {variant!r}")
    out_info = info.loc[keep].reset_index(drop=True)
    out_matrix = matrix[keep]
    if out_matrix.shape[0] < 4:
        raise ValueError(f"{variant}: only {out_matrix.shape[0]} feature rows after filtering")
    return out_matrix, out_info


def _run_wgcna(matrix: np.ndarray, feature_info: pd.DataFrame, cfg: WgcnaModelConfig) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    _check_rscript()
    feature_ids = feature_info["feature_id"].astype(str).tolist()
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        input_csv = tmp / "feature_matrix.csv"
        output_json = tmp / "wgcna_result.json"
        df = pd.DataFrame(matrix, index=feature_ids)
        df.index.name = "feature_id"
        df.to_csv(input_csv)
        cmd = [
            "Rscript", str(_RUNNER_R),
            f"input={input_csv}",
            f"output={output_json}",
            f"power={'auto' if cfg.power is None else cfg.power}",
            f"power_range={','.join(str(p) for p in cfg.power_range)}",
            f"sft_r2={cfg.sft_r2_threshold}",
            f"min_module_size={cfg.min_module_size}",
            f"merge_cut_height={cfg.merge_cut_height}",
            f"deep_split={cfg.deep_split}",
            f"network_type={cfg.network_type}",
            f"seed={cfg.random_state}",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=cfg.timeout_seconds)
        if proc.returncode != 0:
            raise RuntimeError(f"WGCNA R script failed (exit {proc.returncode}):\n{proc.stderr}")
        result = json.loads(output_json.read_text())

    feature_to_gene = feature_info.set_index("feature_id")["gene_id"].astype(str).to_dict()
    modules = pd.DataFrame(result.get("modules", []))
    if modules.empty:
        module_table = pd.DataFrame(columns=["gene_id", "module_id"])
    else:
        module_table = modules.rename(columns={"gene_id": "feature_id"})
        module_table["gene_id"] = module_table["feature_id"].map(feature_to_gene)
        module_table = (
            module_table.dropna(subset=["gene_id"])[["gene_id", "module_id"]]
            .drop_duplicates()
            .reset_index(drop=True)
        )

    edges = pd.DataFrame(result.get("edges", []))
    if edges.empty:
        edge_table = pd.DataFrame(columns=["source", "target", "weight"])
    else:
        edge_table = edges.rename(columns={"source": "source_feature_id", "target": "target_feature_id"})
        edge_table["source"] = edge_table["source_feature_id"].map(feature_to_gene)
        edge_table["target"] = edge_table["target_feature_id"].map(feature_to_gene)
        edge_table = edge_table.dropna(subset=["source", "target"])
        edge_table = edge_table.loc[edge_table["source"] != edge_table["target"]].reset_index(drop=True)
    return module_table, edge_table, result.get("calibration", {})


def _eigengene_table(module_table: pd.DataFrame, feature_scores: pd.DataFrame) -> pd.DataFrame:
    sample_cols = feature_sample_columns(feature_scores)
    if module_table.empty:
        return pd.DataFrame(columns=["module_id", *sample_cols])
    rows = {}
    for module_id in sorted(module_table["module_id"].unique()):
        genes = module_table.loc[module_table["module_id"] == module_id, "gene_id"]
        subset = feature_scores.loc[feature_scores["gene_id"].isin(genes)]
        rows[module_id] = subset[sample_cols].to_numpy(dtype=float).mean(axis=0)
    return pd.DataFrame(rows, index=sample_cols).T.reset_index().rename(columns={"index": "module_id"})


def _fit_artifacts(bundle, transcript_counts: np.ndarray, transcript_table: pd.DataFrame,
                   covariates: list[str], variant: str, cfg: WgcnaModelConfig) -> FitArtifacts:
    matrix, info = gene_feature_channels(transcript_counts, transcript_table)
    if matrix.size:
        design = build_design_matrix(bundle.sample_table, covariates)
        matrix = residualize_rows(matrix, design)
    matrix, info = _feature_subset(matrix, info, variant)
    module_table, edge_table, calibration = _run_wgcna(matrix, info, cfg)
    feature_scores = make_feature_scores(matrix, info, bundle.sample_table)
    eigengenes = _eigengene_table(module_table, feature_scores)
    roles = compute_module_gene_roles(module_table, feature_scores, bundle.sample_table)
    return FitArtifacts(
        module_table=module_table,
        edge_table=edge_table,
        trait_table=pd.DataFrame(columns=["module_id", "trait", "effect", "pvalue"]),
        feature_scores=feature_scores,
        eigengene_table=eigengenes,
        module_gene_roles=roles,
        calibration=calibration,
    )


def _cfg(seed: int, timeout_seconds: int) -> WgcnaModelConfig:
    return WgcnaModelConfig(
        min_module_size=30,
        random_state=seed,
        timeout_seconds=timeout_seconds,
        residualize_covariates=[],
    )


def run_brainseq_aging(regions: list[str], variants: list[str], seed: int, timeout_seconds: int) -> None:
    for region in regions:
        bundle = load_dataset_bundle(rel("inputs", "bundles", "brainseq_v1", region))
        tc, tt = _filter_expressed_transcripts(bundle.matrices["transcript_counts"], bundle.feature_tables["transcript"])
        for variant in variants:
            print(f"[brainseq-aging/{region}/{variant}] matched-feature WGCNA", flush=True)
            art = _fit_artifacts(bundle, tc, tt, BRAINSEQ_COVARIATES, variant, _cfg(seed, timeout_seconds))
            out = ensure_dir(region_store("brainseq", region, VARIANT_DIRS[variant]))
            _save_age_artifacts(
                art, out, bundle.sample_table, BRAINSEQ_COVARIATES, age_col="Age",
                label=f"brainseq/{region}/{variant}", qc_table=_rnaseqc_covariate_table(region),
            )


def run_gtex_aging(regions: list[str], variants: list[str], seed: int, timeout_seconds: int) -> None:
    for region in regions:
        bundle = load_dataset_bundle(rel("inputs", "bundles", "gtex_v11_brain", region))
        for variant in variants:
            print(f"[gtex-aging/{region}/{variant}] matched-feature WGCNA", flush=True)
            art = _fit_artifacts(
                bundle, bundle.matrices["transcript_counts"], bundle.feature_tables["transcript"],
                GTEX_COVARIATES, variant, _cfg(seed, timeout_seconds),
            )
            out = ensure_dir(region_store("gtex", region, VARIANT_DIRS[variant]))
            _save_age_artifacts(
                art, out, bundle.sample_table, GTEX_COVARIATES, age_col="AGE",
                label=f"gtex/{region}/{variant}", qc_table=_gtex_qc_covariate_table(bundle.sample_table),
            )


def run_brainseq_sczd(variants: list[str], seed: int, timeout_seconds: int) -> None:
    bundle = load_dataset_bundle(rel("inputs", "bundles", "brainseq_sczd", "caudate"))
    for variant in variants:
        print(f"[brainseq-sczd/caudate/{variant}] matched-feature WGCNA", flush=True)
        art = _fit_artifacts(
            bundle, bundle.matrices["transcript_counts"], bundle.feature_tables["transcript"],
            BRAINSEQ_COVARIATES, variant, _cfg(seed, timeout_seconds),
        )
        out = ensure_dir(region_store("brainseq", "caudate_sczd", VARIANT_DIRS[variant]))
        _save_diagnosis_artifacts(
            art, out, bundle, covariate_cols=SCZD_COVARIATES, label=f"brainseq_sczd/{variant}",
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis", choices=["brainseq-aging", "gtex-aging", "brainseq-sczd"])
    parser.add_argument("--region", action="append", dest="regions")
    parser.add_argument("--variant", choices=sorted(VARIANT_DIRS), action="append", dest="variants")
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--timeout-seconds", type=int, default=7200)
    args = parser.parse_args()

    variants = args.variants or ["switch", "multiplex"]
    if args.analysis == "brainseq-aging":
        run_brainseq_aging(args.regions or BRAINSEQ_AGING_REGIONS, variants, args.seed, args.timeout_seconds)
    elif args.analysis == "gtex-aging":
        run_gtex_aging(args.regions or GTEX_REGIONS, variants, args.seed, args.timeout_seconds)
    else:
        run_brainseq_sczd(variants, args.seed, args.timeout_seconds)


if __name__ == "__main__":
    main()
