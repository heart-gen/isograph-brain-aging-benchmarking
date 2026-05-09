from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import pandas as pd

from isograph.evaluation.metrics import module_recovery_score
from isograph.io.artifacts import load_dataset_bundle
from isograph.models.baseline import BaselineNetworkModel
from isograph.models.graph import GraphNetworkModel
from isograph.models.gpu_latent import GpuLatentNetworkModel
from isograph.models.latent import LatentNetworkModel
from isograph.models.vae import VaeNetworkModel
from isograph.models.wgcna import WgcnaNetworkModel
from isograph.workflow.config import (
    BaselineModelConfig,
    GpuLatentModelConfig,
    GraphModelConfig,
    LatentModelConfig,
    VaeModelConfig,
    WgcnaModelConfig,
)

from isograph_benchmark.benchmark.synthetic_data import ensure_dataset
from isograph_benchmark.benchmark.telemetry import (
    hardware_info,
    measured_run,
    reset_torch_peak_memory,
    slurm_info,
    software_versions,
    torch_info,
    torch_peak_memory,
    write_json,
)
from isograph_benchmark.paths import ensure_dir, rel


ARTIFACT_TABLES = {
    "module_table": "modules.parquet",
    "edge_table": "edges.parquet",
    "trait_table": "traits.parquet",
    "feature_scores": "feature_scores.parquet",
    "eigengene_table": "eigengenes.parquet",
    "module_gene_roles": "module_gene_roles.parquet",
}


def _row_to_dict(row: pd.Series) -> dict[str, Any]:
    return {key: None if pd.isna(value) else value for key, value in row.to_dict().items()}


def _wgcna_threads(row: pd.Series) -> int:
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus:
        return int(slurm_cpus)
    requested = row.get("requested_cpus")
    if pd.notna(requested):
        return int(requested)
    return 1


def _set_thread_env(threads: int) -> None:
    for key in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        os.environ[key] = str(threads)


def build_model(row: pd.Series):
    method = str(row["method"])
    seed = int(row["seed"])
    if method == "isograph_baseline":
        return BaselineNetworkModel(BaselineModelConfig(alpha=0.10, min_module_size=2))
    if method == "isograph_latent":
        return LatentNetworkModel(LatentModelConfig(alpha=0.10, min_module_size=2, n_components_cv_folds=3))
    if method == "isograph_graph":
        return GraphNetworkModel(GraphModelConfig(alpha=0.10, min_module_size=2, n_components_cv_folds=3))
    if method == "isograph_gpu_latent":
        return GpuLatentNetworkModel(GpuLatentModelConfig(alpha=0.10, min_module_size=2, random_state=seed))
    if method == "isograph_cpu_latent":
        return GpuLatentNetworkModel(GpuLatentModelConfig(alpha=0.10, min_module_size=2, random_state=seed, device="cpu"))
    if method == "isograph_vae":
        return VaeNetworkModel(
            VaeModelConfig(
                latent_dim_grid=[2, 4, 6, 8, 12],
                hidden_dim=128 if int(row["n_genes"]) <= 1000 else 256,
                n_epochs=300,
                patience=35,
                alpha=0.70,
                min_module_size=2,
                random_state=seed,
                device="cpu",
            )
        )
    if method == "isograph_vae_gpu":
        return VaeNetworkModel(
            VaeModelConfig(
                latent_dim_grid=[2, 4, 6, 8, 12],
                hidden_dim=128 if int(row["n_genes"]) <= 1000 else 256,
                n_epochs=300,
                patience=35,
                alpha=0.70,
                min_module_size=2,
                random_state=seed,
            )
        )
    if method == "isograph_vae_multiplex":
        return VaeNetworkModel(
            VaeModelConfig(
                latent_dim_grid=[2, 4, 6, 8, 12],
                hidden_dim=128 if int(row["n_genes"]) <= 1000 else 256,
                n_epochs=300,
                patience=35,
                alpha=0.70,
                alpha_switch=0.70,
                allow_abundance_abundance=True,
                alpha_abundance_grid=[0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90],
                min_module_size=2,
                random_state=seed,
                device="cpu",
            )
        )
    if method == "wgcna_gene":
        threads = _wgcna_threads(row)
        _set_thread_env(threads)
        return WgcnaNetworkModel(WgcnaModelConfig(min_module_size=2, random_state=seed))
    raise ValueError(f"Unknown method: {method}")


def write_artifacts(out_dir: Path, artifacts) -> None:
    for attr, filename in ARTIFACT_TABLES.items():
        table = getattr(artifacts, attr, None)
        if table is not None and not table.empty:
            table.to_parquet(out_dir / filename, index=False, compression="zstd")
    if artifacts.calibration is not None:
        write_json(out_dir / "calibration.json", artifacts.calibration)


def compute_metrics(artifacts, bundle) -> dict[str, Any]:
    truth_modules = bundle.truth_tables.get("truth_modules.parquet", pd.DataFrame())
    truth_switch = bundle.truth_tables.get("truth_switch.parquet", pd.DataFrame())
    truth_abundance = bundle.truth_tables.get("truth_abundance.parquet", pd.DataFrame())
    predicted_genes = set(artifacts.module_table["gene_id"]) if not artifacts.module_table.empty else set()
    switching_genes = set(truth_switch.loc[truth_switch["has_switch"], "gene_id"]) if not truth_switch.empty else set()
    nonswitching_genes = set(truth_switch.loc[~truth_switch["has_switch"], "gene_id"]) if not truth_switch.empty else set()
    abundance_genes = set(truth_abundance.loc[truth_abundance["has_abundance"], "gene_id"]) if not truth_abundance.empty else set()

    metrics: dict[str, Any] = {
        "module_recovery": module_recovery_score(artifacts.module_table, truth_modules),
        "n_predicted_modules": int(artifacts.module_table["module_id"].nunique()) if not artifacts.module_table.empty else 0,
        "n_edges": int(len(artifacts.edge_table)),
        "switch_gene_detection_rate": (
            len(predicted_genes & switching_genes) / len(switching_genes) if switching_genes else None
        ),
        "nonswitch_gene_module_rate": (
            len(predicted_genes & nonswitching_genes) / len(nonswitching_genes) if nonswitching_genes else None
        ),
        "abundance_gene_detection_rate": (
            len(predicted_genes & abundance_genes) / len(abundance_genes) if abundance_genes else None
        ),
    }

    roles = artifacts.module_gene_roles
    if roles is not None and not roles.empty and "module_role" in roles.columns:
        calibration = artifacts.calibration or {}
        metrics["selected_alpha_abundance"] = calibration.get("selected_alpha_abundance")
        for role in ("switch_only", "abundance_only", "coupled", "discordant"):
            role_genes = set(roles.loc[roles["module_role"] == role, "gene_id"])
            metrics[f"role_{role}_n"] = int(len(role_genes))
        if switching_genes:
            switch_role = set(roles.loc[roles["module_role"].isin(("switch_only", "coupled")), "gene_id"])
            metrics["role_switch_recall"] = len(switch_role & switching_genes) / len(switching_genes)
        if abundance_genes:
            abund_role = set(roles.loc[roles["module_role"].isin(("abundance_only", "coupled")), "gene_id"])
            metrics["role_abundance_recall"] = len(abund_role & abundance_genes) / len(abundance_genes)

    return metrics


def run(row: pd.Series, grid_path: Path, output_root: Path, dataset_root: Path, force: bool = False) -> Path:
    out_dir = ensure_dir(output_root / str(row["run_id"]))
    done = out_dir / "done.json"
    if done.exists() and not force:
        return done

    telemetry: dict[str, Any] = {
        "run": _row_to_dict(row),
        "hardware": hardware_info(),
        "slurm": slurm_info(),
        "software": software_versions(include_r=str(row["method"]) == "wgcna_gene"),
        "torch": torch_info(),
        "grid_path": str(grid_path),
    }
    if str(row["method"]) == "wgcna_gene":
        telemetry["wgcna"] = {
            "wgcna_threads": _wgcna_threads(row),
            "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
            "requested_cpus": int(row["requested_cpus"]),
        }

    status = "completed"
    metrics: dict[str, Any] = {}
    try:
        with measured_run() as measurement:
            dataset_path = ensure_dataset(row, dataset_root)
            bundle = load_dataset_bundle(dataset_path)
            reset_torch_peak_memory()
            model = build_model(row)
            artifacts = model.fit(
                transcript_counts=bundle.matrices["transcript_counts"],
                transcript_table=bundle.feature_tables["transcript"],
                sample_table=bundle.sample_table,
            )
            write_artifacts(out_dir, artifacts)
            metrics = compute_metrics(artifacts, bundle)
            telemetry["dataset_path"] = str(dataset_path)
        telemetry["measurement"] = measurement
        telemetry["torch"].update(torch_peak_memory())
    except Exception as exc:
        status = "failed"
        telemetry["error"] = {"type": type(exc).__name__, "message": str(exc)}

    telemetry["status"] = status
    telemetry["metrics"] = metrics
    write_json(out_dir / "telemetry.json", telemetry)
    if status == "failed":
        raise RuntimeError(f"Run {row['run_id']} failed: {telemetry['error']['message']}")
    write_json(done, {"status": status, "run_id": str(row["run_id"]), "metrics": metrics})
    return done


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", default="benchmark/00_design/_m/synthetic_run_grid.parquet")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-root", default="benchmark/01_synthetic/_o/runs")
    parser.add_argument("--dataset-root", default="benchmark/01_synthetic/_m/datasets")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    grid_path = Path(args.grid)
    if not grid_path.is_absolute():
        grid_path = rel(args.grid)
    grid = pd.read_parquet(grid_path)
    matches = grid.loc[grid["run_id"].astype(str) == str(args.run_id)]
    if len(matches) != 1:
        raise ValueError(f"Expected one row for run_id={args.run_id}, found {len(matches)}")
    run(
        matches.iloc[0],
        grid_path=grid_path,
        output_root=rel(args.output_root),
        dataset_root=rel(args.dataset_root),
        force=args.force,
    )


if __name__ == "__main__":
    main()
