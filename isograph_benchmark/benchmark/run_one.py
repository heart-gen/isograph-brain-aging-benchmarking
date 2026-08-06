from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from isograph.evaluation.metrics import module_recovery_score
from isograph.io.artifacts import load_dataset_bundle
from isograph.models.baseline import BaselineNetworkModel
from isograph.models.graph import GraphNetworkModel
from isograph.models.latent import LatentNetworkModel
from isograph.models.vae import VaeNetworkModel
from isograph.models.wgcna import WgcnaNetworkModel
from isograph.workflow.config import (
    BaselineModelConfig,
    GraphModelConfig,
    LatentModelConfig,
    VaeModelConfig,
    WgcnaModelConfig,
)

from isograph_benchmark.benchmark.spearman_leiden import SpearmanLeidenConfig, SpearmanLeidenModel
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
        # Same VAE model/config as isograph_vae, run on GPU. Kept for the
        # supplementary "what VAE GPU does" figure; not a core comparator.
        return VaeNetworkModel(
            VaeModelConfig(
                latent_dim_grid=[2, 4, 6, 8, 12],
                hidden_dim=128 if int(row["n_genes"]) <= 1000 else 256,
                n_epochs=300,
                patience=35,
                alpha=0.70,
                min_module_size=2,
                random_state=seed,
                device="cuda",
            )
        )
    if method == "isograph_vae_residual":
        # Gap #6 ablation: identical to isograph_vae but residualizes the recorded
        # nuisance covariates (RNA degradation, cell composition, batch, library
        # depth). The isograph_vae vs isograph_vae_residual contrast on confounded
        # scenarios is the WITH/WITHOUT residualization headline result.
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
                residualize_covariates=["RIN", "neuron_frac", "batch", "library_size"],
            )
        )
    if method == "isograph_vae_reliability":
        # Approach #3: degradation-aware switch reliability driven by an OBSERVED
        # sample-level 3' coverage covariate (median TIN) -- a sharper proxy of the
        # degradation artifact than raw RIN, while staying one vector per sample.
        # Per-gene switch reliability downweights switch-switch edges from
        # degradation-aligned genes so they fall back to the abundance channel.
        # CRITICAL: this requires the abundance channel to be enabled (multiplex);
        # without an abundance fallback, downweighting switch edges only fragments
        # modules and can HURT under heavy degradation (verified). The contrast vs
        # isograph_vae_multiplex therefore isolates the reliability contribution on
        # top of an enabled abundance channel.
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
                switch_reliability_weighting=True,
                degradation_covariate="median_tin",
                # Match isograph_vae_multiplex so the contrast isolates reliability
                # weighting (not clustering). Leiden also prevents the giant-module
                # collapse the dense abundance channel triggers under connected components.
                leiden_resolution=2.0,
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
                # Dense abundance edges can fuse distinct modules into one giant
                # component under connected-components clustering (verified giant-module
                # collapse, ~62% on some seeds). Resolution-controlled Leiden splits
                # them (the same fix WGCNA's dynamic tree cut applies), lifting coupled-
                # scenario recovery 0.53->0.68 with no effect when no giant forms.
                leiden_resolution=2.0,
            )
        )
    if method == "isograph_spearman_leiden":
        return SpearmanLeidenModel(SpearmanLeidenConfig(min_r=0.30, leiden_resolution=1.0, min_module_size=2))
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


_FEATURE_META_COLS = ("feature_id", "gene_id", "feature_type", "n_transcripts")


def _uniform_module_eigengenes(module_table, feature_scores, sample_ids):
    """Module eigengenes computed identically for every method.

    eigengene(module) = mean, over that module's gene feature-scores, of the per-sample
    score (the same formula IsoGraph's own base model uses). Crucially we compute it
    from each method's RAW ``feature_scores`` + its module partition rather than reading
    the model's emitted ``eigengene_table`` -- some baselines (e.g. Spearman+Leiden) do
    not emit one, and persisted feature_scores are raw for every method (VAE denoising
    lives in the embedding it clusters, not the scores). So this isolates the one thing
    the A1 comparison is about: did the method recover the co-switching module?
    """
    if module_table is None or module_table.empty or feature_scores is None or feature_scores.empty:
        return [], np.empty((0, 0))
    ids = [s for s in sample_ids if s in feature_scores.columns]
    rows: list[np.ndarray] = []
    mods: list[str] = []
    for module_id, genes in module_table.groupby("module_id")["gene_id"]:
        sub = feature_scores.loc[feature_scores["gene_id"].isin(set(genes)), ids]
        if sub.empty:
            continue
        mods.append(module_id)
        rows.append(sub.to_numpy(dtype=float).mean(axis=0))
    return mods, (np.vstack(rows) if rows else np.empty((0, len(ids))))


def genetic_recovery_metrics(artifacts, bundle) -> dict[str, Any]:
    """A1: does the method recover the planted cis-genetic anchoring of modules?

    Ground truth (``truth_genetics.parquet`` + ``genotypes`` matrix) plants, per anchored
    module, one bi-allelic variant whose dosage shifts that module's switch latent, plus
    an equal number of unlinked null variants. The per-allele effect is modest, so the
    genotype->PSI signal is only detectable by pooling the module's genes into an
    eigengene -- i.e. only after the co-switching module is recovered. For each planted
    variant we take the method's best (most-associated) recovered-module eigengene and
    test it against the genotype (Bonferroni over the predicted modules). Reported:

    * ``genetic_anchor_recall``      -- fraction of ANCHORED variants recovered;
    * ``genetic_anchor_fpr``         -- fraction of NULL variants spuriously recovered;
    * ``genetic_anchor_best_r2_mean``-- mean best eigengene-genotype R^2 over anchored
                                        variants (graded recovery, robust to threshold).
    Returns ``{}`` when the dataset carries no planted genetics (all other scenarios).
    """
    truth_gen = bundle.truth_tables.get("truth_genetics.parquet", pd.DataFrame())
    if truth_gen.empty or "genotypes" not in bundle.matrices:
        return {}
    genotypes = np.asarray(bundle.matrices["genotypes"], dtype=float)  # (n_snps, n_samples)
    n_anchored = int(truth_gen["anchored"].sum())
    n_null = int((~truth_gen["anchored"]).sum())
    base = {
        "n_planted_mqtl": n_anchored,
        "n_null_variants": n_null,
        "genetic_anchor_recall": 0.0,
        "genetic_anchor_fpr": 0.0,
        "genetic_anchor_best_r2_mean": 0.0,
    }

    sample_ids = list(bundle.sample_table["sample_id"])
    mods, E = _uniform_module_eigengenes(artifacts.module_table, artifacts.feature_scores, sample_ids)
    if E.shape[0] == 0:
        return base  # no modules recovered -> zero anchoring by construction
    # Restrict/reorder genotype columns to the eigengene sample order.
    pos = {s: i for i, s in enumerate(sample_ids)}
    used_ids = [s for s in sample_ids if s in set(artifacts.feature_scores.columns)]
    g_cols = [pos[s] for s in used_ids]
    G = genotypes[:, g_cols]
    n_pred = E.shape[0]

    recalls: list[float] = []
    fprs: list[float] = []
    best_r2s: list[float] = []
    for i, anchored in enumerate(truth_gen["anchored"].tolist()):
        g = G[i]
        if np.std(g) == 0:
            continue
        best_r2 = 0.0
        best_p = 1.0
        for j in range(n_pred):
            e = E[j]
            if not np.all(np.isfinite(e)) or np.std(e) == 0:
                continue
            r, p = stats.pearsonr(e, g)
            if r * r > best_r2:
                best_r2 = float(r * r)
                best_p = float(p)
        detected = 1.0 if min(1.0, best_p * max(n_pred, 1)) < 0.05 else 0.0
        if anchored:
            recalls.append(detected)
            best_r2s.append(best_r2)
        else:
            fprs.append(detected)

    base["genetic_anchor_recall"] = float(np.mean(recalls)) if recalls else 0.0
    base["genetic_anchor_fpr"] = float(np.mean(fprs)) if fprs else 0.0
    base["genetic_anchor_best_r2_mean"] = float(np.mean(best_r2s)) if best_r2s else 0.0
    return base


def compute_metrics(artifacts, bundle) -> dict[str, Any]:
    truth_modules = bundle.truth_tables.get("truth_modules.parquet", pd.DataFrame())
    truth_switch = bundle.truth_tables.get("truth_switch.parquet", pd.DataFrame())
    truth_abundance = bundle.truth_tables.get("truth_abundance.parquet", pd.DataFrame())
    predicted_genes = set(artifacts.module_table["gene_id"]) if not artifacts.module_table.empty else set()
    switching_genes = set(truth_switch.loc[truth_switch["has_switch"], "gene_id"]) if not truth_switch.empty else set()
    nonswitching_genes = set(truth_switch.loc[~truth_switch["has_switch"], "gene_id"]) if not truth_switch.empty else set()
    abundance_genes = set(truth_abundance.loc[truth_abundance["has_abundance"], "gene_id"]) if not truth_abundance.empty else set()

    # module_recovery_score: mean best-Jaccard over truth modules.
    # For each ground-truth module T_i, find the predicted module P_j that
    # maximises |T_i ∩ P_j| / |T_i ∪ P_j|, then average over all T_i.
    # Score ∈ [0, 1]; 1 = perfect recovery of all truth modules.
    # See benchmark/README.md (Methods → Metrics) for the full formula and interpretation.
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

    # A1: genetic-anchoring recovery (no-op dict for non-genetic scenarios).
    metrics.update(genetic_recovery_metrics(artifacts, bundle))

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
