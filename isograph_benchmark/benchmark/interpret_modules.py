from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from isograph.io.artifacts import load_dataset_bundle
from isograph_benchmark.interpretation import explain_artifact_modules
from isograph_benchmark.paths import ensure_dir, rel


DEFAULT_METHODS = ["isograph_vae", "isograph_vae_gpu", "isograph_vae_multiplex"]
METRIC_COLUMNS = [
    "gene_driver_truth_module_auroc",
    "gene_driver_switch_auroc",
    "top10_truth_module_precision",
    "top10_switch_precision",
    "switch_strength_auroc",
    "opposite_polarity_rate",
    "switch_driver_recall",
    "abundance_driver_recall",
    "explanation_fidelity",
]


def _field(row: pd.Series, name: str) -> Any:
    run_name = f"run_{name}"
    if run_name in row.index:
        return row[run_name]
    return row[name]


def _auc(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    mask = np.isfinite(scores)
    labels = labels[mask]
    scores = scores[mask]
    n_pos = int(labels.sum())
    n_neg = int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = rankdata(scores, method="average")
    pos_rank_sum = ranks[labels].sum()
    return float((pos_rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def _precision_at(gene_table: pd.DataFrame, truth_genes: set[str], k: int = 10) -> float:
    if gene_table.empty or not truth_genes:
        return np.nan
    ranked = gene_table.copy()
    ranked["_score"] = pd.to_numeric(ranked["r"], errors="coerce").abs()
    top = ranked.sort_values("_score", ascending=False).head(k)["gene_id"].astype(str)
    if top.empty:
        return np.nan
    return float(top.isin(truth_genes).mean())


def _best_truth_module(predicted_genes: set[str], truth_modules: pd.DataFrame) -> tuple[str | None, set[str]]:
    if not predicted_genes or truth_modules.empty:
        return None, set()
    best_module = None
    best_genes: set[str] = set()
    best_overlap = -1
    for module_id, frame in truth_modules.groupby("module_id"):
        truth_genes = set(frame["gene_id"].astype(str))
        overlap = len(predicted_genes & truth_genes)
        if overlap > best_overlap:
            best_overlap = overlap
            best_module = str(module_id)
            best_genes = truth_genes
    return best_module, best_genes


def _role_attribution_metrics(
    module_gene_roles_path: Path,
    switching_genes: set[str],
    abundance_genes: set[str],
    predicted_genes: set[str],
) -> dict[str, float]:
    """Compute role attribution fidelity using module_gene_roles.parquet if available."""
    metrics: dict[str, float] = {
        "switch_driver_recall": np.nan,
        "abundance_driver_recall": np.nan,
        "explanation_fidelity": np.nan,
    }
    if not module_gene_roles_path.exists():
        return metrics
    roles = pd.read_parquet(module_gene_roles_path)
    if roles.empty or "module_role" not in roles.columns or "gene_id" not in roles.columns:
        return metrics
    roles = roles[roles["gene_id"].isin(predicted_genes)].copy()

    switch_role_genes = set(roles.loc[roles["module_role"].isin(("switch_only", "coupled")), "gene_id"].astype(str))
    abund_role_genes = set(roles.loc[roles["module_role"].isin(("abundance_only", "coupled")), "gene_id"].astype(str))

    switch_in_pred = switching_genes & predicted_genes
    abund_in_pred = abundance_genes & predicted_genes

    if switch_in_pred:
        metrics["switch_driver_recall"] = len(switch_role_genes & switch_in_pred) / len(switch_in_pred)
    if abund_in_pred:
        metrics["abundance_driver_recall"] = len(abund_role_genes & abund_in_pred) / len(abund_in_pred)

    correct = len(switch_role_genes & switch_in_pred) + len(abund_role_genes & abund_in_pred)
    total = len(switch_in_pred) + len(abund_in_pred)
    if total > 0:
        metrics["explanation_fidelity"] = correct / total
    return metrics


def evaluate_interpret_output(
    artifact_dir: Path,
    explain_dir: Path,
    dataset_path: Path,
    run_info: dict[str, Any],
) -> pd.DataFrame:
    bundle = load_dataset_bundle(dataset_path)
    truth_modules = bundle.truth_tables.get("truth_modules.parquet", pd.DataFrame())
    truth_switch = bundle.truth_tables.get("truth_switch.parquet", pd.DataFrame())
    truth_abundance = bundle.truth_tables.get("truth_abundance.parquet", pd.DataFrame())
    switching_genes = set(truth_switch.loc[truth_switch["has_switch"], "gene_id"].astype(str))
    abundance_genes = set(truth_abundance.loc[truth_abundance["has_abundance"], "gene_id"].astype(str)) if not truth_abundance.empty else set()

    modules = pd.read_parquet(artifact_dir / "modules.parquet")
    if modules.empty or "module_id" not in modules.columns:
        return pd.DataFrame()
    module_sizes = modules.groupby("module_id")["gene_id"].nunique().to_dict()

    # Load module_gene_roles once for the whole run (used for role attribution metrics)
    roles_path = artifact_dir / "module_gene_roles.parquet"
    all_predicted_genes = set(modules["gene_id"].astype(str))
    run_role_metrics = _role_attribution_metrics(roles_path, switching_genes, abundance_genes, all_predicted_genes)

    rows: list[dict[str, Any]] = []
    for module_id in sorted(modules["module_id"].astype(str).unique()):
        module_dir = explain_dir / module_id
        gene_path = module_dir / "gene_driver_table.parquet"
        tx_path = module_dir / "transcript_polarity_table.parquet"
        if not gene_path.exists() or not tx_path.exists():
            continue

        gene_table = pd.read_parquet(gene_path)
        tx_table = pd.read_parquet(tx_path)
        predicted_genes = set(modules.loc[modules["module_id"].astype(str) == module_id, "gene_id"].astype(str))
        truth_module_id, truth_module_genes = _best_truth_module(predicted_genes, truth_modules)

        gene_scores = pd.to_numeric(gene_table.get("r", pd.Series(dtype=float)), errors="coerce").abs().to_numpy()
        gene_ids = gene_table.get("gene_id", pd.Series(dtype=str)).astype(str).to_numpy()
        truth_module_labels = np.array([gene in truth_module_genes for gene in gene_ids], dtype=bool)
        switch_labels = np.array([gene in switching_genes for gene in gene_ids], dtype=bool)

        switch_strength_auc = np.nan
        opposite_polarity_rate = np.nan
        if not tx_table.empty and {"gene_id", "r", "switch_strength"}.issubset(tx_table.columns):
            tx_tmp = tx_table.copy()
            tx_tmp["abs_switch_strength"] = pd.to_numeric(tx_tmp["switch_strength"], errors="coerce").abs()
            per_gene_strength = tx_tmp.groupby("gene_id")["abs_switch_strength"].max()
            strength_genes = per_gene_strength.index.astype(str).to_numpy()
            strength_labels = np.array([gene in switching_genes for gene in strength_genes], dtype=bool)
            switch_strength_auc = _auc(strength_labels, per_gene_strength.to_numpy(dtype=float))

            polarity_rows = []
            for gene_id, frame in tx_tmp.groupby("gene_id"):
                if str(gene_id) not in switching_genes or len(frame) < 2:
                    continue
                r = pd.to_numeric(frame["r"], errors="coerce").to_numpy(dtype=float)
                r = r[np.isfinite(r)]
                if len(r) >= 2:
                    polarity_rows.append(bool(r.min() < 0 and r.max() > 0))
            if polarity_rows:
                opposite_polarity_rate = float(np.mean(polarity_rows))

        rows.append(
            {
                **run_info,
                "module_id": module_id,
                "truth_module_id": truth_module_id,
                "n_predicted_module_genes": int(module_sizes.get(module_id, 0)),
                "n_truth_module_genes": int(len(truth_module_genes)),
                "n_truth_module_overlap": int(len(predicted_genes & truth_module_genes)),
                "gene_driver_truth_module_auroc": _auc(truth_module_labels, gene_scores),
                "gene_driver_switch_auroc": _auc(switch_labels, gene_scores),
                "top10_truth_module_precision": _precision_at(gene_table, truth_module_genes, k=10),
                "top10_switch_precision": _precision_at(gene_table, switching_genes, k=10),
                "switch_strength_auroc": switch_strength_auc,
                "opposite_polarity_rate": opposite_polarity_rate,
                **run_role_metrics,
            }
        )
    return pd.DataFrame(rows)


def _aggregate_run_metrics(module_metrics: pd.DataFrame, run_info: dict[str, Any], status: str) -> dict[str, Any]:
    row = {**run_info, "status": status, "n_modules_evaluated": int(len(module_metrics))}
    for metric in METRIC_COLUMNS:
        if metric in module_metrics.columns:
            values = pd.to_numeric(module_metrics[metric], errors="coerce").to_numpy(dtype=float)
            row[metric] = float(np.nanmean(values)) if np.isfinite(values).any() else np.nan
        else:
            row[metric] = np.nan
    return row


def run_interpretation(
    results_path: Path,
    run_root: Path,
    dataset_root: Path,
    output_root: Path,
    out_dir: Path,
    methods: list[str],
    force: bool = False,
    limit: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    results = pd.read_parquet(results_path)
    if "status" in results.columns:
        results = results.loc[results["status"] == "completed"].copy()
    method_col = "run_method" if "run_method" in results.columns else "method"
    results = results.loc[results[method_col].isin(methods)].copy()
    if limit is not None:
        results = results.head(limit)

    ensure_dir(output_root)
    ensure_dir(out_dir)
    run_rows: list[dict[str, Any]] = []
    module_parts: list[pd.DataFrame] = []

    for _, row in results.iterrows():
        run_id = str(_field(row, "run_id"))
        dataset_id = str(_field(row, "dataset_id"))
        run_info = {
            "run_id": run_id,
            "dataset_id": dataset_id,
            "scenario": str(_field(row, "scenario")),
            "method": str(_field(row, "method")),
            "replicate": int(_field(row, "replicate")),
            "seed": int(_field(row, "seed")),
        }
        artifact_dir = run_root / run_id
        dataset_path = dataset_root / dataset_id
        explain_dir = output_root / run_id
        modules_path = artifact_dir / "modules.parquet"
        if not modules_path.exists() or not dataset_path.exists():
            run_rows.append(_aggregate_run_metrics(pd.DataFrame(), run_info, "missing_inputs"))
            continue
        modules = pd.read_parquet(modules_path)
        if modules.empty or "module_id" not in modules.columns:
            run_rows.append(_aggregate_run_metrics(pd.DataFrame(), run_info, "no_modules"))
            continue

        manifest = explain_dir / "module_explanation_manifest.json"
        if force or not manifest.exists():
            explain_artifact_modules(
                artifact_dir=artifact_dir,
                bundle_path=dataset_path,
                output_dir=explain_dir,
            )

        module_metrics = evaluate_interpret_output(artifact_dir, explain_dir, dataset_path, run_info)
        if not module_metrics.empty:
            module_parts.append(module_metrics)
        run_rows.append(_aggregate_run_metrics(module_metrics, run_info, "completed"))

    run_df = pd.DataFrame(run_rows)
    module_df = pd.concat(module_parts, ignore_index=True) if module_parts else pd.DataFrame()
    run_df.to_parquet(out_dir / "synthetic_interpret_results.parquet", index=False, compression="zstd")
    module_df.to_parquet(out_dir / "synthetic_interpret_module_metrics.parquet", index=False, compression="zstd")
    return run_df, module_df


def _bootstrap_ci(values: np.ndarray, n_iter: int, alpha: float, seed: int = 0) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    means = np.empty(n_iter)
    for i in range(n_iter):
        means[i] = rng.choice(values, size=len(values), replace=True).mean()
    return tuple(np.quantile(means, [alpha / 2, 1 - alpha / 2]))


def summarize_interpretation(
    run_df: pd.DataFrame,
    out_dir: Path,
    bootstrap_iterations: int = 10000,
    confidence_level: float = 0.95,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    alpha = 1 - confidence_level
    if run_df.empty or "status" not in run_df.columns:
        summary = pd.DataFrame(columns=["scenario", "method", "metric", "n", "mean", "median", "ci_low", "ci_high"])
        summary.to_parquet(out_dir / "synthetic_interpret_summary.parquet", index=False, compression="zstd")
        return summary
    completed = run_df.loc[run_df["status"] == "completed"].copy()
    for (scenario, method), group in completed.groupby(["scenario", "method"], dropna=False):
        for metric in METRIC_COLUMNS:
            values = pd.to_numeric(group[metric], errors="coerce").to_numpy(dtype=float)
            lo, hi = _bootstrap_ci(values, n_iter=bootstrap_iterations, alpha=alpha)
            rows.append(
                {
                    "scenario": scenario,
                    "method": method,
                    "metric": metric,
                    "n": int(np.isfinite(values).sum()),
                    "mean": float(np.nanmean(values)) if np.isfinite(values).any() else np.nan,
                    "median": float(np.nanmedian(values)) if np.isfinite(values).any() else np.nan,
                    "ci_low": lo,
                    "ci_high": hi,
                }
            )
    summary = pd.DataFrame(rows)
    summary.to_parquet(out_dir / "synthetic_interpret_summary.parquet", index=False, compression="zstd")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate IsoGraph module interpretation on synthetic truth.")
    parser.add_argument("--results", default="benchmark/01_synthetic/_m/synthetic_results.parquet")
    parser.add_argument("--run-root", default="benchmark/01_synthetic/_o/runs")
    parser.add_argument("--dataset-root", default="benchmark/01_synthetic/_m/datasets")
    parser.add_argument("--output-root", default="benchmark/03_interpret/_o/runs")
    parser.add_argument("--out-dir", default="benchmark/03_interpret/_m")
    parser.add_argument("--method", action="append", dest="methods", help="Method to evaluate. Can be repeated.")
    parser.add_argument("--force", action="store_true", help="Regenerate existing explain-module outputs.")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of completed runs for smoke tests.")
    parser.add_argument("--bootstrap-iterations", type=int, default=10000)
    parser.add_argument("--confidence-level", type=float, default=0.95)
    args = parser.parse_args()

    methods = args.methods or DEFAULT_METHODS
    run_df, _ = run_interpretation(
        results_path=rel(args.results),
        run_root=rel(args.run_root),
        dataset_root=rel(args.dataset_root),
        output_root=rel(args.output_root),
        out_dir=ensure_dir(rel(args.out_dir)),
        methods=methods,
        force=args.force,
        limit=args.limit,
    )
    summary = summarize_interpretation(
        run_df,
        out_dir=ensure_dir(rel(args.out_dir)),
        bootstrap_iterations=args.bootstrap_iterations,
        confidence_level=args.confidence_level,
    )
    print(f"Wrote {len(run_df):,} run rows and {len(summary):,} summary rows to {rel(args.out_dir)}")


if __name__ == "__main__":
    main()
