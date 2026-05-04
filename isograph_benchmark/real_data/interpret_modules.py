from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.interpretation import explain_artifact_modules, summarize_explain_results
from isograph_benchmark.paths import ensure_dir, rel


DEFAULT_GTF_PATH = Path(
    "/ocean/projects/bio250020p/shared/resources/genomes/human/gencode-v47/gtf/"
    "gencode.v47.primary_assembly.annotation.gtf"
)
DEFAULT_GTF_CACHE = rel("real_data", "_m", "tmp", "gencode.v47.primary_assembly.annotation.gtf_cache.parquet")


@dataclass(frozen=True)
class InterpretTarget:
    collection: str
    dataset: str
    bundle_path: Path
    artifact_dir: Path
    association_paths: tuple[Path, ...]


def _existing_gtex_regions() -> list[str]:
    bundle_root = rel("inputs", "bundles", "gtex_v11_brain")
    if not bundle_root.exists():
        return []
    return sorted(path.name for path in bundle_root.iterdir() if path.is_dir())


def iter_targets(scope: str) -> list[InterpretTarget]:
    targets: list[InterpretTarget] = []
    if scope in {"all", "brainseq"}:
        for region in ["caudate", "hippocampus", "dlpfc"]:
            artifact_dir = rel("real_data", "brainseq", region, "_m", "isograph_vae")
            targets.append(
                InterpretTarget(
                    collection="brainseq",
                    dataset=region,
                    bundle_path=rel("inputs", "bundles", "brainseq_v1", region),
                    artifact_dir=artifact_dir,
                    association_paths=(
                        artifact_dir / "age_linear.parquet",
                        artifact_dir / "age_spline.parquet",
                    ),
                )
            )
        sczd_artifact_dir = rel("real_data", "brainseq", "caudate_sczd", "_m", "isograph_vae")
        targets.append(
            InterpretTarget(
                collection="brainseq",
                dataset="caudate_sczd",
                bundle_path=rel("inputs", "bundles", "brainseq_sczd", "caudate"),
                artifact_dir=sczd_artifact_dir,
                association_paths=(sczd_artifact_dir / "diagnosis_assoc.parquet",),
            )
        )

    if scope in {"all", "gtex"}:
        for region in _existing_gtex_regions():
            artifact_dir = rel("real_data", "gtex", region, "_m", "isograph_vae")
            targets.append(
                InterpretTarget(
                    collection="gtex",
                    dataset=region,
                    bundle_path=rel("inputs", "bundles", "gtex_v11_brain", region),
                    artifact_dir=artifact_dir,
                    association_paths=(
                        artifact_dir / "age_linear.parquet",
                        artifact_dir / "age_spline.parquet",
                    ),
                )
            )
    return targets


def _read_modules(artifact_dir: Path) -> pd.DataFrame:
    modules_path = artifact_dir / "modules.parquet"
    if not modules_path.exists():
        return pd.DataFrame()
    modules = pd.read_parquet(modules_path)
    if modules.empty or "module_id" not in modules.columns:
        return pd.DataFrame()
    return modules


def select_modules(
    target: InterpretTarget,
    fdr_threshold: float = 0.10,
    top_n_fallback: int = 5,
    all_modules: bool = False,
) -> tuple[list[str], pd.DataFrame]:
    modules = _read_modules(target.artifact_dir)
    if modules.empty:
        return [], pd.DataFrame()
    all_ids = sorted(modules["module_id"].astype(str).unique())
    if all_modules:
        return all_ids, pd.DataFrame({"module_id": all_ids, "selection_reason": "all_modules"})

    assoc_parts: list[pd.DataFrame] = []
    for path in target.association_paths:
        if not path.exists():
            continue
        assoc = pd.read_parquet(path)
        if assoc.empty or "module_id" not in assoc.columns:
            continue
        assoc = assoc.copy()
        assoc["source"] = path.name
        assoc_parts.append(assoc)
    if not assoc_parts:
        fallback = all_ids[:top_n_fallback]
        return fallback, pd.DataFrame({"module_id": fallback, "selection_reason": "fallback_no_association"})

    assoc_all = pd.concat(assoc_parts, ignore_index=True)
    assoc_all["module_id"] = assoc_all["module_id"].astype(str)
    score_col = "fdr" if "fdr" in assoc_all.columns else "pvalue"
    assoc_all[score_col] = pd.to_numeric(assoc_all[score_col], errors="coerce")

    selected = assoc_all.loc[assoc_all[score_col] <= fdr_threshold].copy()
    reason = f"fdr<={fdr_threshold}" if score_col == "fdr" else f"pvalue<={fdr_threshold}"
    if selected.empty:
        selected = assoc_all.sort_values(score_col, na_position="last").head(top_n_fallback).copy()
        reason = f"top{top_n_fallback}_fallback"

    selected = selected.drop_duplicates("module_id", keep="first")
    selected = selected.loc[selected["module_id"].isin(all_ids)].copy()
    if selected.empty:
        fallback = all_ids[:top_n_fallback]
        return fallback, pd.DataFrame({"module_id": fallback, "selection_reason": "fallback_module_intersection"})
    selected["selection_reason"] = reason
    return selected["module_id"].tolist(), selected


def _summarize_existing_output(target: InterpretTarget, output_dir: Path, module_ids: list[str]) -> pd.DataFrame:
    modules = _read_modules(target.artifact_dir)
    module_sizes = modules.groupby("module_id")["gene_id"].nunique().to_dict() if not modules.empty else {}
    rows = []
    for module_id in module_ids:
        module_dir = output_dir / module_id
        gene_path = module_dir / "gene_driver_table.parquet"
        tx_path = module_dir / "transcript_polarity_table.parquet"
        if not gene_path.exists() or not tx_path.exists():
            continue
        gene_drivers = pd.read_parquet(gene_path)
        transcript_polarity = pd.read_parquet(tx_path)
        top_gene = None
        top_gene_abs_r = np.nan
        if not gene_drivers.empty and "r" in gene_drivers.columns:
            idx = pd.to_numeric(gene_drivers["r"], errors="coerce").abs().fillna(0).to_numpy().argmax()
            top = gene_drivers.iloc[idx]
            top_gene = top.get("gene_id")
            top_gene_abs_r = abs(float(top.get("r", np.nan)))
        max_switch_strength = np.nan
        if not transcript_polarity.empty and "switch_strength" in transcript_polarity.columns:
            max_switch_strength = float(pd.to_numeric(transcript_polarity["switch_strength"], errors="coerce").max())
        rows.append(
            {
                "collection": target.collection,
                "dataset": target.dataset,
                "module_id": module_id,
                "n_module_genes": int(module_sizes.get(module_id, 0)),
                "n_gene_drivers": int(len(gene_drivers)),
                "n_transcript_features": int(len(transcript_polarity)),
                "top_gene_id": top_gene,
                "top_gene_abs_r": top_gene_abs_r,
                "max_switch_strength": max_switch_strength,
            }
        )
    return pd.DataFrame(rows)


def run_target(
    target: InterpretTarget,
    fdr_threshold: float,
    top_n_fallback: int,
    all_modules: bool,
    annotation_table: Path | None,
    gtf_path: Path | None,
    gtf_cache: Path | None,
    max_transcripts_per_gene: int,
    max_pairs_per_gene: int,
    force: bool,
    plot: bool,
) -> pd.DataFrame:
    if not target.bundle_path.exists() or not (target.artifact_dir / "modules.parquet").exists():
        return pd.DataFrame()
    module_ids, selection = select_modules(
        target,
        fdr_threshold=fdr_threshold,
        top_n_fallback=top_n_fallback,
        all_modules=all_modules,
    )
    output_dir = ensure_dir(target.artifact_dir / "module_interpret")
    if selection is not None and not selection.empty:
        selection.to_parquet(output_dir / "module_selection.parquet", index=False, compression="zstd")
    if not module_ids:
        return pd.DataFrame()

    manifest = output_dir / "module_explanation_manifest.json"
    manifest_modules: set[str] = set()
    if manifest.exists():
        manifest_data = json.loads(manifest.read_text())
        manifest_modules = set(map(str, manifest_data.get("module_ids", [])))
    needs_run = force or not manifest.exists() or not set(module_ids).issubset(manifest_modules)
    if needs_run:
        results = explain_artifact_modules(
            artifact_dir=target.artifact_dir,
            bundle_path=target.bundle_path,
            output_dir=output_dir,
            module_ids=module_ids,
            annotation_table=annotation_table,
            gtf_path=gtf_path if annotation_table is None else None,
            gtf_cache=gtf_cache,
            max_transcripts_per_gene=max_transcripts_per_gene,
            max_pairs_per_gene=max_pairs_per_gene,
            plot=plot,
        )
        return summarize_explain_results(results, dataset=target.dataset, collection=target.collection)
    return _summarize_existing_output(target, output_dir, module_ids)


def run_real_data_interpretation(
    scope: str = "all",
    fdr_threshold: float = 0.10,
    top_n_fallback: int = 5,
    all_modules: bool = False,
    annotation_table: Path | None = None,
    gtf_path: Path | None = DEFAULT_GTF_PATH,
    gtf_cache: Path | None = DEFAULT_GTF_CACHE,
    max_transcripts_per_gene: int = 6,
    max_pairs_per_gene: int = 15,
    force: bool = False,
    plot: bool = False,
) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for target in iter_targets(scope):
        summary = run_target(
            target,
            fdr_threshold=fdr_threshold,
            top_n_fallback=top_n_fallback,
            all_modules=all_modules,
            annotation_table=annotation_table,
            gtf_path=gtf_path,
            gtf_cache=gtf_cache,
            max_transcripts_per_gene=max_transcripts_per_gene,
            max_pairs_per_gene=max_pairs_per_gene,
            force=force,
            plot=plot,
        )
        if not summary.empty:
            parts.append(summary)
            collection_out = ensure_dir(rel("real_data", target.collection, "_m"))
            current = pd.concat(
                [part for part in parts if (part["collection"] == target.collection).all()],
                ignore_index=True,
            )
            current.to_parquet(
                collection_out / "module_interpret_summary.parquet",
                index=False,
                compression="zstd",
            )
    result = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if not result.empty:
        for collection, frame in result.groupby("collection"):
            out = ensure_dir(rel("real_data", collection, "_m"))
            frame.to_parquet(out / "module_interpret_summary.parquet", index=False, compression="zstd")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Run IsoGraph module interpretation on real-data modules.")
    parser.add_argument("--scope", choices=["all", "brainseq", "gtex"], default="all")
    parser.add_argument("--fdr-threshold", type=float, default=0.10)
    parser.add_argument("--top-n-fallback", type=int, default=5)
    parser.add_argument("--all-modules", action="store_true")
    parser.add_argument("--annotation-table", default=None)
    parser.add_argument("--gtf", default=str(DEFAULT_GTF_PATH))
    parser.add_argument("--gtf-cache", default=str(DEFAULT_GTF_CACHE))
    parser.add_argument("--no-gtf-annotation", action="store_true")
    parser.add_argument("--max-transcripts-per-gene", type=int, default=6)
    parser.add_argument("--max-pairs-per-gene", type=int, default=15)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    annotation_table = Path(args.annotation_table) if args.annotation_table else None
    gtf_path = None if args.no_gtf_annotation else Path(args.gtf)
    gtf_cache = None if args.no_gtf_annotation or not args.gtf_cache else Path(args.gtf_cache)
    result = run_real_data_interpretation(
        scope=args.scope,
        fdr_threshold=args.fdr_threshold,
        top_n_fallback=args.top_n_fallback,
        all_modules=args.all_modules,
        annotation_table=annotation_table,
        gtf_path=gtf_path,
        gtf_cache=gtf_cache,
        max_transcripts_per_gene=args.max_transcripts_per_gene,
        max_pairs_per_gene=args.max_pairs_per_gene,
        force=args.force,
        plot=args.plot,
    )
    print(f"Wrote {len(result):,} interpreted module summary rows")


if __name__ == "__main__":
    main()
