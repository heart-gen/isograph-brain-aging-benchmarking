"""Module sizes for every method in every cohort x region store, as plain text tables.

IsoGraph and WGCNA both produce very large modules in some regions (IsoGraph GTEx at resolution 5.0;
classical gene-level WGCNA as a rule), so module size is something to report and compare across methods,
not a property of one method. This writes one row per module and one summary row per fit.

Outputs (tab-separated, under ``02_module_discovery/_m/module_sizes/``):
    module_sizes.tsv         cohort, region, method, module_id, n_genes, size_rank, frac_of_assigned
    module_size_summary.tsv  cohort, region, method, n_modules, n_genes_assigned, largest, median, mean,
                             largest_frac, n_modules_ge_900, partition_sha256

``n_genes_assigned`` counts genes in a module; every method's ``modules.parquet`` lists assigned genes
only (WGCNA's unassigned grey genes are not rows). ``partition_sha256`` hashes the sorted
(gene_id, module_id) pairs, so a size table can be matched to the fit it was read from.

    python -m isograph_benchmark.real_data.module_sizes
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd

from isograph_benchmark.paths import COHORTS, ensure_dir, stage_out

METHODS = ("isograph_vae", "wgcna_gene", "wgcna_switch_only", "wgcna_multiplex")
GIANT_GENES = 900  # the phenotype-blind giant-module criterion behind resolution 5.0


def _partition_sha256(modules: pd.DataFrame) -> str:
    pairs = modules[["gene_id", "module_id"]].astype(str).sort_values(["gene_id", "module_id"])
    return hashlib.sha256(pairs.to_csv(index=False, header=False).encode()).hexdigest()


def collect(modules_root: Path, methods: tuple[str, ...] = METHODS) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-module and per-fit size tables for every ``<cohort>/<region>/_m/<method>/modules.parquet``."""
    rows, summary = [], []
    for cohort in COHORTS:
        for store in sorted((modules_root / cohort).glob("*/_m")):
            region = store.parent.name
            for method in methods:
                path = store / method / "modules.parquet"
                if not path.exists():
                    continue
                modules = pd.read_parquet(path, columns=["gene_id", "module_id"])
                sizes = (modules.groupby(modules["module_id"].astype(str)).size()
                         .sort_values(ascending=False, kind="stable"))
                n_assigned = int(sizes.sum())
                for rank, (module_id, n) in enumerate(sizes.items(), start=1):
                    rows.append({"cohort": cohort, "region": region, "method": method,
                                 "module_id": module_id, "n_genes": int(n), "size_rank": rank,
                                 "frac_of_assigned": round(n / n_assigned, 4)})
                summary.append({"cohort": cohort, "region": region, "method": method,
                                "n_modules": int(len(sizes)), "n_genes_assigned": n_assigned,
                                "largest": int(sizes.iloc[0]), "median": float(sizes.median()),
                                "mean": round(float(sizes.mean()), 1),
                                "largest_frac": round(float(sizes.iloc[0]) / n_assigned, 4),
                                "n_modules_ge_900": int((sizes >= GIANT_GENES).sum()),
                                "partition_sha256": _partition_sha256(modules)})
    return pd.DataFrame(rows), pd.DataFrame(summary)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=None, help="output directory (default: modules.sizes)")
    args = ap.parse_args()
    per_module, summary = collect(stage_out("modules"))
    if summary.empty:
        raise SystemExit("no modules.parquet found under 02_module_discovery/<cohort>/<region>/_m/")
    out = ensure_dir(args.out or stage_out("modules.sizes"))
    per_module.to_csv(out / "module_sizes.tsv", sep="\t", index=False)
    summary.to_csv(out / "module_size_summary.tsv", sep="\t", index=False)
    print(f"{len(per_module)} modules across {len(summary)} fits -> {out}")
    print(summary.drop(columns="partition_sha256").to_string(index=False))


if __name__ == "__main__":
    main()
