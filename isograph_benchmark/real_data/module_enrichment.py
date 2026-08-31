"""Module-level GO enrichment + network metrics — the module-network narrative.

IsoGraph is a network method, so the story is told at the MODULE level: what is
each module's biology (GO:BP), how is it wired (intramodular hubs/connectivity),
and is it phenotype-associated. This runs per-module GO enrichment for BOTH
IsoGraph and WGCNA partitions (a fair method comparison: how much of each method's
module structure is functionally coherent), plus network metrics for IsoGraph
(WGCNA's saved artifacts are partition-only, no edge list).

Per analysis it writes, under <region>/_m/module_enrichment/:
  isograph_modules.parquet / wgcna_modules.parquet
      one row per module: n_genes, n_go_terms, top_go_terms, min_go_fdr,
      phenotype fdr (joined from the method's trait/diagnosis output), and for
      IsoGraph the network metrics (n_intra_edges, intra_density, mean_abs_weight,
      frac_pos_edges, hub_genes, hub_degrees).
  summary.json — method comparison (fraction of modules GO-enriched, etc.).

Complements IsoGraph's built-in transcript-level interpretation
(interpret_modules.py / 04.interpret_modules.sh), which explains the switch
events inside the phenotype-significant modules.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir
from isograph_benchmark.real_data.partition_provenance import write_with_fingerprint
from isograph_benchmark.real_data.go_enrichment import GoAnnotations, HAS_GOATOOLS
from isograph_benchmark.real_data.sweep_leiden import _DEFAULT_GO_CACHE, _artifact_dir

_TOP_GO = 8
_TOP_HUBS = 10


METHOD_DIRS = {
    "wgcna": "wgcna_gene",
    "wgcna_switch": "wgcna_switch_only",
    "wgcna_multiplex": "wgcna_multiplex",
}


def _method_dir(analysis: str, region: str | None, method: str, variant: str):
    """Artifact dir for a method. IsoGraph honours --variant; WGCNA dirs are explicit."""
    iso_dir = _artifact_dir(analysis, region, variant)   # .../_m/isograph_vae[_with_abundance]
    if method == "isograph":
        return iso_dir
    if method in METHOD_DIRS:
        return iso_dir.parent / METHOD_DIRS[method]
    raise ValueError(f"unknown method {method!r}")


def _phenotype_fdr(analysis: str, method_dir) -> pd.Series:
    """module_id -> phenotype FDR from the method's saved association table."""
    if analysis == "brainseq-sczd":
        f = method_dir / "diagnosis_assoc.parquet"
        col = "fdr"
    else:
        f = method_dir / "age_spline.parquet"
        col = "fdr_ftest"
    if not f.exists():
        return pd.Series(dtype=float)
    d = pd.read_parquet(f)
    if "module_id" not in d.columns or col not in d.columns:
        return pd.Series(dtype=float)
    return d.drop_duplicates("module_id").set_index("module_id")[col]


def _network_metrics(edges: pd.DataFrame, modules: pd.DataFrame) -> dict:
    """Per-module intramodular network metrics from the IsoGraph edge list."""
    gene_mod = dict(zip(modules["gene_id"], modules["module_id"]))
    e = edges[edges["weight"] != 0].copy()
    e["ms"] = e["source"].map(gene_mod)
    e["mt"] = e["target"].map(gene_mod)
    intra = e[(e["ms"].notna()) & (e["ms"] == e["mt"])]
    out: dict = {}
    for mid, grp in intra.groupby("ms"):
        n = int((modules["module_id"] == mid).sum())
        deg: dict[str, int] = {}
        for a, b in zip(grp["source"], grp["target"]):
            deg[a] = deg.get(a, 0) + 1
            deg[b] = deg.get(b, 0) + 1
        hubs = sorted(deg, key=lambda g: -deg[g])[:_TOP_HUBS]
        n_edges = len(grp)
        out[mid] = {
            "n_intra_edges": n_edges,
            "intra_density": round(2 * n_edges / (n * (n - 1)), 5) if n > 1 else 0.0,
            "mean_abs_weight": round(float(grp["weight"].abs().mean()), 4),
            "frac_pos_edges": round(float((grp["weight"] > 0).mean()), 4),
            "hub_genes": hubs,
            "hub_degrees": [deg[g] for g in hubs],
        }
    return out


_GO_LONG_COLS = ["module_id", "term_id", "term_name", "p_value", "p_fdr_bh", "study_count", "study_n"]


def _module_go(modules: pd.DataFrame, helper: GoAnnotations | None) -> tuple[dict, pd.DataFrame]:
    """Per-module GO summary plus the full enriched-term long table.

    Returns (summary_by_module, full_go_long).  ``summary_by_module`` keeps the
    compact per-module row (n_genes, n_go_terms = full enriched count, top_go_terms
    capped for readability, min_go_fdr).  ``full_go_long`` has one row per
    (module, enriched BP term) so downstream analyses (e.g. cross-cohort GO
    overlap) can use the complete enriched set rather than the top-N names.
    """
    out: dict = {}
    full_rows: list[pd.DataFrame] = []
    for mid, grp in modules.groupby("module_id"):
        genes = grp["gene_id"].tolist()
        rec = {"n_genes": len(genes), "n_go_terms": 0, "top_go_terms": [], "min_go_fdr": np.nan}
        if helper is not None:
            terms = helper.enrich_gene_set(genes, max_terms=None)   # full enriched BP set
            if not terms.empty:
                rec["n_go_terms"] = int(len(terms))
                rec["top_go_terms"] = terms["term_name"].tolist()[:_TOP_GO]
                rec["min_go_fdr"] = float(terms["p_fdr_bh"].min())
                t = terms.copy()
                t.insert(0, "module_id", mid)
                full_rows.append(t)
        out[mid] = rec
    full = (
        pd.concat(full_rows, ignore_index=True)[_GO_LONG_COLS]
        if full_rows else pd.DataFrame(columns=_GO_LONG_COLS)
    )
    return out, full


def characterize_method(analysis, region, method, variant, helper) -> pd.DataFrame | None:
    md = _method_dir(analysis, region, method, variant)
    mfile = md / "modules.parquet"
    if not mfile.exists():
        print(f"  [{method}] no modules.parquet at {md} — skipping")
        return None
    modules = pd.read_parquet(mfile)
    go, go_full = _module_go(modules, helper)
    rows = []
    pheno = _phenotype_fdr(analysis, md)
    net = {}
    if method == "isograph":
        epath = md / "edges.parquet"
        if epath.exists():
            net = _network_metrics(pd.read_parquet(epath), modules)
    for mid in sorted(modules["module_id"].unique()):
        row = {"method": method, "module_id": mid, **go[mid]}
        row["pheno_fdr"] = float(pheno.get(mid, np.nan))
        row.update(net.get(mid, {}))
        rows.append(row)
    df = pd.DataFrame(rows)
    out = ensure_dir(md.parent / "module_enrichment")
    # Stamp the fit this table describes. Module ids are re-assigned on every fit, so
    # a downstream join on module_id is only meaningful against THIS partition; the
    # fingerprint lets consumers reject a stale table instead of silently scrambling
    # pheno_fdr / n_go_terms across modules.
    write_with_fingerprint(df, out / f"{method}_modules.parquet", modules)
    go_full.to_parquet(out / f"{method}_module_go.parquet", index=False, compression="zstd")
    return df


def run_analysis(analysis, region, variant, methods, skip_go) -> dict:
    label = f"{analysis}/{region}" if region else analysis
    helper = None
    if not skip_go and HAS_GOATOOLS:
        # background = the IsoGraph gene universe (genes both methods are scored on)
        iso_modules = _artifact_dir(analysis, region, variant) / "modules.parquet"
        fs = _artifact_dir(analysis, region, variant) / "feature_scores.parquet"
        bg = pd.read_parquet(fs)["gene_id"].astype(str).unique().tolist() if fs.exists() \
            else pd.read_parquet(iso_modules)["gene_id"].astype(str).unique().tolist()
        helper = GoAnnotations(_DEFAULT_GO_CACHE)
        helper.prepare(sorted(bg))
    elif not skip_go:
        print(f"[{label}] goatools missing — GO skipped")

    frames, summary = [], {"analysis": analysis, "region": region, "variant": variant, "methods": {}}
    for method in methods:
        df = characterize_method(analysis, region, method, variant, helper)
        if df is None:
            continue
        frames.append(df)
        n = len(df)
        n_go = int((df["n_go_terms"] > 0).sum())
        n_ph = int((df["pheno_fdr"] <= 0.10).sum())
        summary["methods"][method] = {
            "n_modules": n, "n_go_enriched": n_go,
            "frac_go_enriched": round(n_go / n, 3) if n else 0.0,
            "n_pheno_sig_fdr10": n_ph,
        }
        print(f"[{label}] {method}: {n} modules | GO-enriched {n_go}/{n} "
              f"({summary['methods'][method]['frac_go_enriched']:.0%}) | pheno-sig {n_ph}")

    if frames:
        out = ensure_dir(_method_dir(analysis, region, "isograph", variant).parent / "module_enrichment")
        pd.concat(frames, ignore_index=True).to_parquet(out / "all_modules.parquet", index=False, compression="zstd")
        (out / "summary.json").write_text(json.dumps(summary, indent=2))
        print(f"[{label}] written to {out}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis", choices=["brainseq-sczd", "brainseq-aging", "gtex-aging"])
    parser.add_argument("--region", action="append", dest="regions")
    parser.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    parser.add_argument(
        "--method",
        choices=["isograph", "wgcna", "wgcna_switch", "wgcna_multiplex", "both", "all"],
        default="both",
    )
    parser.add_argument("--no-go", action="store_true")
    args = parser.parse_args()
    if args.method == "both":
        methods = ["isograph", "wgcna"]
    elif args.method == "all":
        methods = ["isograph", "wgcna", "wgcna_switch", "wgcna_multiplex"]
    else:
        methods = [args.method]

    if args.analysis == "brainseq-sczd":
        run_analysis("brainseq-sczd", None, args.variant, methods, args.no_go)
    elif args.analysis == "gtex-aging":
        from isograph_benchmark.real_data.run_models import GTEX_REGIONS
        for region in (args.regions or GTEX_REGIONS):
            run_analysis("gtex-aging", region, args.variant, methods, args.no_go)
    else:
        for region in (args.regions or ["caudate", "hippocampus", "dlpfc"]):
            run_analysis("brainseq-aging", region, args.variant, methods, args.no_go)


if __name__ == "__main__":
    main()
