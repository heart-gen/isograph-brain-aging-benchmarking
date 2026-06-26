"""Characterize the composition-unique gene sets from the de-confounded test.

For genes whose isoform COMPOSITION is phenotype-associated beyond their own
abundance (category == 'composition_unique' in incremental_association/
gene_level.parquet), describe what they are — this is IsoGraph's complementary,
otherwise-invisible contribution (DTU-without-DGE):

  - go_terms.parquet           : GO:BP enrichment of the set vs the tested-gene
                                 background.
  - module_concentration.parquet: hypergeometric enrichment of the set in each
                                 IsoGraph module (do these genes organize into
                                 specific switch modules?).
  - genes.parquet              : the composition-unique genes + their stats.

--combine reads the per-analysis gene lists and writes the cross-analysis
overlap (e.g. SCZD ∩ caudate-aging) to real_data/brainseq/_m/.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.go_enrichment import GoAnnotations, HAS_GOATOOLS
from isograph_benchmark.real_data.sweep_leiden import _artifact_dir

_DEFAULT_GO_CACHE = rel("inputs", "go_annotations")
_ANALYSES = [
    ("brainseq-sczd", None),
    ("brainseq-aging", "caudate"),
    ("brainseq-aging", "hippocampus"),
    ("brainseq-aging", "dlpfc"),
]


def _label(analysis: str, region: str | None) -> str:
    return f"{analysis}:{region}" if region else analysis


def _wait_for_upstream(path, label: str, timeout_s: float = 600.0, poll_s: float = 10.0):
    """Block until an upstream parquet materializes, absorbing afterok launch skew.

    This stage consumes incremental_association/gene_level.parquet. If the SLURM
    DAG is wired to depend on the wrong parent (afterok on 04.interpret rather than
    08.incremental) or simply races a still-flushing writer, a bare existence check
    crashes the whole array. Poll up to ``timeout_s`` so a small ordering skew is
    tolerated; otherwise fail loudly with the correct dependency named.
    """
    deadline = time.monotonic() + timeout_s
    waited = False
    while not (path.exists() and path.stat().st_size > 0):
        if time.monotonic() >= deadline:
            raise FileNotFoundError(
                f"[{label}] missing {path} after waiting {timeout_s:.0f}s; run "
                "incremental_association first (wire afterok on 08.incremental, "
                "not 04.interpret)."
            )
        if not waited:
            print(f"[{label}] waiting for upstream {path.name} (afterok skew)...", flush=True)
            waited = True
        time.sleep(poll_s)
    if waited:
        print(f"[{label}] upstream {path.name} appeared; proceeding", flush=True)


def _module_concentration(unique_genes: set[str], background: set[str], modules: pd.DataFrame) -> pd.DataFrame:
    """Hypergeometric enrichment of the composition-unique set per IsoGraph module."""
    modules = modules[modules["gene_id"].isin(background)].copy()
    N = len(background)
    K = len(unique_genes & background)
    rows = []
    for mid, grp in modules.groupby("module_id"):
        mod_genes = set(grp["gene_id"])
        n = len(mod_genes)
        k = len(mod_genes & unique_genes)
        if k == 0:
            continue
        p = stats.hypergeom.sf(k - 1, N, K, n)
        rows.append({"module_id": mid, "module_size": n, "n_composition_unique": k,
                     "frac_module": round(k / n, 4), "hyperg_p": p})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["hyperg_fdr"] = stats.false_discovery_control(out["hyperg_p"], method="bh")
        out = out.sort_values("hyperg_p").reset_index(drop=True)
    return out


def characterize(analysis: str, region: str | None, variant: str,
                 go_cache_dir=None, skip_go: bool = False) -> dict:
    label = _label(analysis, region)
    artifact_dir = _artifact_dir(analysis, region, variant)
    inc_dir = artifact_dir / "incremental_association"
    gl_path = inc_dir / "gene_level.parquet"
    _wait_for_upstream(gl_path, label)
    gene_level = pd.read_parquet(gl_path)
    modules = pd.read_parquet(artifact_dir / "modules.parquet")

    background = set(gene_level["gene_id"])
    uniq = gene_level[gene_level["category"] == "composition_unique"].copy()
    unique_set = set(uniq["gene_id"])
    print(f"[{label}] composition-unique genes: {len(unique_set)} / {len(background)} tested", flush=True)

    out = ensure_dir(inc_dir.parent / "composition_unique")
    uniq.to_parquet(out / "genes.parquet", index=False, compression="zstd")

    conc = _module_concentration(unique_set, background, modules)
    conc.to_parquet(out / "module_concentration.parquet", index=False, compression="zstd")

    go = pd.DataFrame()
    if not skip_go and HAS_GOATOOLS and len(unique_set) >= 5:
        try:
            helper = GoAnnotations(go_cache_dir or _DEFAULT_GO_CACHE)
            helper.prepare(sorted(background))
            go = helper.enrich_gene_set(sorted(unique_set))
        except Exception as exc:  # pragma: no cover - GO is optional
            print(f"[{label}] WARNING: GO enrichment failed ({exc})", flush=True)
    go.to_parquet(out / "go_terms.parquet", index=False, compression="zstd")

    n_conc_sig = int((conc["hyperg_fdr"] <= 0.10).sum()) if "hyperg_fdr" in conc.columns else 0
    summary = {
        "analysis": analysis, "region": region, "variant": variant,
        "n_composition_unique": len(unique_set), "n_tested": len(background),
        "n_modules_enriched_fdr10": n_conc_sig, "n_go_terms": int(len(go)),
        "top_go_terms": go.head(10)["term_name"].tolist() if not go.empty else [],
        "top_concentrated_modules": conc.head(5)["module_id"].tolist() if not conc.empty else [],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[{label}] modules enriched (FDR<=.10)={n_conc_sig} | GO:BP terms={len(go)}")
    if not go.empty:
        print(f"[{label}] top GO: {', '.join(go.head(5)['term_name'])}")
    print(f"[{label}] written to {out}")
    return summary


def combine(variant: str) -> None:
    """Cross-analysis overlap of composition-unique gene sets."""
    sets = {}
    for analysis, region in _ANALYSES:
        p = _artifact_dir(analysis, region, variant) / "composition_unique" / "genes.parquet"
        if p.exists():
            sets[_label(analysis, region)] = set(pd.read_parquet(p)["gene_id"])
    labels = list(sets)
    rows = []
    for i, a in enumerate(labels):
        for b in labels[i:]:
            inter = sets[a] & sets[b]
            rows.append({"set_a": a, "set_b": b, "n_a": len(sets[a]), "n_b": len(sets[b]),
                         "n_overlap": len(inter),
                         "jaccard": round(len(inter) / len(sets[a] | sets[b]), 4) if (sets[a] | sets[b]) else 0.0})
    out = ensure_dir(rel("real_data", "brainseq", "_m"))
    df = pd.DataFrame(rows)
    df.to_parquet(out / "composition_unique_overlap.parquet", index=False, compression="zstd")
    print("Cross-analysis composition-unique overlap:")
    print(df.to_string(index=False))
    print(f"written to {out / 'composition_unique_overlap.parquet'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis", nargs="?", choices=["brainseq-sczd", "brainseq-aging"])
    parser.add_argument("--region", action="append", dest="regions")
    parser.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    parser.add_argument("--no-go", action="store_true")
    parser.add_argument("--combine", action="store_true",
                        help="Cross-analysis overlap (ignores positional analysis).")
    args = parser.parse_args()

    if args.combine:
        combine(args.variant)
        return
    if args.analysis == "brainseq-sczd":
        characterize("brainseq-sczd", None, args.variant, skip_go=args.no_go)
    elif args.analysis == "brainseq-aging":
        for region in (args.regions or ["caudate", "hippocampus", "dlpfc"]):
            characterize("brainseq-aging", region, args.variant, skip_go=args.no_go)
    else:
        parser.error("provide an analysis or --combine")


if __name__ == "__main__":
    main()
