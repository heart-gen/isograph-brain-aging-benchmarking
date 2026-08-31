"""Cross-cohort GO consistency for replicated, age-associated modules.

Extends the BrainSEQ-vs-GTEx replication (replication.py) with the functional
layer: for module pairs that are both cross-cohort *preserved* and *age-associated*,
do the two cohorts' matched modules enrich for the SAME biology?

For each preserved + age-associated source module (best-match target in the other
cohort), we compare the full enriched GO:BP term sets of the source and target
modules (persisted by module_enrichment.py as ``<method>_module_go.parquet``):

* ``go_jaccard`` = |shared terms| / |union of terms| over enriched BP term IDs,
* the shared term names, and
* a pooled permutation test: is the mean GO Jaccard of the preserved+aging matched
  pairs higher than when each source module is paired with a *random* target
  module from the same cohort/region? (1,000 permutations -> empirical p.)

"Age-associated" uses the eigengene-age linear FDR carried in the replication
match table (``source_age_fdr``), which is available for both methods and both
cohorts.

Usage::

    python -m isograph_benchmark.real_data.replication_go            # both methods
    python -m isograph_benchmark.real_data.replication_go --methods isograph_vae

Outputs (under ``03_module_trust/_m/replication/``):
    <method>_go_overlap.parquet   per preserved+aging pair: GO overlap + shared terms
    replication_go_summary.parquet / .json   per-method pooled test
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from isograph_benchmark.real_data.partition_provenance import load_region_enrichment
from isograph_benchmark.paths import ensure_dir, region_store, stage_out
from isograph_benchmark.real_data.replication import REGION_PAIRS

DEFAULT_METHODS = ["isograph_vae", "wgcna_gene"]
AGE_FDR_THRESHOLD = 0.10
N_PERMUTATIONS = 1000
SEED = 13

# replication method name -> module_enrichment file prefix
_PREFIX = {
    "isograph_vae": "isograph",
    "wgcna_gene": "wgcna",
    "wgcna_switch_only": "wgcna_switch",
    "wgcna_multiplex": "wgcna_multiplex",
}


def _region_map() -> dict[str, tuple[str, str]]:
    return {label: (bs, gt) for bs, gt, label in REGION_PAIRS}


def _go_sets(cohort: str, region: str, method: str) -> dict[str, set]:
    """module_id -> set of enriched BP term IDs (empty set for modules w/o terms)."""
    prefix = _PREFIX[method]
    base = region_store(cohort, region, "module_enrichment")
    go_path = base / f"{prefix}_module_go.parquet"
    mod_path = base / f"{prefix}_modules.parquet"
    if not go_path.exists() or not mod_path.exists():
        return {}
    all_modules = load_region_enrichment(
        cohort, region, prefix, context=f"replication_go {cohort}/{region} [{method}]"
    )["module_id"].unique()
    sets: dict[str, set] = {m: set() for m in all_modules}
    go = pd.read_parquet(go_path)
    if not go.empty:
        for mid, grp in go.groupby("module_id"):
            sets[mid] = set(grp["term_id"])
    return sets


def _jaccard(a: set, b: set) -> float:
    u = a | b
    return len(a & b) / len(u) if u else np.nan


def run_method(method: str, rng: np.random.Generator) -> tuple[pd.DataFrame, dict]:
    rmap = _region_map()
    match_path = stage_out("trust.replication", f"{method}_module_match.parquet")
    if not match_path.exists():
        print(f"  [{method}] no module_match table — run replication.py first; skipping")
        return pd.DataFrame(), {}
    match = pd.read_parquet(match_path)
    aging = match[match["preserved"] & (match["source_age_fdr"] < AGE_FDR_THRESHOLD)].copy()

    # GO-set caches keyed by (cohort, region)
    cache: dict[tuple[str, str], dict[str, set]] = {}

    def go_for(cohort, region):
        key = (cohort, region)
        if key not in cache:
            cache[key] = _go_sets(cohort, region, method)
        return cache[key]

    rows = []
    for _, r in aging.iterrows():
        bs, gt = rmap[r["region"]]
        if r["direction"] == "brainseq_to_gtex":
            sc, sr, tc, tr = "brainseq", bs, "gtex", gt
        else:
            sc, sr, tc, tr = "gtex", gt, "brainseq", bs
        src = go_for(sc, sr).get(r["source_module"], set())
        tgt = go_for(tc, tr).get(r["best_target_module"], set())
        jac = _jaccard(src, tgt)
        rows.append({
            "method": method, "region": r["region"], "direction": r["direction"],
            "source_module": r["source_module"], "target_module": r["best_target_module"],
            "module_jaccard": float(r["jaccard"]), "source_age_fdr": float(r["source_age_fdr"]),
            "age_sign_concordant": r["age_sign_concordant"],
            "n_source_go": len(src), "n_target_go": len(tgt),
            "n_shared_go": len(src & tgt), "go_jaccard": jac,
            "shared_go_terms": sorted(_shared_names(method, sc, sr, src & tgt)),
        })
    per_pair = pd.DataFrame(rows)

    # ---- pooled permutation test over pairs with a defined go_jaccard ----------
    usable = per_pair.dropna(subset=["go_jaccard"])
    summary = {
        "method": method,
        "n_preserved_aging_pairs": int(len(per_pair)),
        "n_pairs_with_go": int(len(usable)),
        "mean_go_jaccard": float(usable["go_jaccard"].mean()) if len(usable) else np.nan,
        "median_go_jaccard": float(usable["go_jaccard"].median()) if len(usable) else np.nan,
        "perm_p": np.nan,
    }
    if len(usable) >= 3:
        obs = usable["go_jaccard"].mean()
        null = np.empty(N_PERMUTATIONS)
        for k in range(N_PERMUTATIONS):
            vals = []
            for _, r in usable.iterrows():
                bs, gt = rmap[r["region"]]
                tc, tr = ("gtex", gt) if r["direction"] == "brainseq_to_gtex" else ("brainseq", bs)
                sc, sr = ("brainseq", bs) if r["direction"] == "brainseq_to_gtex" else ("gtex", gt)
                pool = list(go_for(tc, tr).values())
                src = go_for(sc, sr).get(r["source_module"], set())
                rand_tgt = pool[rng.integers(len(pool))] if pool else set()
                vals.append(_jaccard(src, rand_tgt))
            null[k] = np.nanmean(vals) if np.any(np.isfinite(vals)) else np.nan
        valid_null = null[np.isfinite(null)]
        if len(valid_null):
            summary["perm_p"] = float((1 + np.sum(valid_null >= obs)) / (len(valid_null) + 1))
            summary["null_mean_go_jaccard"] = float(valid_null.mean())
    return per_pair, summary


def _shared_names(method, cohort, region, term_ids: set) -> list[str]:
    """Resolve shared term IDs to names from the source cohort's GO table."""
    if not term_ids:
        return []
    prefix = _PREFIX[method]
    go_path = region_store(cohort, region, "module_enrichment", f"{prefix}_module_go.parquet")
    if not go_path.exists():
        return [str(t) for t in term_ids]
    go = pd.read_parquet(go_path).drop_duplicates("term_id").set_index("term_id")["term_name"]
    return [str(go.get(t, t)) for t in term_ids]


def main() -> None:
    ap = argparse.ArgumentParser(description="Cross-cohort GO consistency for replicated aging modules.")
    ap.add_argument("--methods", nargs="+", default=DEFAULT_METHODS)
    args = ap.parse_args()

    out_dir = ensure_dir(stage_out("trust.replication"))
    rng = np.random.default_rng(SEED)

    summaries = []
    for method in args.methods:
        print(f"GO consistency for {method} ...")
        per_pair, summ = run_method(method, rng)
        if per_pair.empty:
            continue
        per_pair.to_parquet(out_dir / f"{method}_go_overlap.parquet", index=False, compression="zstd")
        summaries.append(summ)
        print(f"  {summ['n_preserved_aging_pairs']} preserved+aging pairs | "
              f"{summ['n_pairs_with_go']} with GO | mean GO Jaccard "
              f"{summ['mean_go_jaccard']:.3f} | perm p {summ['perm_p']}")

    if not summaries:
        raise SystemExit("No GO-overlap results — are the module_go enrichment tables present?")
    sdf = pd.DataFrame(summaries)
    sdf.to_parquet(out_dir / "replication_go_summary.parquet", index=False, compression="zstd")
    (out_dir / "replication_go_summary.json").write_text(json.dumps(summaries, indent=2, default=str))
    print(f"\nWrote {out_dir}/replication_go_summary.parquet")
    print(sdf.to_string(index=False))


if __name__ == "__main__":
    main()
