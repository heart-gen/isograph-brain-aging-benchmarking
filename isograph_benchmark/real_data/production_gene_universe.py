"""Emit the production gene universe for each cohort x region.

Every production fit applies ``run_models.filter_production_transcripts`` to the bundle
transcript matrix; the genes that survive it are the universe IsoGraph models. The classical
gene-level WGCNA baseline (``02_module_discovery/_h/01d,01e,01f.wgcna_gene_*.R``) builds its
expression matrix in R and cannot call that filter, so before 2026-09-15 it took the bundle's
full gene list instead -- a wider universe than IsoGraph's, which made the two module sets
non-comparable in every downstream head-to-head (MAGMA GSA, module trust, replication).

This CLI writes that universe once, as an artifact both methods read:

    02_module_discovery/<cohort>/<region>/_m/production_gene_universe.parquet
    02_module_discovery/<cohort>/<region>/_m/production_gene_universe.json

The filter itself stays defined in exactly one place (``filter_production_transcripts``); this
module only records what it returns, so an R consumer never re-implements the rule.
"""
from __future__ import annotations

import argparse
import json

import pandas as pd

from isograph.io.artifacts import load_dataset_bundle

from isograph_benchmark.paths import ensure_dir, region_store, rel
from isograph_benchmark.real_data.run_models import (
    GTEX_REGIONS,
    PRODUCTION_TRANSCRIPT_FILTER_LABEL,
    filter_production_transcripts,
)

BRAINSEQ_AGING_REGIONS = ["caudate", "hippocampus", "dlpfc"]

#: (cohort, region, bundle store, bundle region) -- the region dir a fit writes to is not always
#: the bundle region: the SCZD caudate bundle lands in the ``caudate_sczd`` store.
COLLECTIONS: list[tuple[str, str, str, str]] = (
    [("brainseq", r, "brainseq_v1", r) for r in BRAINSEQ_AGING_REGIONS]
    + [("brainseq", "caudate_sczd", "brainseq_sczd", "caudate")]
    + [("gtex", r, "gtex_v11_brain", r) for r in GTEX_REGIONS]
)


def build_one(cohort: str, region: str, store: str, bundle_region: str) -> pd.DataFrame:
    bundle = load_dataset_bundle(rel("inputs", "bundles", store, bundle_region))
    tx_table = bundle.feature_tables["transcript"]
    n_genes_bundle = int(tx_table["gene_id"].nunique())
    _, kept = filter_production_transcripts(bundle.matrices["transcript_counts"], tx_table)

    universe = (
        kept.groupby("gene_id", sort=True)
        .size()
        .rename("n_transcripts_kept")
        .reset_index()
    )
    universe["gene_id"] = universe["gene_id"].astype(str)

    out = ensure_dir(region_store(cohort, region))
    universe.to_parquet(out / "production_gene_universe.parquet", index=False)
    (out / "production_gene_universe.json").write_text(
        json.dumps(
            {
                "cohort": cohort,
                "region": region,
                "bundle": f"{store}/{bundle_region}",
                "filter": PRODUCTION_TRANSCRIPT_FILTER_LABEL,
                "n_genes_bundle": n_genes_bundle,
                "n_genes_kept": int(len(universe)),
                "n_transcripts_kept": int(len(kept)),
                "n_multi_transcript_genes": int((universe["n_transcripts_kept"] > 1).sum()),
            },
            indent=2,
        )
        + "\n"
    )
    print(
        f"[{cohort}/{region}] {len(universe)}/{n_genes_bundle} genes kept "
        f"-> {out / 'production_gene_universe.parquet'}",
        flush=True,
    )
    return universe


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort", choices=["brainseq", "gtex"], help="limit to one cohort")
    ap.add_argument("--region", help="limit to one region (store dir name)")
    args = ap.parse_args()

    todo = [
        c for c in COLLECTIONS
        if (args.cohort is None or c[0] == args.cohort)
        and (args.region is None or c[1] == args.region)
    ]
    if not todo:
        raise SystemExit(f"no collection matches cohort={args.cohort} region={args.region}")
    for cohort, region, store, bundle_region in todo:
        build_one(cohort, region, store, bundle_region)


if __name__ == "__main__":
    main()
