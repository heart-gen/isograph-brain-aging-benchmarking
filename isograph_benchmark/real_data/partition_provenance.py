"""Guard against joining a module table to an enrichment table from a DIFFERENT fit.

WHY THIS EXISTS
---------------
``module_enrichment`` keys its output on ``module_id``. Leiden module ids are
size-rank labels re-assigned on every fit: after a refit, ``M002`` is a different
gene set than it was before. Median same-id Jaccard between two fits of GTEx
cortex measured 0.000 — the ids are essentially a random relabelling.

Every downstream consumer joins the enrichment table to ``modules.parquet`` on
``module_id``. Nothing checked that the two came from the same fit, so a stale
enrichment table did not fail; it silently handed each module some OTHER module's
``pheno_fdr`` and ``n_go_terms``. That happened once, in commit 8376556
(2026-06-29): ``qtl_anchoring`` ran against a refreshed fit but the pre-refresh
enrichment table, and the resulting GO-invisible / GO-visible partition — and the
splicing-specificity headline built on it — was scrambled rather than merely out
of date. It was caught two months later only because the numbers were re-derived.

WHAT IS CHECKED
---------------
Two tiers, cheapest first:

1. ``partition_fingerprint`` — sha256 over the sorted ``(gene_id, module_id)``
   assignment. Written into the enrichment parquet's file metadata by
   ``module_enrichment`` and compared here. Exact, but only present on tables
   written after this guard landed.

2. Structural equivalence — the module-id SETS must be equal and the per-module
   ``n_genes`` must match the actual group sizes in ``modules.parquet``. This
   needs no regeneration, so it protects the tables already committed. On the
   2026-06-29 stale pairing it flags 29 of 33 modules; on the correct pairing it
   flags none, and it passes on all 66 committed (enrichment, fit) pairs.

Tier 2 alone is not airtight — a relabelling that happened to preserve every
group size would slip through — which is why tier 1 exists and why regenerating
the enrichment tables is worth doing. Tier 2 is the floor, not the ceiling.

Consumers should call :func:`load_enrichment` instead of ``pd.read_parquet``.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

#: Parquet file-metadata key holding the fit fingerprint.
FINGERPRINT_KEY = b"isograph_partition_sha256"


class StalePartitionError(RuntimeError):
    """An enrichment table does not describe the fit it is being joined to."""


def partition_fingerprint(modules: pd.DataFrame) -> str:
    """sha256 over the sorted ``(gene_id, module_id)`` assignment.

    Order-independent and index-independent, so it is stable across the different
    ways callers happen to load ``modules.parquet``.
    """
    pairs = (
        modules[["gene_id", "module_id"]]
        .astype(str)
        .sort_values(["gene_id", "module_id"], kind="mergesort")
        .itertuples(index=False, name=None)
    )
    h = hashlib.sha256()
    for gene, module in pairs:
        h.update(gene.encode())
        h.update(b"\t")
        h.update(module.encode())
        h.update(b"\n")
    return h.hexdigest()


def read_fingerprint(path: Path) -> str | None:
    """The fingerprint stamped into a parquet file, or ``None`` for older files."""
    meta = pq.read_schema(path).metadata or {}
    raw = meta.get(FINGERPRINT_KEY)
    return raw.decode() if raw is not None else None


def write_with_fingerprint(df: pd.DataFrame, path: Path, modules: pd.DataFrame) -> str:
    """Write ``df`` to parquet, stamping the fingerprint of ``modules`` into it."""
    import pyarrow as pa

    sha = partition_fingerprint(modules)
    table = pa.Table.from_pandas(df, preserve_index=False)
    meta = {**(table.schema.metadata or {}), FINGERPRINT_KEY: sha.encode()}
    pq.write_table(table.replace_schema_metadata(meta), path)
    return sha


def check_partition(
    modules: pd.DataFrame,
    enrich: pd.DataFrame,
    *,
    context: str,
    enrich_path: Path | None = None,
) -> None:
    """Raise :class:`StalePartitionError` unless ``enrich`` describes ``modules``.

    ``context`` names the caller and region so the message says which re-run is
    needed, e.g. ``"qtl_anchoring gtex-aging/cortex [isograph]"``.
    """
    fit_ids = modules["module_id"].astype(str)
    enrich_ids = enrich["module_id"].astype(str)

    # Tier 1: exact fingerprint, when the table carries one.
    if enrich_path is not None:
        stamped = read_fingerprint(enrich_path)
        if stamped is not None and stamped != partition_fingerprint(modules):
            raise StalePartitionError(
                f"{context}: enrichment table {enrich_path} was built from a "
                f"different fit (fingerprint mismatch). Re-run module_enrichment "
                f"for this region, then re-run this analysis."
            )

    # Tier 2: structural equivalence, which works on pre-fingerprint tables.
    only_enrich = sorted(set(enrich_ids) - set(fit_ids))
    only_fit = sorted(set(fit_ids) - set(enrich_ids))
    if only_enrich or only_fit:
        raise StalePartitionError(
            f"{context}: module-id sets differ between the fit and its enrichment "
            f"table ({len(only_fit)} module(s) in the fit with no enrichment row: "
            f"{only_fit[:5]}; {len(only_enrich)} enrichment row(s) naming a module "
            f"absent from the fit: {only_enrich[:5]}). The enrichment table is "
            f"stale — re-run module_enrichment for this region."
        )

    if "n_genes" in enrich.columns:
        actual = fit_ids.value_counts()
        stated = enrich.set_index(enrich_ids)["n_genes"]
        mismatched = stated.index[stated.values != stated.index.map(actual).values]
        if len(mismatched):
            raise StalePartitionError(
                f"{context}: {len(mismatched)} of {len(stated)} modules disagree on "
                f"gene count between the fit and its enrichment table "
                f"(e.g. {sorted(mismatched)[:5]}). Module ids are re-assigned on "
                f"every fit, so this join would silently scramble pheno_fdr and "
                f"n_go_terms across modules. Re-run module_enrichment."
            )


#: ``module_enrichment`` table prefix -> the fit directory it describes.
ENRICH_PREFIX_FIT_DIR = {
    "isograph": "isograph_vae",
    "wgcna": "wgcna_gene",
    "wgcna_switch": "wgcna_switch_only",
    "wgcna_multiplex": "wgcna_multiplex",
}


def load_region_enrichment(
    cohort: str, region: str, prefix: str, *, context: str
) -> pd.DataFrame | None:
    """Load ``<region>/_m/module_enrichment/<prefix>_modules.parquet``, verified.

    For the consumers that address the store by (cohort, region, method) rather
    than by an already-open fit directory. Returns ``None`` when the enrichment
    table is absent, matching those callers' existing "no data for this region"
    handling. A missing FIT cannot be verified against, so it warns loudly rather
    than passing silently — an unverifiable join is the state this module exists
    to make visible.
    """
    from isograph_benchmark.paths import region_store

    enrich_path = region_store(cohort, region, "module_enrichment", f"{prefix}_modules.parquet")
    if not enrich_path.exists():
        return None
    fit_dir = ENRICH_PREFIX_FIT_DIR.get(prefix)
    fit_path = region_store(cohort, region, fit_dir, "modules.parquet") if fit_dir else None
    if fit_path is None or not fit_path.exists():
        print(
            f"WARNING [{context}]: cannot verify {enrich_path.name} against a fit "
            f"({fit_path} absent) — module-id join is UNVERIFIED.",
            flush=True,
        )
        return pd.read_parquet(enrich_path)
    return load_enrichment(
        enrich_path, pd.read_parquet(fit_path), context=context, required=True
    )


def load_enrichment(
    enrich_path: Path,
    modules: pd.DataFrame,
    *,
    context: str,
    required: bool = False,
) -> pd.DataFrame | None:
    """Load an enrichment table and verify it describes ``modules``.

    Returns ``None`` when the table is absent and ``required`` is False — several
    consumers treat a missing table as "this region has no module-set breakdown"
    and degrade gracefully. A table that is PRESENT but stale always raises: that
    is the failure this module exists to stop.
    """
    enrich_path = Path(enrich_path)
    if not enrich_path.exists():
        if required:
            raise FileNotFoundError(f"{context}: missing enrichment table {enrich_path}")
        return None
    enrich = pd.read_parquet(enrich_path)
    check_partition(modules, enrich, context=context, enrich_path=enrich_path)
    return enrich
