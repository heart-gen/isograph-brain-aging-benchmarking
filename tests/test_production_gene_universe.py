"""The gene universe both methods must share, and the R baselines' contract with it.

The classical gene-level WGCNA baseline is fit in R and cannot call the production
transcript filter, so it reads the universe this module writes. If the two drift apart the
IsoGraph-vs-WGCNA head-to-heads compare module sets drawn from different gene pools.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.real_data import production_gene_universe as pgu
from isograph_benchmark.real_data import run_models

R_BASELINES = [
    "02_module_discovery/_h/01d.wgcna_gene_brainseq_aging.R",
    "02_module_discovery/_h/01e.wgcna_gene_brainseq_sczd.R",
    "02_module_discovery/_h/01f.wgcna_gene_gtex.R",
]


def _toy():
    n = 10
    rows = {
        "t1": ("A", [100] * n), "t2": ("A", [50, 50] + [0] * (n - 2)), "t3": ("A", [5] * n),
        "t6": ("B", [1] * n), "t7": ("B", [2] * n),
        "t4": ("C", [1000] * n), "t5": ("C", [20] * n),
    }
    table = pd.DataFrame({"transcript_id": list(rows), "gene_id": [g for g, _ in rows.values()]})
    counts = np.array([c for _, c in rows.values()], dtype=float)
    return counts, table


def test_universe_is_the_genes_surviving_the_production_filter():
    counts, table = _toy()
    _, kept = run_models.filter_production_transcripts(counts, table)
    # gene B never clears the gene-level floor, so it is not in the universe
    assert set(kept["gene_id"]) == {"A", "C"}
    per_gene = kept.groupby("gene_id").size().to_dict()
    assert per_gene == {"A": 2, "C": 1}


def test_collections_cover_every_production_region():
    covered = {(c, r) for c, r, _, _ in pgu.COLLECTIONS}
    expected = (
        {("brainseq", r) for r in pgu.BRAINSEQ_AGING_REGIONS}
        | {("brainseq", "caudate_sczd")}
        | {("gtex", r) for r in run_models.GTEX_REGIONS}
    )
    assert covered == expected
    # the SCZD store dir and its bundle region differ; the mapping must not be lost
    sczd = [c for c in pgu.COLLECTIONS if c[1] == "caudate_sczd"]
    assert sczd == [("brainseq", "caudate_sczd", "brainseq_sczd", "caudate")]


def test_universe_written_beside_the_backend_dirs_not_inside_one():
    # The R baselines read dirname(out_dir) -- the region _m dir, one level above
    # wgcna_gene/ -- so the universe must sit there and not inside any backend dir.
    from isograph_benchmark.paths import region_store

    out = region_store("brainseq", "caudate") / "production_gene_universe.parquet"
    assert out.parent.name == "_m"
    assert out.parts[-4:] == ("brainseq", "caudate", "_m", "production_gene_universe.parquet")


def test_r_baselines_read_the_universe_and_never_the_bundle_gene_list():
    for path in R_BASELINES:
        text = Path(path).read_text()
        assert "production_gene_universe.parquet" in text, path
        # the bundle's genes.parquet is the wider list this replaced
        assert not re.search(r'file\.path\(bundle_dir,\s*"genes\.parquet"\)', text), path
        # a missing universe must stop the fit, never silently fall back
        assert "missing gene universe" in text, path


def test_stage_dag_runs_the_universe_before_the_wgcna_gene_fits():
    dag = Path("02_module_discovery/_h/run_stage.sh").read_text()
    assert re.search(r'^step 00a ""\s+\$H/00a\.production_gene_universe\.sh', dag, re.M)
    for sid in ("01d", "01e", "01f"):
        assert re.search(rf'^step {sid} "00a"', dag, re.M), sid
