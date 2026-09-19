"""Tests for the phenotype-blind resolution sweep on split-half partitions.

The sweep exists to answer a standing reviewer objection: the canonical Leiden resolution
was chosen while looking at a trait. Answering it needs the SAME split halves clustered
across a resolution grid, which must not disturb the committed canonical partitions and
must not cost one VAE fit per resolution.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data import stability as st


# --------------------------------------------------------------------------- #
# method tagging: the canonical arm must stay byte-identical
# --------------------------------------------------------------------------- #
def test_canonical_resolution_is_unsuffixed():
    """An empty tag is what keeps existing filenames, the resume-skip and the committed
    aggregates working untouched."""
    assert st._res_token(st.CANONICAL_LEIDEN_RESOLUTION) == ""


def test_noncanonical_resolutions_get_distinct_tags():
    tags = {st._res_token(r) for r in (1.0, 5.0, 2.5, 12.0)}
    assert tags == {"_res1", "_res5", "_res2p5", "_res12"}
    assert "" not in tags


def test_sweep_partitions_cannot_collide_with_canonical(tmp_path, monkeypatch):
    """Two resolutions of the same (cohort, region, seed, half) must be different files."""
    monkeypatch.setattr(st, "_partitions_dir", lambda: tmp_path)
    mods = pd.DataFrame({"gene_id": ["g1", "g2"], "module_id": ["M000", "M000"]})
    for res in (st.CANONICAL_LEIDEN_RESOLUTION, 5.0):
        st._write_partition(mods, "gtex", "cortex", "isograph" + st._res_token(res), 0, "A")
    written = sorted(p.name for p in tmp_path.glob("*.parquet"))
    assert written == [
        "isograph__gtex__cortex__seed0__A.parquet",
        "isograph_res5__gtex__cortex__seed0__A.parquet",
    ]


# --------------------------------------------------------------------------- #
# the resolution actually reaches the model config
# --------------------------------------------------------------------------- #
def test_resolution_reaches_the_vae_config():
    spec = st.COHORTS["gtex"]
    assert st._vae_config(spec).leiden_resolution == st.CANONICAL_LEIDEN_RESOLUTION
    assert st._vae_config(spec, leiden_resolution=2.0).leiden_resolution == 2.0


def test_config_default_matches_production():
    """The split-half baseline must cluster where production clusters, or the trust funnel
    validates a partition the paper does not ship."""
    from isograph_benchmark.real_data.run_models import CANONICAL_LEIDEN_RESOLUTION as prod
    assert st.CANONICAL_LEIDEN_RESOLUTION == prod


# --------------------------------------------------------------------------- #
# sweep over saved graphs
# --------------------------------------------------------------------------- #
def _planted_edges(n_modules: int = 3, per_module: int = 25) -> pd.DataFrame:
    """Disjoint cliques -- a graph with an unambiguous community structure."""
    rows = []
    for m in range(n_modules):
        genes = [f"g{m}_{i}" for i in range(per_module)]
        for i, a in enumerate(genes):
            for b in genes[i + 1:]:
                rows.append({"source": a, "target": b, "weight": 1.0})
    return pd.DataFrame(rows)


def test_sweep_reclusters_saved_graphs_without_refitting(tmp_path, monkeypatch):
    edges_dir = tmp_path / "edges"
    parts_dir = tmp_path / "partitions"
    edges_dir.mkdir(); parts_dir.mkdir()
    monkeypatch.setattr(st, "_edges_dir", lambda: edges_dir)
    monkeypatch.setattr(st, "_partitions_dir", lambda: parts_dir)

    edges = _planted_edges()
    for half in ("A", "B"):
        st._write_edges(edges, "gtex", "cortex", "isograph", 0, half)

    st.sweep("gtex", "cortex", [st.CANONICAL_LEIDEN_RESOLUTION, 1.0])

    names = sorted(p.name for p in parts_dir.glob("*.parquet"))
    # the canonical arm is the committed baseline; a sweep must not rewrite it
    assert names == [
        "isograph_res1__gtex__cortex__seed0__A.parquet",
        "isograph_res1__gtex__cortex__seed0__B.parquet",
    ]
    out = pd.read_parquet(parts_dir / names[0])
    assert out["module_id"].nunique() == 3, "planted cliques should recover as 3 modules"
    assert set(out.columns) >= {"gene_id", "module_id", "method", "cohort", "region"}
    assert out["method"].iloc[0] == "isograph_res1"


def test_sweep_without_saved_graphs_says_what_to_run(tmp_path, monkeypatch):
    monkeypatch.setattr(st, "_edges_dir", lambda: tmp_path)
    monkeypatch.setattr(st, "_partitions_dir", lambda: tmp_path)
    with pytest.raises(SystemExit, match="--save-edges"):
        st.sweep("gtex", "cortex", [2.0])


def test_sweep_resumes_and_does_not_rewrite(tmp_path, monkeypatch):
    edges_dir = tmp_path / "edges"
    parts_dir = tmp_path / "partitions"
    edges_dir.mkdir(); parts_dir.mkdir()
    monkeypatch.setattr(st, "_edges_dir", lambda: edges_dir)
    monkeypatch.setattr(st, "_partitions_dir", lambda: parts_dir)
    st._write_edges(_planted_edges(), "gtex", "cortex", "isograph", 0, "A")

    st.sweep("gtex", "cortex", [5.0])
    target = parts_dir / "isograph_res5__gtex__cortex__seed0__A.parquet"
    stamp = target.stat().st_mtime_ns
    st.sweep("gtex", "cortex", [5.0])
    assert target.stat().st_mtime_ns == stamp, "an existing partition was rewritten"


def test_saved_edges_keep_only_what_reclustering_needs(tmp_path, monkeypatch):
    """The graphs are a cache; carrying the annotation columns would multiply their size."""
    monkeypatch.setattr(st, "_edges_dir", lambda: tmp_path)
    fat = _planted_edges().head(10).assign(
        source_feature_id="x", target_feature_id="y",
        source_feature_type="switch", target_feature_type="switch")
    st._write_edges(fat, "gtex", "cortex", "isograph", 0, "A")
    back = pd.read_parquet(st._edges_path("gtex", "cortex", "isograph", 0, "A"))
    assert list(back.columns) == ["source", "target", "weight"]


# --------------------------------------------------------------------------- #
# resume must not silently no-op the edge-saving pass
# --------------------------------------------------------------------------- #
def test_resume_key_includes_the_edges_when_saving(tmp_path, monkeypatch):
    """--save-edges on an already-fitted region must still fit; the partition alone is not
    a complete record of what the invocation produces.

    Without this the sweep's prerequisite pass is a silent no-op on exactly the regions that
    matter -- the ones already fitted -- and its skip message reads like a normal resume.
    """
    parts = tmp_path / "partitions"; edges = tmp_path / "edges"
    parts.mkdir(); edges.mkdir()
    monkeypatch.setattr(st, "_partitions_dir", lambda: parts)
    monkeypatch.setattr(st, "_edges_dir", lambda: edges)

    st._write_partition(pd.DataFrame({"gene_id": ["g1"], "module_id": ["M000"]}),
                        "gtex", "cortex", "isograph", 0, "A")
    part = parts / "isograph__gtex__cortex__seed0__A.parquet"
    edge = st._edges_path("gtex", "cortex", "isograph", 0, "A")
    assert part.exists() and not edge.exists()

    def resume_done(save_edges):
        return part.exists() and (not save_edges or edge.exists())

    assert resume_done(save_edges=False), "a plain resume should still skip"
    assert not resume_done(save_edges=True), "an edge-saving pass must not skip"

    st._write_edges(_planted_edges().head(3), "gtex", "cortex", "isograph", 0, "A")
    assert resume_done(save_edges=True), "once both exist, resume should skip"


# --------------------------------------------------------------------------- #
# sweep clustering = production's edge-weighted detector
# --------------------------------------------------------------------------- #
def test_sweep_clustering_uses_edge_weights():
    """Two groups joined by every cross pair, but only weakly: the topology is one complete
    graph, so only weighted Leiden can recover the groups. Unweighted Leiden (the pre-
    2026-09-17 sweep) sees K40 and cannot split it along the weights."""
    from itertools import combinations

    from isograph_benchmark.real_data.sweep_leiden import _build_module_table

    genes = [f"g{i:02d}" for i in range(40)]
    group = {g: i // 20 for i, g in enumerate(genes)}
    edges = pd.DataFrame(
        [(a, b, 1.0 if group[a] == group[b] else 0.01) for a, b in combinations(genes, 2)],
        columns=["source", "target", "weight"])
    mods = _build_module_table(edges, genes, 1.0, seed=13, min_module_size=20)
    got = mods.groupby("module_id")["gene_id"].apply(frozenset)
    assert set(got) == {frozenset(genes[:20]), frozenset(genes[20:])}
