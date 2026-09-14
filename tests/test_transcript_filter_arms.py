"""Switching filter semantics, arm layout (never inside a production store), fit comparison."""
import inspect

import numpy as np
import pandas as pd

from isograph_benchmark.paths import region_store
from isograph_benchmark.real_data import run_models
from isograph_benchmark.real_data import transcript_filter_arms as tfa


def _toy():
    n = 10
    rows = {
        # gene A: dominant t1; t2 used (share 1/3) in 2/10 samples only; t3 never reaches 10
        "t1": ("A", [100] * n), "t2": ("A", [50, 50] + [0] * (n - 2)), "t3": ("A", [5] * n),
        # gene B: never expressed enough for proportions to be estimable
        "t6": ("B", [1] * n), "t7": ("B", [2] * n),
        # gene C: t5 is abundant but always < 10% of its gene
        "t4": ("C", [1000] * n), "t5": ("C", [20] * n),
    }
    table = pd.DataFrame({"transcript_id": list(rows), "gene_id": [g for g, _ in rows.values()]})
    counts = np.array([c for _, c in rows.values()], dtype=float)
    return counts, table


def test_switching_filter_keeps_minority_used_isoforms_and_drops_trivial_shares():
    counts, table = _toy()
    c, t = tfa.switching_filter(counts, table)
    assert set(t["transcript_id"]) == {"t1", "t2", "t4"}
    assert c.shape == (3, 10)
    _, prod = run_models._filter_expressed_transcripts(counts, table)
    assert set(prod["transcript_id"]) == {"t1", "t4", "t5"}


def test_production_defaults_are_unchanged():
    region = inspect.signature(run_models.run_brainseq_region).parameters
    assert region["transcript_filter"].default is run_models._filter_expressed_transcripts
    assert region["out"].default is None
    sczd = inspect.signature(run_models.run_brainseq_caudate_sczd).parameters
    assert sczd["transcript_filter"].default is None and sczd["out"].default is None
    assert region["random_state"].default == 13 == sczd["random_state"].default == tfa.PRODUCTION_SEED


def test_every_target_has_a_floor_and_a_switching_arm_outside_production_stores():
    for region, analysis in tfa.TARGETS.items():
        arms = {a for r, a in tfa.FITS if r == region}
        assert {"switching", tfa.PRODUCTION_ARM[analysis]} <= arms
    assert ("caudate_sczd", "abundance") in tfa.FITS
    assert len(set(tfa.FITS)) == len(tfa.FITS)
    dirs = [tfa.arm_dir(r, a, s).resolve() for grid in tfa.GRIDS.values() for r, a, s in grid]
    assert len(set(dirs)) == len(dirs)  # no seed arm lands on another fit's directory
    for (region, _, _), d in zip([f for g in tfa.GRIDS.values() for f in g], dirs):
        assert region_store("brainseq", region).resolve() not in d.parents
        assert "isograph_vae" not in d.parts


def test_seed_grid_covers_every_arm_at_every_extra_seed():
    assert tfa.PRODUCTION_SEED not in tfa.EXTRA_SEEDS
    assert set(tfa.SEED_FITS) == {(r, a, s) for r, a in tfa.FITS for s in tfa.EXTRA_SEEDS}
    assert tfa.arm_dir("caudate", "switching", tfa.PRODUCTION_SEED).name == "switching"
    assert tfa.arm_dir("caudate", "switching", 14).name == "switching_seed14"


def _fit(region, arm, assign, effects, fdrs):
    modules = pd.DataFrame({"gene_id": list(assign), "module_id": list(assign.values())})
    trait = pd.DataFrame({"module_id": list(effects), "effect": list(effects.values()),
                          "fdr": [fdrs[k] for k in effects]})
    return {"region": region, "arm": arm, "modules": modules, "trait": trait, "kept": None,
            "meta": {}, "alpha": np.nan}


def test_compare_fits_is_label_invariant_and_requires_the_same_sign():
    genes = [f"g{i}" for i in range(40)]
    a = {g: ("M0" if i < 20 else "M1") for i, g in enumerate(genes)}
    b = {g: ("X9" if i < 20 else "X3") for i, g in enumerate(genes)}
    ref = _fit("caudate", "production", a, {"M0": 0.5, "M1": -0.2}, {"M0": 0.01, "M1": 0.5})
    same = _fit("caudate", "switching", b, {"X9": 0.4, "X3": 0.1}, {"X9": 0.001, "X3": 0.9})
    flip = _fit("caudate", "switching", b, {"X9": -0.4, "X3": 0.1}, {"X9": 0.001, "X3": 0.9})
    r = tfa.compare_fits(ref, same)
    assert r["ari"] == 1.0 and r["n_ref_sig_retained_in_qry"] == 1
    assert r["n_qry_sig_retained_in_ref"] == 1 and r["trait_gene_jaccard"] == 1.0
    assert tfa.compare_fits(ref, flip)["n_ref_sig_retained_in_qry"] == 0
    assert r["frac_ref_sig_genes_kept_same_sign"] == 1.0 and r["qry_sig_gene_rate"] == 0.5
    assert tfa.compare_fits(ref, flip)["frac_ref_sig_genes_kept_same_sign"] == 0.0
    s = tfa.fit_summary(ref)
    assert s["n_modules"] == 2 and s["n_trait_sig"] == 1 and s["n_genes_in_trait_sig"] == 20
