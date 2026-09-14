"""Switch-axis orientation and frozen-eigengene projection."""
import numpy as np
import pandas as pd

from isograph_benchmark.real_data import eigengene_projection as ep


def test_alignment_flips_an_inverted_axis_and_drops_unorientable_genes():
    src = pd.DataFrame({"gene": ["g1"] * 3 + ["g2"] * 3 + ["g3"] * 2,
                        "transcript": ["a", "b", "c", "d", "e", "f", "g", "h"],
                        "loading": [0.8, -0.5, -0.3, 0.7, -0.7, 0.0, 0.9, -0.4]})
    tgt = src.copy()
    tgt.loc[tgt["gene"] == "g1", "loading"] *= -1          # same axis, opposite pivot
    tgt.loc[tgt["gene"] == "g2", "loading"] = [0.0, 0.0, 1.0]  # orthogonal: a different axis
    tgt = tgt[tgt["transcript"] != "h"]                     # g3 keeps one shared transcript
    t = ep.align_switch_axes(src, tgt).set_index("gene")
    assert t.loc["g1", "usable"] and t.loc["g1", "sign"] == -1.0
    assert not t.loc["g2", "usable"]
    assert not t.loc["g3", "usable"]


def test_projecting_onto_identical_data_preserves_the_eigengene():
    rng = np.random.default_rng(0)
    factor = rng.normal(size=60)
    x = ep._zrows(np.outer(rng.uniform(0.5, 1.0, 20), factor) + rng.normal(scale=0.5, size=(20, 60)))
    w = ep.frozen_weights(x)
    e = ep.project(w, x)
    assert abs(np.corrcoef(e, factor)[0, 1]) > 0.9
    assert ep.signed_kme(w, x, e) > 0.5


def test_a_cohort_wide_age_factor_sets_the_raw_sign_but_not_the_null_centred_z():
    rng = np.random.default_rng(2)
    n_samples, n_features = 120, 400
    age = rng.uniform(20, 80, n_samples)
    cohort_factor = (age - age.mean()) / age.std() + rng.normal(scale=0.3, size=n_samples)
    # every feature loads on the same age-correlated factor; no module-specific aging
    x = ep._zrows(np.outer(rng.uniform(0.4, 0.8, n_features), cohort_factor)
                  + rng.normal(size=(n_features, n_samples)))
    module = x[:25]
    w = ep.frozen_weights(module)
    r_obs = ep._age_r(ep.project(w, module), age)[0]
    null = np.array([ep._age_r(ep.project(w, x[rng.choice(n_features, 25, replace=False)]), age)[0]
                     for _ in range(300)])
    z, p = ep.age_vs_null(r_obs, null)
    assert abs(r_obs) > 0.5          # the raw projection looks strongly age-associated
    assert abs(z) < 2.5 and p > 0.01  # but not beyond what any same-weight projection shows


def test_an_unaligned_sign_flip_destroys_preservation_and_alignment_restores_it():
    rng = np.random.default_rng(1)
    factor = rng.normal(size=80)
    x = ep._zrows(np.outer(np.ones(12), factor) + rng.normal(scale=0.3, size=(12, 80)))
    w = ep.frozen_weights(x)
    flip = np.array([-1.0] * 6 + [1.0] * 6)
    flipped = x * flip[:, None]
    unaligned = ep.signed_kme(w, flipped, ep.project(w, flipped))
    realigned = flipped * flip[:, None]
    aligned = ep.signed_kme(w, realigned, ep.project(w, realigned))
    assert aligned > 0.8 and abs(unaligned) < 0.2
