"""Selection-symmetric switch-vs-abundance effect sizes (brainseq_switch_qtl --stage effect_size)."""
import numpy as np
import pandas as pd

from isograph_benchmark.real_data import brainseq_switch_qtl as bsq


def _paired(n: int = 40) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    return pd.DataFrame({
        "phenotype_id": [f"g{i}" for i in range(n)],
        "num_var_ab": rng.integers(100, 2000, n),
        "variant_id_sw": [f"s{i}" for i in range(n)],
        "variant_id_ab": [f"a{i}" for i in range(n)],
        "slope_sw": rng.normal(0, 0.3, n), "slope_se_sw": np.full(n, 0.1),
        "slope_ab": rng.normal(0, 0.3, n), "slope_se_ab": np.full(n, 0.1),
        "qval_sw": rng.uniform(0, 1, n), "qval_ab": rng.uniform(0, 1, n),
    })


def test_cross_lead_values_come_from_the_other_axis_lead_variant():
    m = _paired(3)
    sw_pairs = pd.DataFrame({"phenotype_id": m["phenotype_id"], "variant_id": m["variant_id_ab"],
                             "slope": [0.2, -0.4, 0.6], "slope_se": [0.1, 0.1, 0.2],
                             "pval_nominal": [0.1, 0.1, 0.1]})
    ab_pairs = pd.DataFrame({"phenotype_id": m["phenotype_id"], "variant_id": m["variant_id_sw"],
                             "slope": [-0.1, 0.3, 0.0], "slope_se": [0.1, 0.1, 0.1],
                             "pval_nominal": [0.1, 0.1, 0.1]})
    g = bsq.symmetric_effects(m, sw_pairs, ab_pairs).set_index("phenotype_id")
    assert np.allclose(g["z_sw_at_ablead"], [2.0, 4.0, 3.0])
    assert np.allclose(g["z_ab_at_swlead"], [1.0, 3.0, 0.0])
    assert np.allclose(g["z_sw_own"], (m["slope_sw"] / 0.1).abs().to_numpy())


def test_lead_variant_frequency_flags_rarer_switch_leads_and_their_inflated_slopes():
    rng = np.random.default_rng(3)
    n = 400
    maf_sw = rng.uniform(0.01, 0.10, n)
    maf_ab = rng.uniform(0.15, 0.45, n)
    inv_sd = 1 / np.sqrt(2 * maf_sw * (1 - maf_sw))
    m = pd.DataFrame({
        "af_sw": maf_sw, "af_ab": 1 - maf_ab,                 # allele orientation must not matter
        "slope_sw": 0.1 * inv_sd * rng.choice([-1, 1], n), "slope_se_sw": 0.05 * inv_sd,
        "slope_se_ab": 0.05 / np.sqrt(2 * maf_ab * (1 - maf_ab)),
    })
    r = bsq.lead_variant_frequency(m, "caudate", "ea_only")
    assert r["median_lead_maf_switch"] < r["median_lead_maf_abundance"]
    assert r["frac_lead_maf_below_005_switch"] > r["frac_lead_maf_below_005_abundance"]
    assert r["median_se_ratio_switch_over_abundance"] > 1
    assert r["corr_abs_slope_switch_vs_inv_genotype_sd"] > 0.99


def test_cli_dispatches_effect_size_stage_without_mapping(monkeypatch):
    calls = {}
    monkeypatch.setattr(bsq, "run_effect_size", lambda arm, regions: calls.update(arm=arm, regions=regions))
    monkeypatch.setattr(bsq, "run_map", lambda *a, **k: (_ for _ in ()).throw(AssertionError("mapped")))
    bsq.main(["--stage", "effect_size", "--arm", "ea_only", "--region", "caudate"])
    assert calls == {"arm": "ea_only", "regions": ["caudate"]}


def test_summary_reports_all_tested_and_both_significant_subsets_and_power_deciles():
    m = _paired(60)
    sw_pairs = pd.DataFrame({"phenotype_id": m["phenotype_id"], "variant_id": m["variant_id_ab"],
                             "slope": m["slope_sw"] * 0.5, "slope_se": 0.1, "pval_nominal": 0.1})
    ab_pairs = pd.DataFrame({"phenotype_id": m["phenotype_id"], "variant_id": m["variant_id_sw"],
                             "slope": m["slope_ab"] * 0.5, "slope_se": 0.1, "pval_nominal": 0.1})
    g = bsq.symmetric_effects(m, sw_pairs, ab_pairs)
    s = bsq.summarize_effects(g, "caudate", "ea_only", fdr=0.05)
    overall = s[(s["power_decile"] == "all") & (s["subset"] == "all_tested")]
    assert set(overall["contrast"]) == {c[0] for c in bsq.EFFECT_CONTRASTS}
    assert (overall["n"] == 60).all()
    assert (s["power_decile"] != "all").any()
