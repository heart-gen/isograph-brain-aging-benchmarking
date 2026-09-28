"""The funnel claims: the denominators that make each fraction mean something."""
import json

import pandas as pd
import pytest

mtt = pytest.importorskip("isograph_benchmark.real_data.module_trust_tables")


def _stab(method, trusted):
    return pd.DataFrame({
        "cohort": "brainseq", "region": "caudate", "method": method,
        "module_id": [f"M{i:03d}" for i in range(len(trusted))],
        "coassign_density": 0.4, "null_mean": 0.1, "best_match_jaccard": 0.2,
        "trusted": trusted,
    })


def _pairs(method, both_sig, concordant):
    return pd.DataFrame({
        "cohort": "brainseq", "region": "caudate", "method": method,
        "both_age_sig": both_sig, "sign_concordant": concordant,
        "driver_load_rho": [0.8] * len(both_sig),
    })


def _proj_summary():
    return pd.DataFrame({
        "method": ["isograph", "wgcna"],
        "direction": ["brainseq_to_gtex", "brainseq_to_gtex"],
        # n_modules is NOT the denominator: a module with no testable age effect is
        # dropped from the sign claim but still counted as a module.
        "n_modules": [40, 52], "n_age_testable": [38, 52], "sign_match": [33, 28],
    })


def test_sign_concordance_is_counted_over_both_significant_pairs_only():
    # Three pairs agree on direction but only two are detectable in both halves. The
    # claim is 2/2, not 3/4: a pair whose age effect is not significant in both halves
    # carries no direction to agree on, and putting it in the denominator would dilute
    # the fraction with pairs that were never tested.
    pairs = _pairs("isograph", [True, True, False, False], [True, True, True, False])
    claims = mtt.funnel_claims(_stab("isograph", [True]), pairs, _proj_summary())
    row = claims[claims["claim"] == "split_half_age_sign_concordance"].iloc[0]
    assert (int(row["num"]), int(row["den"])) == (2, 2)


def test_each_claim_carries_its_own_denominator():
    stab = pd.concat([_stab("isograph", [True, True, False]),
                      _stab("wgcna", [True, True])])
    pairs = pd.concat([_pairs("isograph", [True], [True]),
                       _pairs("wgcna", [True, True], [True, False])])
    claims = mtt.funnel_claims(stab, pairs, _proj_summary())
    got = {(r["claim"], r["method"]): (int(r["num"]), int(r["den"]))
           for _, r in claims.iterrows()}
    assert got[("modules_chance_trusted", "isograph")] == (2, 3)
    assert got[("modules_chance_trusted", "wgcna")] == (2, 2)
    assert got[("split_half_age_sign_concordance", "wgcna")] == (1, 2)
    # The transfer rows come from the projection summary's age-testable count, not from
    # its module count.
    assert got[("aging_axis_transfers__brainseq_to_gtex", "isograph")] == (33, 38)


def test_claims_keep_the_order_the_results_text_makes_them_in():
    claims = mtt.funnel_claims(_stab("isograph", [True]), _pairs("isograph", [True], [True]),
                               _proj_summary())
    assert list(dict.fromkeys(claims["claim"].astype(str))) == [
        "modules_chance_trusted", "split_half_age_sign_concordance",
        "aging_axis_transfers__brainseq_to_gtex"]


def test_a_measure_with_no_computable_pair_is_written_as_untested_not_as_zero(tmp_path):
    # An absent upstream input must survive into the table as a null row. Dropping it, or
    # writing 0.0, would turn "we could not test this" into "we tested it and found
    # nothing" -- the one substitution the stage's own verdict section forbids.
    stats = {
        "method": "isograph", "model": "linear", "n_pairs": 38, "n_concordant": 1,
        "n_permutations": 1000, "seed": 13, "median_gene_jaccard": 0.04,
        "measures": {
            "go_jaccard": {"n_finite": 24, "mean_matched": 0.097, "null_mean": 0.030,
                           "null_sd": 0.019, "p_emp": 0.001, "mean_concordant": 0.43,
                           "mean_discordant": 0.083},
            "structure_r": {"n_finite": 0, "mean_matched": None, "null_mean": None,
                            "null_sd": None, "p_emp": None, "mean_concordant": None,
                            "mean_discordant": None},
        },
    }
    f = tmp_path / "functional_preservation__isograph__linear__stats.json"
    f.write_text(json.dumps(stats))

    df = mtt.functional_preservation_rows(tmp_path).set_index("measure")
    assert df.loc["structure_r", "n_finite"] == 0
    assert pd.isna(df.loc["structure_r", "p_emp"])
    assert df.loc["go_jaccard", "mean_matched"] == pytest.approx(0.097)


def test_region_funnel_keeps_the_methods_apart():
    # The complementarity ledgers carry no method column of their own, so a region slice
    # that forgets to filter on method takes the median over both methods' modules.
    stab = pd.concat([_stab("isograph", [True, False]), _stab("wgcna", [True])])
    pairs = pd.concat([_pairs("isograph", [True], [True]), _pairs("wgcna", [False], [False])])
    comp = pd.DataFrame({
        "cohort": "brainseq", "region": "caudate",
        "method": ["isograph", "isograph", "wgcna", "wgcna"],
        "frac_dtu_without_dge": [0.2, 0.2, 0.9, 0.9],
        "frac_in_wgcna_age_modules": [0.4, 0.4, 0.9, 0.9],
        "age_sig": True,
        "drv_cds_changed": 1.0, "drv_utr_changed": 1.0,
        "drv_biotype_switch": 1.0, "drv_coding_status_change": 1.0,
    })
    out = mtt.region_funnel(stab, pairs, comp).set_index("method")
    assert out.loc["isograph", "median_frac_dtu_without_dge"] == pytest.approx(0.2)
    # WGCNA rows get no complementarity columns at all, rather than IsoGraph's values.
    assert pd.isna(out.loc["wgcna"].get("median_frac_dtu_without_dge", float("nan")))


def test_driver_structure_summarises_age_significant_isograph_modules_only():
    comp = pd.DataFrame({
        "method": ["isograph", "isograph", "isograph", "wgcna"],
        "age_sig": [True, True, False, True],
        "drv_cds_changed": [1.0, 0.5, 0.0, 0.0],
        "drv_utr_changed": [1.0, 1.0, 0.0, 0.0],
        "drv_biotype_switch": [0.5, 0.5, 0.0, 0.0],
        "drv_coding_status_change": [0.0, 0.0, 0.0, 0.0],
    })
    out = mtt.driver_structure(comp).set_index("switch_class")
    assert out.loc["CDS change", "n_modules"] == 2
    assert out.loc["CDS change", "mean_frac"] == pytest.approx(0.75)


def _proj_rows(method, pair, raw_both, raw_match, both, match):
    return pd.DataFrame({
        "method": method, "direction": "brainseq_to_gtex", "pair": pair,
        "raw_both_sig": raw_both, "raw_sign_match": raw_match,
        "both_sig": both, "sign_match": match,
        "age_r_target": 0.1, "age_z_target": 2.0,
    })


def test_raw_sign_agreement_is_reported_against_both_denominators():
    # Four modules: three agree on the raw sign but only two of those are significant in
    # both cohorts. The raw arm's whole point is that its denominator is small and its
    # agreement near-unanimous, so the table carries 2/2 AND 3/4 rather than picking one.
    proj = _proj_rows("isograph", "caudate",
                      raw_both=[True, True, False, False],
                      raw_match=[True, True, True, False],
                      both=[True, False, False, False],
                      match=[True, True, False, False])
    r = mtt.projection_sign_scale(proj).iloc[0]
    assert (int(r["raw_sign_match_both_sig"]), int(r["raw_n_both_sig"])) == (2, 2)
    assert int(r["raw_sign_match_all"]) == 3
    assert (int(r["std_sign_match_all"]), int(r["n_modules"])) == (2, 4)


def test_switch_axis_orientation_is_attached_to_isograph_rows_only():
    # WGCNA projects abundance features, which have no switch axis to flip. Its rows must
    # stay empty rather than inherit IsoGraph's orientation numbers for the same region.
    proj = pd.concat([
        _proj_rows("isograph", "caudate", [True], [True], [True], [True]),
        _proj_rows("wgcna", "caudate", [True], [True], [True], [True]),
    ])
    orient = pd.DataFrame({"pair": ["caudate"], "n_shared_genes": [100],
                           "n_orientable": [50], "n_flipped": [23],
                           "frac_flipped": [0.46], "median_abs_cosine": [0.50],
                           "median_abs_cosine_orientable": [0.80]})
    out = mtt.projection_sign_scale(proj, orient).set_index("method")
    assert out.loc["isograph", "n_orientable"] == 50
    assert pd.isna(out.loc["wgcna", "n_orientable"])


def test_resolution_sweep_keeps_the_wgcna_baseline_with_no_resolution():
    # The sweep is only interpretable beside the fixed baseline, and WGCNA has no Leiden
    # resolution at all -- so its rows carry a null resolution instead of being dropped or
    # given a placeholder number that would sort in among the IsoGraph arms.
    summary = pd.DataFrame({
        "method": ["isograph", "isograph_res0p5", "isograph_res12", "wgcna"],
        "cohort": "brainseq", "region": "caudate",
        "comparison": "within_cohort_split_half", "n_seeds": 5,
        "mean_ari": [0.45, 0.34, 0.46, 0.40], "sd_ari": 0.02,
        "mean_nmi": [0.39, 0.34, 0.50, 0.36], "sd_nmi": 0.01,
        "mean_n_common": 5911.4,
    })
    out = mtt.resolution_sensitivity(summary).set_index("method")
    assert out.loc["isograph", "resolution"] == mtt.PRODUCTION_RESOLUTION
    assert out.loc["isograph", "is_production"]
    assert out.loc["isograph_res0p5", "resolution"] == 0.5
    assert not out.loc["isograph_res0p5", "is_production"]
    assert out.loc["isograph_res12", "resolution"] == 12.0
    assert pd.isna(out.loc["wgcna", "resolution"])
    assert out.loc["wgcna", "method_family"] == "wgcna"
