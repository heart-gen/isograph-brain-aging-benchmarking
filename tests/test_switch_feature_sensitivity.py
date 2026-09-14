"""Cohort-specific published settings and the cross-region worst-case summary."""
import pandas as pd

from isograph_benchmark.real_data import switch_feature_sensitivity as sfs


def test_both_cohorts_share_the_production_transcript_filter():
    from isograph_benchmark.real_data import run_models, stability

    assert sfs.PUBLISHED_BY_COHORT["gtex"] == sfs.PUBLISHED_BY_COHORT["brainseq"]
    d = run_models._filter_expressed_transcripts.__defaults__
    assert (sfs.PUBLISHED_BY_COHORT["gtex"]["min_count"],
            sfs.PUBLISHED_BY_COHORT["gtex"]["min_fraction"]) == d
    assert all(spec["filter_transcripts"] for spec in stability.COHORTS.values())


def test_an_off_grid_published_setting_leads_the_expression_axis():
    unfiltered = dict(sfs.PUBLISHED, min_count=0.0, min_fraction=0.0)
    settings = sfs._settings("expression", unfiltered)
    assert sfs._label(settings[0], "expression") == "no filter"
    assert sum(all(s[k] == unfiltered[k] for k in unfiltered) for s in settings) == 1


def test_brainseq_expression_axis_is_unchanged():
    settings = sfs._settings("expression", sfs.PUBLISHED_BY_COHORT["brainseq"])
    assert [(s["min_count"], s["min_fraction"]) for s in settings] == list(sfs._EXPRESSION_GRID)


def test_region_worst_case_ignores_the_published_row_and_handles_no_significant_modules():
    summary = pd.DataFrame([
        {"cohort": "gtex", "region": "a", "setting": "pub", "is_published": True,
         "effect_pearson_vs_published": 1.0, "n_fdr_sig_published": 4,
         "n_published_sig_retained": 4, "n_sign_flips_among_published_sig": 0,
         "median_abs_feature_r_vs_published": 1.0},
        {"cohort": "gtex", "region": "a", "setting": "x", "is_published": False,
         "effect_pearson_vs_published": 0.8, "n_fdr_sig_published": 4,
         "n_published_sig_retained": 3, "n_sign_flips_among_published_sig": 1,
         "median_abs_feature_r_vs_published": 0.9},
        {"cohort": "gtex", "region": "b", "setting": "x", "is_published": False,
         "effect_pearson_vs_published": 0.95, "n_fdr_sig_published": 0,
         "n_published_sig_retained": 0, "n_sign_flips_among_published_sig": 0,
         "median_abs_feature_r_vs_published": 0.97},
    ])
    worst = sfs.region_worst_case(summary).set_index("region")
    assert worst.loc["a", "n_settings"] == 1
    assert worst.loc["a", "min_retained_frac"] == 0.75
    assert worst.loc["a", "max_sign_flips_among_published_sig"] == 1
    assert pd.isna(worst.loc["b", "min_retained_frac"])


def test_gtex_donor_ids_strip_the_tissue_sample_suffix():
    samples = ["GTEX-1117F-0011-R10a-SM-AHZ7F", "GTEX-1117F-3226-SM-5N9CT", "GTEX-111CU-0126"]
    assert sfs._donor_ids(pd.DataFrame(), samples, "gtex") == {"GTEX-1117F", "GTEX-111CU"}
