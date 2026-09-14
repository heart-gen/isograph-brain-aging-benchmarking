"""Cohort-specific published settings and the cross-region worst-case summary."""
import numpy as np
import pandas as pd

from isograph_benchmark.real_data import switch_feature_sensitivity as sfs


def test_both_cohorts_share_the_production_transcript_filter():
    from isograph_benchmark.real_data import run_models, stability

    assert sfs.PUBLISHED_BY_COHORT["gtex"] == sfs.PUBLISHED_BY_COHORT["brainseq"]
    pub = sfs.PUBLISHED_BY_COHORT["gtex"]
    assert pub["transcript_filter"] == "switching"
    d = run_models.filter_switching_transcripts.__kwdefaults__ or dict(zip(
        ("min_gene_count", "min_gene_fraction", "min_tx_count", "min_tx_prop", "min_tx_fraction"),
        run_models.filter_switching_transcripts.__defaults__))
    assert (pub["min_tx_prop"], pub["min_tx_fraction"]) == (d["min_tx_prop"], d["min_tx_fraction"])
    assert all(spec["filter_transcripts"] for spec in stability.COHORTS.values())


def test_published_setting_applies_exactly_the_production_filter():
    from isograph_benchmark.real_data import run_models

    rng = np.random.default_rng(0)
    counts = rng.poisson(rng.gamma(0.4, 40.0, size=(60, 1)), size=(60, 30)).astype(float)
    table = pd.DataFrame({"transcript_id": [f"t{i}" for i in range(60)],
                          "gene_id": [f"g{i // 3}" for i in range(60)]})
    _, via_setting = sfs.apply_transcript_filter(counts, table, sfs.PUBLISHED)
    _, via_production = run_models.filter_production_transcripts(counts, table)
    assert via_setting["transcript_id"].tolist() == via_production["transcript_id"].tolist()


def test_expression_axis_leads_with_the_published_filter_and_holds_it_once():
    settings = sfs._settings("expression", sfs.PUBLISHED)
    assert all(settings[0][k] == sfs.PUBLISHED[k] for k in sfs.PUBLISHED)
    assert sum(all(s[k] == sfs.PUBLISHED[k] for k in sfs.PUBLISHED) for s in settings) == 1
    labels = [sfs._label(s, "expression") for s in settings]
    assert "no filter" in labels and any(lab.startswith("legacy") for lab in labels)
    assert len(set(labels)) == len(labels)


def test_an_off_grid_published_setting_leads_the_expression_axis():
    unfiltered = dict(sfs.PUBLISHED, transcript_filter="none", min_tx_prop=0.0, min_tx_fraction=0.0)
    settings = sfs._settings("expression", unfiltered)
    assert sfs._label(settings[0], "expression") == "no filter"
    assert sum(all(s[k] == unfiltered[k] for k in unfiltered) for s in settings) == 1


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
