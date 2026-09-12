"""BrainSEQ coloc: the arm gate, the donors the in-sample LD comes from, and the reuse of the
GTEx layer's cell logic on the S_g / A_g axes.

The reuse is the point to pin. The GTEx layer's hierarchy, gene collapse and paired tests are
called on relabelled columns rather than re-implemented, so a relabelling slip would silently
compare the wrong columns while every number still looked plausible.
"""
import pandas as pd
import pytest

from isograph_benchmark.real_data import coloc_brainseq as cb
from isograph_benchmark.real_data import coloc_modality_contrast as cmc
from isograph_benchmark.real_data.coloc_brainseq import (
    brainseq_nominations,
    build_hierarchy,
    check_arm,
    collapse_regions,
    contrast_table,
    gtex_crossref,
    pivot_axes,
    region_table,
    work_list,
)

K = dict(analysis="aging__lbd", trait="lbd", LOCUS_ID="locus08_chr4", symbol="SNCA",
         module_id="M1", go_invisible=True)
SWEEP = (1e-6, 5e-6, 1e-5, 5e-5, 1e-4)


def test_only_ea_only_with_both_checks_passed_is_accepted():
    with pytest.raises(SystemExit, match="ea_only"):
        check_arm("all_samples", failures=[])
    with pytest.raises(SystemExit, match="checks"):
        check_arm("ea_only", failures=["dlpfc: positive control FAIL (pi1 0.3 < 0.5)"])
    assert check_arm("ea_only", failures=[]) is None


def test_region_table_takes_the_donors_each_region_was_mapped_in(tmp_path, monkeypatch):
    monkeypatch.setattr(cb, "qtl_dir", lambda arm, region: tmp_path / region)
    q = tmp_path / "caudate" / "qtl"
    q.mkdir(parents=True)
    cov = pd.DataFrame({"PC1": [0.1, 0.2, 0.3]}, index=["Br2", "Br1", "Br3"])
    for ft in ("switch", "abundance"):
        cov.to_csv(q / f"covariates_used_{ft}.txt", sep="\t")
    pd.DataFrame({"modality": ["switch", "abundance"], "n_donors": [3, 3]}).to_csv(
        q / "map_summary.tsv", sep="\t", index=False)
    out = tmp_path / "out"
    out.mkdir()
    t = region_table("ea_only", regions=("caudate",), dest=out)
    assert t["n_donors"].tolist() == [3]
    assert t["gtex_tissue"].tolist() == ["Brain_Caudate_basal_ganglia"]
    assert (out / "donors_caudate.txt").read_text().split() == ["Br1", "Br2", "Br3"]

    cov.iloc[:2].to_csv(q / "covariates_used_abundance.txt", sep="\t")
    with pytest.raises(SystemExit, match="different donors"):
        region_table("ea_only", regions=("caudate",), dest=out)


def test_one_task_per_analysis_and_region_in_a_fixed_order():
    w = work_list(pd.DataFrame({"analysis": ["b", "a", "a"]}))
    assert w.values.tolist() == [["a", "caudate"], ["a", "dlpfc"], ["a", "hippocampus"],
                                 ["b", "caudate"], ["b", "dlpfc"], ["b", "hippocampus"]]


def test_hierarchy_scores_susie_where_both_fine_map_and_abf_elsewhere_with_no_gtex_filter():
    cell = {**K, "gene": "G1", "tissue": "caudate"}
    abf = pd.DataFrame(
        [{**cell, "modality": "S_g", "p12": p, "PP3": 0.05, "PP4": v, "nsnps": 3000}
         for p, v in zip(SWEEP, (0.5, 0.8, 0.9, 0.95, 0.97))]
        + [{**cell, "modality": "A_g", "p12": p, "PP3": 0.05, "PP4": 0.1, "nsnps": 3000}
           for p in SWEEP])
    # No cs_matches_gtex column at all: in-sample LD, so there is nothing to agree with.
    pairs = pd.DataFrame(
        [{**cell, "modality": "S_g", "phenotype_id": "G1.3", "p12": p, "idx1": 1, "idx2": 1,
          "n_shared": 2500, "PP3": 0.01, "PP4": v}
         for p, v in zip(SWEEP, (0.85, 0.9, 0.95, 0.97, 0.98))])
    status = pd.DataFrame([{"analysis": "aging__lbd", "LOCUS_ID": "locus08_chr4",
                            "fitted": True, "n_cs": 1, "reason": "ok"}])
    _, merged = build_hierarchy(pairs, abf, status)
    m = merged.set_index("modality")
    assert m.at["S_g", "estimator"] == "susie"
    assert m.at["S_g", "PP4"] == pytest.approx(0.95)
    assert m.at["S_g", "prior_robustness"] == "robust"
    assert m.at["A_g", "estimator"] == "abf"
    assert m.at["A_g", "fallback_reason"] == "no_qtl_credible_set"


def _merged_rows() -> pd.DataFrame:
    rows = []
    for g, (s, a) in {"G1": (0.95, 0.10), "G2": (0.05, 0.90), "G3": (0.02, 0.03)}.items():
        for region, bump in (("caudate", 0.0), ("dlpfc", -0.5)):
            rows += [
                {**K, "gene": g, "tissue": region, "modality": "S_g", "PP3": 0.02,
                 "PP4": max(s + bump, 0.01), "nsnps": 3000, "estimator": "susie"},
                {**K, "gene": g, "tissue": region, "modality": "A_g", "PP3": 0.02,
                 "PP4": max(a + bump, 0.01), "nsnps": 3000, "estimator": "abf"}]
    return pd.DataFrame(rows)


def test_axes_are_paired_then_collapsed_and_tested_by_the_gtex_layers_own_functions():
    lone = pd.DataFrame([{**K, "gene": "G4", "tissue": "caudate", "modality": "S_g",
                          "PP3": 0.1, "PP4": 0.99, "nsnps": 3000, "estimator": "abf"}])
    merged = pd.concat([_merged_rows(), lone], ignore_index=True)

    wide = pivot_axes(merged)
    assert "G4" not in set(wide["gene"])                      # no A_g partner: not paired
    assert "G4" in set(pivot_axes(merged, require_both=False)["gene"])
    assert {"PP4_S_g", "PP4_A_g", "cond_S_g", "cond_A_g", "estimator_S_g"} <= set(wide.columns)

    genes = collapse_regions(wide).set_index("gene")
    assert genes.at["G1", "PP4_S_g"] == pytest.approx(0.95)
    assert genes.at["G1", "n_region"] == 2
    assert genes.at["G1", "n_region_S_g_coloc"] == 1
    assert genes.at["G1", "max_region_S_g"] == "caudate"

    c = contrast_table(genes.reset_index())
    pooled = c[c["analysis"] == "POOLED"].iloc[0]
    assert pooled["switch_only"] == 1 and pooled["abundance_only"] == 1
    assert pooled["n_S_g_coloc"] == 1 and pooled["n_A_g_coloc"] == 1
    ref = cmc.paired_tests(genes.reset_index().rename(columns={
        "PP4_S_g": "PP4_sQTL", "PP4_A_g": "PP4_eQTL",
        "cond_S_g": "cond_sQTL", "cond_A_g": "cond_eQTL"}))
    assert pooled["mcnemar_p"] == pytest.approx(ref["mcnemar_p"])


def test_gtex_nominations_are_read_in_every_region_with_the_matched_tissue_call_flagged():
    keys = dict(analysis="aging__lbd", trait="lbd", LOCUS_ID="locus08_chr4", gene="G1")
    gnom = pd.DataFrame([{**keys, "symbol": "SNCA"}])
    gcells = pd.DataFrame([
        {**keys, "tissue": "Brain_Frontal_Cortex_BA9", "PP4_sQTL": 0.95},
        {**keys, "tissue": "Brain_Hippocampus", "PP4_sQTL": 0.30},
    ])
    x = gtex_crossref(gnom, gcells, pivot_axes(_merged_rows(), require_both=False))
    x = x.set_index("region")
    assert list(x.index) == ["caudate", "dlpfc", "hippocampus"]
    assert x.at["dlpfc", "gtex_call_in_matched_tissue"]
    assert not x.at["hippocampus", "gtex_call_in_matched_tissue"]
    assert not x.at["caudate", "gtex_call_in_matched_tissue"]    # no GTEx caudate cell
    assert x.at["caudate", "PP4_S_g"] == pytest.approx(0.95) and x.at["caudate", "S_g_coloc"]
    assert not x.at["hippocampus", "brainseq_tested"]


def test_brainseq_nominations_are_switch_axis_calls_flagged_against_gtex():
    genes = collapse_regions(pivot_axes(_merged_rows()))
    gnom = pd.DataFrame([{"analysis": K["analysis"], "trait": K["trait"],
                          "LOCUS_ID": K["LOCUS_ID"], "gene": "G1", "symbol": "SNCA"}])
    nom = brainseq_nominations(genes, gnom)
    assert nom["gene"].tolist() == ["G1"]
    assert nom["in_gtex_nominations"].tolist() == [True]
