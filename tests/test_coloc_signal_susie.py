"""Signal-level coloc has to earn its place over coloc.abf, and stay honest about LD.

Three properties carry the layer:

  * Where there is genuinely ONE causal variant, `coloc.susie` and `coloc.abf` must
    agree. If they disagreed there, the upgrade would be introducing an artefact rather
    than removing an assumption, and every re-anchored locus would be suspect.
  * Where there are TWO signals and only one is shared, they must diverge -- that is the
    entire reason for the upgrade, and it is what `coloc.abf`'s single-causal-variant
    assumption cannot represent.
  * The estimator hierarchy must be complete: a cell SuSiE cannot speak to falls back to
    abf and is flagged, never silently dropped. Otherwise the denominator shrinks toward
    the loci that happened to fine-map and the yield looks better than it is.

The first two are checked in R against the real `coloc` package rather than mocked,
because the claim is about that package's behaviour, not about our wrapper's.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data.coloc_signal_susie import (
    FALLBACK_REASONS,
    MAX_SNPS_PRIMARY,
    P12_PRIMARY,
    P12_SWEEP,
    PRIOR_ROBUSTNESS,
    apply_hierarchy,
    best_signal_pair,
    fallback_reasons,
    p12_survival,
    prior_robustness,
    results_dir,
    scope_fallback,
    signal_root,
)

R_BIN = "/ocean/projects/bio250020p/shared/opt/env/R_env/bin/Rscript"
_KEYS = ["analysis", "trait", "LOCUS_ID", "gene", "tissue", "modality"]


def _pairs(rows) -> pd.DataFrame:
    base = dict(analysis="aging__lbd", trait="lbd", LOCUS_ID="locus08_chr4",
                gene="ENSG1", tissue="Brain_Cortex", modality="sQTL",
                p12=P12_PRIMARY, PP3=0.0, PP4=0.0, cs_matches_gtex=True,
                idx1=1, idx2=1, nsnps=1000)
    return pd.DataFrame([{**base, **r} for r in rows])


# --------------------------------------------------------------------------- #
# Cell collapse
# --------------------------------------------------------------------------- #
def test_best_signal_pair_keeps_the_maximum_but_records_how_many_were_tested():
    """`coloc.susie` returns a posterior per (GWAS cs x QTL cs) pair. Keeping the best is
    a maximum over several tests, so the count has to travel with it."""
    d = _pairs([{"idx2": 1, "PP4": 0.12}, {"idx2": 2, "PP4": 0.97},
                {"idx2": 3, "PP4": 0.40}])
    out = best_signal_pair(d)
    assert len(out) == 1
    assert out["PP4"].iloc[0] == pytest.approx(0.97)
    assert out["n_signal_pairs"].iloc[0] == 3


def test_n_signal_pairs_counts_pairs_not_rows_across_the_p12_sweep():
    """The R stage emits every signal pair once per prior. A raw row count would report
    3 pairs as 15; counting (idx1, idx2) is exact even when a prior drops a pair."""
    rows = [{"idx2": i, "p12": p, "PP4": 0.5}
            for p in P12_SWEEP for i in (1, 2, 3)]
    # coloc.susie failed for one pair at one prior -- the sweep is not always complete.
    rows = [r for r in rows if not (r["idx2"] == 3 and r["p12"] == P12_SWEEP[0])]
    out = best_signal_pair(_pairs(rows))
    assert out["n_signal_pairs"].iloc[0] == 3


def test_all_introns_pairs_are_identified_by_intron_not_by_credible_set_index():
    """SuSiE numbers credible sets from 1 within each fit, so every intron of a gene has
    an (idx1=1, idx2=1). Identifying a pair by the indices alone collapses a 12-intron
    gene onto the pair count of one intron, and hides that its PP4 is a maximum over 12
    tests -- which is the whole thing the all-introns arm has to be honest about."""
    rows = [{"phenotype_id": f"chr19:{s}:{s + 900}:clu_1_-:ENSG1", "idx2": 1, "PP4": pp}
            for s, pp in ((17630750, 0.31), (17641556, 0.88), (17650000, 0.12))]
    out = best_signal_pair(_pairs(rows))
    assert len(out) == 1
    assert out["n_signal_pairs"].iloc[0] == 3
    assert out["n_phenotypes_tested"].iloc[0] == 3
    # The winning intron must be named: the audit asks whether it IS the curated event.
    assert out["phenotype_id"].iloc[0] == "chr19:17641556:17642456:clu_1_-:ENSG1"


def test_representative_arm_reports_one_phenotype_so_the_arms_stay_comparable():
    """Adding phenotype_id to the pair identity must not move the representative arm: it
    is constant per cell there, so the count is unchanged and n_phenotypes_tested is 1."""
    rows = [{"phenotype_id": "ENSG1.12", "idx2": i, "PP4": pp}
            for i, pp in ((1, 0.12), (2, 0.97), (3, 0.40))]
    out = best_signal_pair(_pairs(rows))
    assert out["n_signal_pairs"].iloc[0] == 3
    assert out["n_phenotypes_tested"].iloc[0] == 1


def test_all_introns_results_land_beside_representative_not_on_top_of_it():
    """Both sqtl modes name their shards <analysis>__<tissue>.parquet, so they must not
    share a directory: an all-introns run would otherwise overwrite the representative
    results it exists to be compared against."""
    base = Path("/tmp/coloc_signal")
    rep = results_dir(base, "representative")
    allx = results_dir(base, "all")
    assert rep == base
    assert allx != rep
    assert rep in allx.parents


def test_a_raised_snp_guard_resolves_to_shards_the_primary_meta_never_reads():
    """The 2026-09-10 aging__ad recovery at MAX_SNPS=30000 was written into the primary
    directories and replaced the 12,000 AD shards in place. A non-default guard must
    resolve to shard dirs disjoint from the primary ones, for both sqtl modes and any arm."""
    prim = signal_root("switch")
    sens = signal_root("switch", 30000)
    assert signal_root("switch", MAX_SNPS_PRIMARY) == prim
    assert sens != prim
    shard_dirs = {results_dir(r, sq) / "susie"
                  for r in (prim, sens) for sq in ("representative", "all")}
    assert len(shard_dirs) == 4
    # `_load_susie` lists one directory, but none may nest inside another either.
    for a in shard_dirs:
        for b in shard_dirs:
            assert a == b or a not in b.parents
    assert signal_root("background", 30000) == sens / "arms" / "background"


def test_the_primary_keeps_every_fallback_but_a_sensitivity_root_only_what_it_reran():
    """Primary: an analysis whose GWAS fine-mapped nowhere still stands on abf, and must.
    Sensitivity: only the re-run analyses, or the other five would be counted as scored in
    an arm they were never run in."""
    abf = pd.DataFrame({"analysis": ["aging__ad", "aging__als", "aging__scz"],
                        "PP4": [0.8, 0.1, 0.2]})
    assert len(scope_fallback(abf, {"aging__ad"})) == 3
    assert scope_fallback(abf, {"aging__ad"}, 30000)["analysis"].tolist() == ["aging__ad"]


def test_gtex_agreement_filter_excludes_signals_gtex_did_not_find():
    """The re-fit uses a REFERENCE LD panel, so a credible set GTEx's own in-sample
    fine-mapping never found is the signature of an LD artefact. The primary arm drops
    it; the sensitivity arm keeps it. Both must be reachable."""
    d = _pairs([{"idx2": 1, "PP4": 0.99, "cs_matches_gtex": False},
                {"idx2": 2, "PP4": 0.30, "cs_matches_gtex": True}])
    primary = best_signal_pair(d, require_gtex_match=True)
    assert primary["PP4"].iloc[0] == pytest.approx(0.30), \
        "the unmatched 0.99 signal must not become the headline"
    loose = best_signal_pair(d, require_gtex_match=False)
    assert loose["PP4"].iloc[0] == pytest.approx(0.99)


def test_no_gtex_matched_signal_yields_no_cell_rather_than_a_false_one():
    d = _pairs([{"idx2": 1, "PP4": 0.99, "cs_matches_gtex": False}])
    assert best_signal_pair(d, require_gtex_match=True).empty


# --------------------------------------------------------------------------- #
# Estimator hierarchy
# --------------------------------------------------------------------------- #
def _abf_cells(rows) -> pd.DataFrame:
    base = dict(analysis="aging__lbd", trait="lbd", LOCUS_ID="locus08_chr4",
                gene="ENSG1", tissue="Brain_Cortex", modality="sQTL",
                PP3=0.5, PP4=0.5, nsnps=1000)
    return pd.DataFrame([{**base, **r} for r in rows])


def test_hierarchy_prefers_susie_and_falls_back_to_abf_without_dropping_cells():
    susie = best_signal_pair(_pairs([{"PP4": 0.97}]))
    abf = _abf_cells([
        {"gene": "ENSG1", "PP4": 0.80},                       # superseded by susie
        {"gene": "ENSG2", "PP4": 0.60},                       # susie could not speak
        {"gene": "ENSG1", "modality": "eQTL", "PP4": 0.10},   # other modality
    ])
    out = apply_hierarchy(susie, abf)
    assert len(out) == 3, "no cell may be dropped by the hierarchy"
    assert set(out["estimator"]) == {"susie", "abf"}
    g1s = out[(out.gene == "ENSG1") & (out.modality == "sQTL")]
    assert len(g1s) == 1 and g1s["estimator"].iloc[0] == "susie"
    assert g1s["PP4"].iloc[0] == pytest.approx(0.97), "abf must not override susie"
    assert out[out.gene == "ENSG2"]["estimator"].iloc[0] == "abf"


def test_hierarchy_is_a_pure_addition_when_susie_is_empty():
    """Every cell falling back is a legitimate outcome (a GWAS that does not fine-map),
    and must leave the abf grid exactly intact."""
    abf = _abf_cells([{"gene": f"ENSG{i}"} for i in range(4)])
    out = apply_hierarchy(pd.DataFrame(columns=[*_KEYS, "PP3", "PP4"]), abf)
    assert len(out) == len(abf)
    assert set(out["estimator"]) == {"abf"}


def test_p12_sweep_brackets_the_coloc_default():
    """The default p12 = 1e-5 must sit strictly inside the sweep, so a hit can be
    reported as surviving a RANGE of priors rather than one arbitrary choice."""
    assert P12_PRIMARY in P12_SWEEP
    assert min(P12_SWEEP) < P12_PRIMARY < max(P12_SWEEP)


# --------------------------------------------------------------------------- #
# Evidence-strength descriptors
# --------------------------------------------------------------------------- #
def test_p12_survival_reports_the_lowest_prior_that_still_calls():
    """PICALM/AD in cortex, as measured: the call holds from the primary prior up, and not
    below it. That is a different claim from a call that holds at 1e-6."""
    rows = [{"p12": p, "PP4": pp}
            for p, pp in zip(P12_SWEEP, (0.303, 0.685, 0.813, 0.956, 0.978))]
    s = p12_survival(_pairs(rows), _KEYS)
    assert s["p12_min_call"].iloc[0] == pytest.approx(P12_PRIMARY)
    assert s["PP4_at_p12_sweep_min"].iloc[0] == pytest.approx(0.303)
    assert prior_robustness(s["p12_min_call"].iloc[0]) == "primary_prior"


def test_prior_robustness_labels_every_point_of_the_sweep():
    assert [prior_robustness(p) for p in P12_SWEEP] == [
        "robust", "intermediate", "primary_prior", "permissive_prior", "permissive_prior"]
    assert prior_robustness(np.nan) == "none"
    assert set(PRIOR_ROBUSTNESS) == {"robust", "intermediate", "primary_prior",
                                     "permissive_prior", "none"}


def test_prior_survival_is_read_on_gtex_matched_pairs_only():
    """An unmatched pair calling at 1e-6 is exactly what the agreement filter removes. It
    must not make the cell look robust to the prior through the back door."""
    rows = []
    for p in P12_SWEEP:
        rows.append({"idx2": 1, "p12": p, "PP4": 0.99, "cs_matches_gtex": False})
        rows.append({"idx2": 2, "p12": p, "PP4": 0.85 if p >= P12_PRIMARY else 0.50,
                     "cs_matches_gtex": True})
    out = best_signal_pair(_pairs(rows))
    assert out["p12_min_call"].iloc[0] == pytest.approx(P12_PRIMARY)
    assert out["PP4_at_p12_sweep_min"].iloc[0] == pytest.approx(0.50)


def test_a_pair_below_the_shared_snp_floor_is_not_a_signal_level_result():
    d = _pairs([{"idx2": 1, "PP4": 0.95, "n_shared": 60},
                {"idx2": 2, "PP4": 0.40, "n_shared": 4000}])
    assert best_signal_pair(d)["PP4"].iloc[0] == pytest.approx(0.40)


def test_fallback_reason_names_the_first_place_a_cell_left_the_pipeline():
    base = dict(analysis="aging__ad", trait="ad", gene="ENSG1", tissue="Brain_Cortex",
                modality="sQTL", estimator="abf", PP4=0.85)
    cells = pd.DataFrame([
        {**base, "LOCUS_ID": "big"},
        {**base, "LOCUS_ID": "thin"},
        {**base, "LOCUS_ID": "flat"},
        {**base, "LOCUS_ID": "absent"},
        {**base, "LOCUS_ID": "ok", "gene": "Q0"},
        {**base, "LOCUS_ID": "ok", "gene": "Q1"},
        {**base, "LOCUS_ID": "ok", "gene": "Q2"},
        {**base, "LOCUS_ID": "ok", "gene": "Q3", "estimator": "susie"},
    ])
    st = dict(analysis="aging__ad")
    status = pd.DataFrame([
        {**st, "LOCUS_ID": "big", "fitted": False, "n_cs": np.nan,
         "reason": "15713 SNPs > MAX_SNPS (long-range LD)"},
        {**st, "LOCUS_ID": "thin", "fitted": False, "n_cs": np.nan, "reason": "<20 usable SNPs"},
        {**st, "LOCUS_ID": "flat", "fitted": True, "n_cs": 0, "reason": "no GWAS credible set"},
        {**st, "LOCUS_ID": "ok", "fitted": True, "n_cs": 1, "reason": "ok"},
    ])
    pk = dict(analysis="aging__ad", trait="ad", LOCUS_ID="ok", tissue="Brain_Cortex",
              modality="sQTL")
    pairs = pd.DataFrame([{**pk, "gene": "Q1", "cs_matches_gtex": False},
                          {**pk, "gene": "Q2", "cs_matches_gtex": True},
                          {**pk, "gene": "Q3", "cs_matches_gtex": True}])
    why = fallback_reasons(cells, status, pairs).tolist()
    assert why == ["gwas_locus_over_max_snps", "gwas_too_few_snps", "gwas_no_credible_set",
                   "gwas_not_in_stage_a", "no_qtl_credible_set", "qtl_cs_not_matching_gtex",
                   "susie_pair_filtered", None]
    assert set(filter(None, why)) == set(FALLBACK_REASONS)
    # Without the agreement filter an unmatched pair is not a reason to fall back.
    loose = fallback_reasons(cells, status, pairs, require_gtex_match=False).tolist()
    assert loose[5] == "susie_pair_filtered"


# --------------------------------------------------------------------------- #
# Against the real coloc package
# --------------------------------------------------------------------------- #
_R_AGREE = r'''
suppressPackageStartupMessages({library(coloc); library(susieR)})
set.seed(13)
p <- 400; n <- 5000
## An LD block with local correlation, so fine-mapping has something to resolve.
pos <- 1:p
R <- 0.9 ^ abs(outer(pos, pos, "-"))
R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
## coloc::annotate_susie -> .susie_setld indexes the LD matrix BY NAME, so an unnamed
## matrix fails with "no 'dimnames' attribute for array". The production scripts get
## names for free (read_ld sets dimnames from the .vars file); the fixture must too.
dimnames(R) <- list(as.character(1:p), as.character(1:p))
L <- chol(R)

mk_z <- function(causal, beta) {
    b <- numeric(p); b[causal] <- beta
    as.vector(R %*% b) * sqrt(n) + as.vector(t(L) %*% rnorm(p)) * 0.35
}
fit <- function(z) {
    f <- susie_rss(z = z, R = R, n = n, L = 10, estimate_residual_variance = FALSE)
    coloc::annotate_susie(f, as.character(1:p), R)
}
d <- function(z) list(beta = z / sqrt(n), varbeta = rep(1 / n, p),
                      snp = as.character(1:p), type = "quant", sdY = 1, N = n)

## --- ONE shared causal variant: susie and abf must agree it colocalizes -----
z1 <- mk_z(200, 0.11); z2 <- mk_z(200, 0.11)
s <- coloc.susie(fit(z1), fit(z2), p12 = 1e-5)
a <- suppressWarnings(coloc.abf(d(z1), d(z2), p12 = 1e-5))
cat(sprintf("SINGLE\tsusie=%.4f\tabf=%.4f\n",
            max(s$summary$PP.H4.abf), a$summary[["PP.H4.abf"]]))

## --- TWO QTL signals, only one shared: this is where they must diverge -----
## GWAS has one causal variant at 200. QTL has 200 AND a stronger private one at 60.
zg <- mk_z(200, 0.11)
bq <- numeric(p); bq[200] <- 0.11; bq[60] <- 0.16
zq <- as.vector(R %*% bq) * sqrt(n) + as.vector(t(L) %*% rnorm(p)) * 0.35
s2 <- coloc.susie(fit(zg), fit(zq), p12 = 1e-5)
a2 <- suppressWarnings(coloc.abf(d(zg), d(zq), p12 = 1e-5))
cat(sprintf("TWO\tsusie=%.4f\tabf=%.4f\tn_cs_qtl=%d\n",
            max(s2$summary$PP.H4.abf), a2$summary[["PP.H4.abf"]],
            length(fit(zq)$sets$cs)))
'''


@pytest.mark.skipif(not shutil.which(R_BIN) and not __import__("os").path.exists(R_BIN),
                    reason="R_env not available")
def test_susie_and_abf_agree_on_a_single_causal_variant_and_diverge_on_two():
    r = subprocess.run([R_BIN, "-e", _R_AGREE], capture_output=True, text=True,
                       timeout=1800)
    assert r.returncode == 0, r.stderr[-3000:]
    got = {}
    for line in r.stdout.splitlines():
        if line.startswith(("SINGLE", "TWO")):
            f = line.split("\t")
            got[f[0]] = {k: float(v) for k, v in (x.split("=") for x in f[1:])}
    assert "SINGLE" in got and "TWO" in got, r.stdout

    # 1. One causal variant: the estimators must not disagree, or the upgrade is a bug.
    s, a = got["SINGLE"]["susie"], got["SINGLE"]["abf"]
    assert s > 0.8 and a > 0.8, f"both should colocalize: susie={s}, abf={a}"

    # 2. Two QTL signals with a stronger private one: SuSiE isolates the shared signal,
    #    while abf -- forced to assume a single causal variant per trait -- is dragged
    #    toward the private one and loses H4. This is the whole point of the upgrade.
    s2, a2, ncs = got["TWO"]["susie"], got["TWO"]["abf"], got["TWO"]["n_cs_qtl"]
    assert ncs >= 2, f"the two-signal fixture did not produce 2 credible sets: {ncs}"
    assert s2 > a2, (
        f"coloc.susie should recover the shared signal that coloc.abf dilutes: "
        f"susie={s2}, abf={a2}")
