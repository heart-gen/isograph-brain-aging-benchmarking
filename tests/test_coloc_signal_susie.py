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
    P12_PRIMARY,
    P12_SWEEP,
    apply_hierarchy,
    best_signal_pair,
    results_dir,
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
