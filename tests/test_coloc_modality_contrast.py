"""The per-gene sQTL-vs-eQTL contrast is only meaningful if the pairing is real.

Every test here pins one of the properties the claim rests on: that a cell enters the
paired analysis only when BOTH modalities were testable, that the conditional posterior
divides out the QTL power difference, that the McNemar test is computed on the
discordant genes and no others, and that the GWAS N / case fraction -- which the module
docstring says does not affect the posteriors -- genuinely does not.
"""
from __future__ import annotations

import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data.coloc_modality_contrast import (
    ANALYSES,
    MIN_SHARED_SNPS,
    P12_PRIMARY,
    PP4_CALL,
    _FIXED_CASE_N,
    collapse_genes,
    pair_cells,
    paired_tests,
)

R_BIN = "/ocean/projects/bio250020p/shared/opt/env/R_env/bin/Rscript"


def _abf(rows) -> pd.DataFrame:
    """Long coloc.abf output: one row per (gene, locus, tissue, modality, p12)."""
    out = []
    for gene, tissue, modality, pp3, pp4, nsnps in rows:
        out.append({
            "analysis": "aging__scz", "trait": "scz", "LOCUS_ID": "locus01_chr1",
            "gene": gene, "symbol": gene, "module_id": "M001", "go_invisible": True,
            "tissue": tissue, "modality": modality, "p12": P12_PRIMARY,
            "PP3": pp3, "PP4": pp4, "nsnps": nsnps,
        })
    return pd.DataFrame(out)


def test_cell_needs_both_modalities_to_be_paired():
    # G1 has both arms; G2 has only sQTL and must not enter the paired table at all.
    abf = _abf([("G1", "Brain_Cortex", "sQTL", 0.1, 0.8, 500),
                ("G1", "Brain_Cortex", "eQTL", 0.5, 0.2, 500),
                ("G2", "Brain_Cortex", "sQTL", 0.1, 0.9, 500)])
    wide = pair_cells(abf)
    assert set(wide["gene"]) == {"G1"}


def test_shared_snp_floor_drops_the_cell_not_just_the_arm():
    # The eQTL arm is below the floor, so the gene is unpaired and the sQTL arm must
    # NOT be carried through alone -- that would credit splicing with a comparison
    # expression was never given the chance to make.
    abf = _abf([("G1", "Brain_Cortex", "sQTL", 0.1, 0.9, 500),
                ("G1", "Brain_Cortex", "eQTL", 0.4, 0.3, MIN_SHARED_SNPS - 1)])
    assert pair_cells(abf).empty


def test_conditional_posterior_divides_out_the_power_difference():
    # Same conditional sharing (PP4/(PP3+PP4) = 0.8) but the eQTL arm has far more
    # total signal. Unconditional PP4 favours eQTL; the conditional statistic ties.
    abf = _abf([("G1", "Brain_Cortex", "sQTL", 0.05, 0.20, 500),
                ("G1", "Brain_Cortex", "eQTL", 0.15, 0.60, 500)])
    w = pair_cells(abf).iloc[0]
    assert w["PP4_eQTL"] > w["PP4_sQTL"]
    assert w["cond_sQTL"] == pytest.approx(0.8)
    assert w["cond_eQTL"] == pytest.approx(0.8)


def test_mcnemar_uses_only_discordant_genes():
    genes = pd.DataFrame({
        "PP4_sQTL": [0.9, 0.9, 0.1, 0.1, 0.9, 0.9, 0.9],
        "PP4_eQTL": [0.1, 0.1, 0.9, 0.1, 0.1, 0.9, 0.1],
        "cond_sQTL": np.nan, "cond_eQTL": np.nan,
    })
    r = paired_tests(genes, call=PP4_CALL)
    assert (r["splicing_only"], r["expression_only"]) == (4, 1)
    assert r["both"] == 1 and r["neither"] == 1
    # Concordant genes carry no information and must not move the P value.
    assert r["mcnemar_p"] == pytest.approx(
        __import__("scipy.stats", fromlist=["stats"]).binomtest(4, 5, 0.5).pvalue)
    assert r["splicing_share_discordant"] == pytest.approx(4 / 5)


def test_tissue_collapse_is_symmetric_across_modalities():
    # Best tissue is taken per modality over the SAME set of paired cells.
    abf = _abf([("G1", "Brain_Cortex", "sQTL", 0.2, 0.7, 500),
                ("G1", "Brain_Cortex", "eQTL", 0.2, 0.3, 500),
                ("G1", "Brain_Amygdala", "sQTL", 0.2, 0.4, 500),
                ("G1", "Brain_Amygdala", "eQTL", 0.2, 0.9, 500)])
    g = collapse_genes(pair_cells(abf))
    assert len(g) == 1
    row = g.iloc[0]
    assert row["n_tissue"] == 2
    assert row["PP4_sQTL"] == pytest.approx(0.7)
    assert row["PP4_eQTL"] == pytest.approx(0.9)


def test_fixed_case_counts_reconcile_with_the_trait_registry():
    """The two traits with no per-SNP N must sum back to the registry's `fixed_n`."""
    from isograph_benchmark.real_data import gwas_traits as gt
    for trait, n_cas in _FIXED_CASE_N.items():
        spec = gt.TRAITS[trait]
        assert spec.fixed_n is not None
        assert 0 < n_cas < spec.fixed_n


def test_registered_analyses_match_the_coloc_layer():
    """Every analysis contrasted here must be one the CLPP layer already reports."""
    from isograph_benchmark.real_data.coloc_modality_contrast import coloc_dir
    for analysis in ANALYSES:
        assert (coloc_dir(analysis) / "susie" / "loci_testable.tsv").exists(), analysis


@pytest.mark.skipif(not shutil.which(R_BIN) and not __import__("os").path.exists(R_BIN),
                    reason="R_env Rscript not available")
def test_coloc_abf_posteriors_do_not_depend_on_gwas_n_or_s():
    """Pins the docstring claim that N and s are provenance-only.

    Both the module docstring and 20.coloc_modality_abf.R assert that supplying
    beta+varbeta for a case-control dataset makes coloc.abf's posteriors independent of
    N and s. If a future coloc release changes that, the LBD/ALS case fractions -- which
    come from publications rather than the sumstats -- would start to matter silently.
    """
    script = """
    suppressMessages(library(coloc)); set.seed(1); n <- 300
    d2 <- list(beta=rnorm(n,0,0.1), varbeta=rep(0.01,n), snp=paste0("s",1:n),
               type="quant", sdY=1)
    mk <- function(s,N) list(beta=rnorm(n,0,0.1), varbeta=rep(0.01,n),
                             snp=paste0("s",1:n), type="cc", s=s, N=N)
    set.seed(2); a <- suppressMessages(coloc.abf(mk(0.30,50000), d2))$summary
    set.seed(2); b <- suppressMessages(coloc.abf(mk(0.05, 1000), d2))$summary
    cat(isTRUE(all.equal(as.numeric(a), as.numeric(b))))
    """
    p = subprocess.run([R_BIN, "-e", script], capture_output=True, text=True, timeout=300)
    assert p.returncode == 0, p.stderr[-2000:]
    assert p.stdout.strip().endswith("TRUE"), p.stdout[-2000:]
