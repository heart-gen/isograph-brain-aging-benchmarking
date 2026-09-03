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

from isograph_benchmark.paths import stage_out
from isograph_benchmark.real_data.coloc_modality_contrast import (
    ANALYSES,
    ARMS_GENE_POOL,
    MIN_SHARED_SNPS,
    P12_PRIMARY,
    PP4_CALL,
    _FIXED_CASE_N,
    _WGCNA_ENRICH,
    _WGCNA_FIT,
    collapse_genes,
    out_dir,
    pair_cells,
    paired_tests,
    run_prep,
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


# --------------------------------------------------------------------------- #
# Arms: the method comparison lives in gene selection, so the arms must differ
# only there
# --------------------------------------------------------------------------- #
def test_arm_out_dirs_are_distinct_and_switch_keeps_the_top_level():
    """`switch` must keep its original directory or the completed run is orphaned."""
    base = stage_out("anchoring", "coloc_modality_contrast")
    assert out_dir("switch") == base
    seen = {arm: out_dir(arm) for arm in ARMS_GENE_POOL}
    assert len(set(seen.values())) == len(ARMS_GENE_POOL), seen
    for arm in ARMS_GENE_POOL:
        if arm != "switch":
            assert seen[arm].parent.name == "arms"
            assert seen[arm].name == arm


def test_unknown_arm_is_rejected():
    with pytest.raises(SystemExit):
        run_prep(analyses=(), arm="not_an_arm")


def test_wgcna_arm_registry_matches_module_enrichment_naming():
    """The WGCNA fit dir and its enrichment file are named differently.

    `wgcna_switch` lives in `wgcna_switch_only/` but is written as
    `wgcna_switch_modules.parquet`. Getting that mapping wrong would silently read
    another method's partition, which is exactly the class of bug
    partition_provenance.py exists to stop.
    """
    from isograph_benchmark.real_data.module_enrichment import METHOD_DIRS
    for arm, fit in _WGCNA_FIT.items():
        assert METHOD_DIRS[arm] == fit, (arm, fit, METHOD_DIRS.get(arm))
        assert _WGCNA_ENRICH[arm] == f"{arm}_modules.parquet"


def test_background_arm_is_a_superset_of_the_switch_arm():
    """The background arm must CONTAIN the switch genes at the same loci.

    This is what makes the switch-vs-non-switch split locus-matched. If a switch
    (locus, gene) pair were missing from the background pool, the two groups would sit
    at different loci and the comparison would be confounded by locus composition.
    """
    base = stage_out("anchoring", "coloc_modality_contrast")
    sw_f = base / "targets.parquet"
    bg_f = base / "arms" / "background" / "targets.parquet"
    if not (sw_f.exists() and bg_f.exists()):
        pytest.skip("arms not prepped in this checkout")
    s = pd.read_parquet(sw_f)
    b = pd.read_parquet(bg_f)
    ks = set(zip(s["analysis"], s["LOCUS_ID"], s["gene"]))
    kb = set(zip(b["analysis"], b["LOCUS_ID"], b["gene"]))
    assert ks <= kb, f"{len(ks - kb)} switch pairs missing from the background arm"
    assert int(b["is_switch"].sum()) == len(ks)


def test_every_arm_reuses_the_same_loci():
    """Loci are held fixed across arms; only the gene pool may vary.

    An arm may drop a locus that contains none of its genes, but it must never
    introduce a locus the switch arm did not have -- that would mean the GWAS side
    changed and the arms were no longer comparable.
    """
    base = stage_out("anchoring", "coloc_modality_contrast")
    sw_f = base / "targets.parquet"
    if not sw_f.exists():
        pytest.skip("switch arm not prepped in this checkout")
    s = pd.read_parquet(sw_f)
    ref = set(zip(s["analysis"], s["LOCUS_ID"]))
    for arm in ARMS_GENE_POOL:
        if arm == "switch":
            continue
        f = base / "arms" / arm / "targets.parquet"
        if not f.exists():
            continue
        t = pd.read_parquet(f)
        got = set(zip(t["analysis"], t["LOCUS_ID"]))
        assert got <= ref, f"{arm} introduced loci absent from the switch arm"


def test_genes_outside_a_module_are_not_silently_dropped():
    """The background arm is 87% genes with no module_id; they must survive pairing.

    `pivot_table` groups with dropna=True, so a NaN in any index key deletes the row.
    module_id / go_invisible are legitimately NaN for genes in no IsoGraph module --
    exactly the population the background arm exists to measure. Without the sentinel
    fill the background arm would collapse to the switch arm and still look valid,
    which is the kind of failure that produces a confident wrong answer.
    """
    rows = []
    for gene, mod in (("G1", "M001"), ("G2", None)):
        for modality, pp3, pp4 in (("sQTL", 0.1, 0.8), ("eQTL", 0.5, 0.2)):
            rows.append({
                "analysis": "aging__scz", "trait": "scz", "LOCUS_ID": "locus01_chr1",
                "gene": gene, "module_id": mod, "go_invisible": None,
                "tissue": "Brain_Cortex", "modality": modality,
                "p12": P12_PRIMARY, "PP3": pp3, "PP4": pp4, "nsnps": 500,
            })
    wide = pair_cells(pd.DataFrame(rows))
    assert set(wide["gene"]) == {"G1", "G2"}, "module-less gene was dropped"
    assert set(collapse_genes(wide)["gene"]) == {"G1", "G2"}


def test_sensitivity_loop_does_not_shadow_the_gene_pool_arm():
    """`run_meta`'s sensitivity loop must not rebind the `arm` parameter.

    `ARMS` (sensitivity: p12 / PP4 call / SNP floor) and `ARMS_GENE_POOL` (which genes
    are tested) are different axes that both got called "arm". A loop written
    `for arm, kw in ARMS:` silently overwrites the function argument, so every use of
    `arm` after the loop -- the background arm's switch/non-switch split, and the gene
    pool named in the report -- would take the value of the LAST sensitivity arm.
    Nothing raises; the split just never runs and the report is mislabelled.
    """
    import inspect
    from isograph_benchmark.real_data import coloc_modality_contrast as m
    src = inspect.getsource(m.run_meta)
    assert "for arm, kw in ARMS" not in src, "sensitivity loop shadows the arm parameter"
    assert "for sens_arm, kw in ARMS" in src
