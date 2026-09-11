"""SMR/HEIDI sits beneath coloc, so its plumbing must not silently flip a sign or a position.

The failure modes that matter here raise no error. An ESD oriented to the wrong allele
reverses `b_SMR`. Positions from a different build than the LD panel mis-select the HEIDI
SNPs. A `.ma` with the columns in the wrong order parses cleanly. Each is pinned against the
formats SMR 1.4.2 actually reads (src/SMR_data.cpp, src/SMR_data_p2.cpp), not from memory.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data.smr_heidi import (
    AGREEMENT,
    ESD_COLS,
    FLIST_COLS,
    HEIDI_REJECT,
    MA_COLS,
    PEQTL_SMR,
    SMR_OUT_COLS,
    build_targets,
    check_qtl_source,
    classify,
    esd_rows,
    flist_rows,
    locus_chr,
    ma_from_loci,
    probe_coloc,
    probe_gene,
    read_smr,
    smr_command,
    smr_threshold,
    write_ma,
)

SNCA_E = "ENSG00000145335.17"
SNCA_I = "chr4:89835692:89836127:clu_43552_-:ENSG00000145335.17"


def _bim(rows):
    return pd.DataFrame(rows, columns=["bchr", "rsid", "cm", "bp", "A1", "A2"])


def _q(variants, af=None, slope=None):
    n = len(variants)
    return pd.DataFrame({"phenotype_id": [SNCA_I] * n, "variant_id": variants,
                         "af": af or [0.3] * n, "pval_nominal": [1e-9] * n,
                         "slope": slope or [0.5] * n, "slope_se": [0.05] * n})


# --------------------------------------------------------------------------- #
# ESD
# --------------------------------------------------------------------------- #
def test_esd_effect_allele_is_gtex_alt_whatever_the_panel_coding():
    """`af` and `slope` refer to GTEx's ALT allele. The panel may code the same SNP either way
    round; the ESD must keep ALT as A1 in both cases or b_SMR changes sign."""
    q = _q(["chr4_89835692_A_G_b38", "chr4_89836000_C_T_b38"], af=[0.3, 0.6],
           slope=[0.5, -0.1])
    br = pd.DataFrame({"variant_id": q["variant_id"], "rsid": ["rs1", "rs2"]})
    bim = _bim([(4, "rs1", 0, 90756843, "G", "A"), (4, "rs2", 0, 90757151, "C", "T")])
    e = esd_rows(q, br, bim)
    assert list(e.columns) == ["phenotype_id", *ESD_COLS]
    e = e.set_index("SNP")
    assert (e.at["rs1", "A1"], e.at["rs1", "A2"]) == ("G", "A")
    assert (e.at["rs2", "A1"], e.at["rs2", "A2"]) == ("T", "C")
    assert e.at["rs1", "Beta"] == pytest.approx(0.5)
    assert e.at["rs2", "Freq"] == pytest.approx(0.6)


def test_esd_position_is_the_hg19_panel_position_not_the_grch38_id():
    q = _q(["chr4_89835692_A_G_b38"])
    br = pd.DataFrame({"variant_id": q["variant_id"], "rsid": ["rs1"]})
    e = esd_rows(q, br, _bim([(4, "rs1", 0, 90756843, "A", "G")]))
    assert e["Bp"].iloc[0] == 90756843


def test_esd_drops_ambiguous_mismatched_unbridged_and_monomorphic_variants():
    q = _q(["chr1_10_A_T_b38",      # strand-ambiguous
            "chr1_20_A_C_b38",      # panel says G/T
            "chr1_30_A_G_b38",      # not in the bridge
            "chr1_40_A_G_b38",      # af = 1: SMR refuses it
            "chr1_50_C_A_b38"],     # kept
           af=[0.3, 0.3, 0.3, 1.0, 0.3])
    br = pd.DataFrame({"variant_id": ["chr1_10_A_T_b38", "chr1_20_A_C_b38",
                                      "chr1_40_A_G_b38", "chr1_50_C_A_b38"],
                       "rsid": ["rs1", "rs2", "rs4", "rs5"]})
    bim = _bim([(1, "rs1", 0, 10, "A", "T"), (1, "rs2", 0, 20, "G", "T"),
                (1, "rs4", 0, 40, "A", "G"), (1, "rs5", 0, 50, "A", "C")])
    assert esd_rows(q, br, bim)["SNP"].tolist() == ["rs5"]


# --------------------------------------------------------------------------- #
# GWAS .ma
# --------------------------------------------------------------------------- #
def test_ma_has_smr_columns_na_frequency_trait_n_and_the_coloc_exclusions(tmp_path):
    g = pd.DataFrame({"rsid": ["rs1", "rs2", "rs3", "rs3", "rs4"],
                      "pos": [100, 30_000_000, 200, 200, 300],
                      "a1": ["a", "C", "C", "C", "G"], "a2": ["g", "T", "T", "T", "A"],
                      "beta": [0.1, 0.2, 0.3, 0.3, 0.4], "se": [0.01, 0.02, 0.03, 0.03, 0.0],
                      "p": [1e-5, 0.3, 0.2, 0.2, 0.1], "z": [10, 10, 10, 10, 10]})
    excl = pd.DataFrame({"chr": [6], "start": [25_000_000], "stop": [34_000_000]})
    ma = ma_from_loci([(6, g)], excl, n=6618)
    assert list(ma.columns) == MA_COLS == ["SNP", "A1", "A2", "freq", "b", "se", "p", "n"]
    # rs2 is inside the MHC exclusion, rs3 is duplicated, rs4 has se = 0.
    assert ma["SNP"].tolist() == ["rs1", "rs3"]

    f = tmp_path / "t.ma"
    write_ma(ma, f)
    head, first = (line.split("\t") for line in f.read_text().splitlines()[:2])
    assert head == MA_COLS
    assert first[1:4] == ["A", "G", "NA"], "effect allele upper-cased; frequency written NA"
    assert first[7] == "6618"


def test_the_exclusion_only_applies_on_its_own_chromosome():
    g = pd.DataFrame({"rsid": ["rs1"], "pos": [30_000_000], "a1": ["A"], "a2": ["G"],
                      "beta": [0.1], "se": [0.01], "p": [1e-5], "z": [10]})
    excl = pd.DataFrame({"chr": [6], "start": [25_000_000], "stop": [34_000_000]})
    assert len(ma_from_loci([(4, g)], excl, n=10)) == 1


# --------------------------------------------------------------------------- #
# flist, probes, command
# --------------------------------------------------------------------------- #
def test_flist_is_the_seven_columns_read_probeinfolst_parses(tmp_path):
    esd = pd.DataFrame({"phenotype_id": [SNCA_I, SNCA_I, "ENSG00000999999.1"],
                        "Chr": [4, 4, 4], "SNP": ["rs1", "rs2", "rs3"],
                        "Bp": [100, 300, 500], "A1": ["A"] * 3, "A2": ["G"] * 3,
                        "Freq": [0.3] * 3, "Beta": [0.1] * 3, "se": [0.01] * 3, "p": [0.1] * 3})
    genes = pd.DataFrame({"gene": ["ENSG00000145335"], "symbol": ["SNCA"], "strand": ["-"],
                          "tss": [90759447]})
    fl = flist_rows(esd, genes, tmp_path)
    assert FLIST_COLS[0] == "Chr" and len(FLIST_COLS) == 7 and list(fl.columns) == FLIST_COLS
    s = fl.set_index("ProbeID")
    assert s.at[SNCA_I, "ProbeBp"] == 90759447 and s.at[SNCA_I, "Orientation"] == "-"
    assert s.at[SNCA_I, "Gene"] == "SNCA"
    assert ":" not in Path(s.at[SNCA_I, "PathOfEsd"]).name
    # No hg19 TSS: the median ESD position, orientation NA.
    assert s.at["ENSG00000999999.1", "ProbeBp"] == 500
    assert s.at["ENSG00000999999.1", "Orientation"] == "NA"


def test_probe_gene_reads_both_gtex_probe_forms():
    assert probe_gene(SNCA_E) == probe_gene(SNCA_I) == "ENSG00000145335"


def test_locus_chr_reads_the_locus_id_and_refuses_anything_else():
    assert locus_chr("locus60_chr11") == 11
    with pytest.raises(SystemExit):
        locus_chr("locus60")


def test_smr_command_passes_every_default_explicitly():
    cmd = smr_command(Path("ref"), Path("g.ma"), Path("b"), Path("p"), Path("o"))
    flags = dict(zip(cmd[1::2], cmd[2::2]))
    assert flags["--peqtl-smr"] == f"{PEQTL_SMR:g}"
    for f in ("--bfile", "--gwas-summary", "--beqtl-summary", "--extract-probe", "--out",
              "--peqtl-heidi", "--heidi-min-m", "--heidi-max-m", "--cis-wind", "--diff-freq",
              "--diff-freq-prop", "--thread-num"):
        assert f in flags


# --------------------------------------------------------------------------- #
# Output, thresholds, classification
# --------------------------------------------------------------------------- #
def test_read_smr_accepts_the_header_smr_writes(tmp_path):
    f = tmp_path / "x.smr"
    f.write_text("\t".join(SMR_OUT_COLS) + "\n"
                 + "\t".join(["P", "4", "SNCA", "1", "rs1", "4", "1", "A", "G", "0.3", "0.1",
                              "0.01", "1e-8", "0.5", "0.05", "1e-9", "0.2", "0.03", "1e-7",
                              "0.4", "12"]) + "\n")
    d = read_smr(f)
    assert d["p_HEIDI"].iloc[0] == pytest.approx(0.4)
    (tmp_path / "bad.smr").write_text("probeID\tb_SMR\n")
    with pytest.raises(SystemExit, match="lacks SMR output columns"):
        read_smr(tmp_path / "bad.smr")


def test_bonferroni_family_is_the_instrumented_probes():
    assert smr_threshold(20) == pytest.approx(0.0025)
    assert np.isnan(smr_threshold(0))


@pytest.mark.parametrize("pp4,p_smr,p_heidi,nsnp,expected", [
    (0.95, 1e-6, 0.30, 12, "coloc_and_smr_heidi_not_rejected"),
    (0.95, 1e-6, HEIDI_REJECT / 2, 12, "coloc_and_smr_heidi_rejected"),
    (0.95, 1e-6, np.nan, 2, "coloc_and_smr_heidi_untestable"),
    (0.95, 0.20, 0.30, 12, "coloc_smr_not_significant"),
    (0.95, np.nan, np.nan, np.nan, "coloc_no_instrument"),
    (0.10, 1e-6, 0.30, 12, "smr_without_coloc"),
    (np.nan, 1e-6, 0.30, 12, "smr_coloc_not_scored"),
    (np.nan, 0.20, 0.30, 12, "neither"),
])
def test_classification(pp4, p_smr, p_heidi, nsnp, expected):
    assert classify(pp4, p_smr, p_heidi, nsnp, threshold=1e-3) == expected
    assert expected in AGREEMENT


def test_a_heidi_rejection_never_turns_a_coloc_call_into_no_call():
    """The rejection is its own category beside the coloc call, never a downgrade of it."""
    v = classify(0.99, 1e-10, 1e-6, 20, threshold=1e-3)
    assert v.startswith("coloc_")


# --------------------------------------------------------------------------- #
# Targets, coloc join, source guard
# --------------------------------------------------------------------------- #
def test_targets_are_nominated_genes_in_tissues_where_the_sqtl_call_holds():
    k = dict(analysis="aging__lbd", trait="lbd", LOCUS_ID="locus08_chr4", gene="ENSG1",
             symbol="SNCA")
    nom = pd.DataFrame([k])
    cells = pd.DataFrame([
        {**k, "tissue": "Brain_Cortex", "PP4_sQTL": 0.97, "PP4_eQTL": 0.1,
         "estimator_sQTL": "susie", "phenotype_id": SNCA_I},
        {**k, "tissue": "Brain_Amygdala", "PP4_sQTL": 0.40, "PP4_eQTL": 0.1,
         "estimator_sQTL": "susie", "phenotype_id": SNCA_I},
        {**k, "gene": "ENSG2", "tissue": "Brain_Cortex", "PP4_sQTL": 0.99, "PP4_eQTL": 0.1,
         "estimator_sQTL": "abf", "phenotype_id": None},
    ])
    t = build_targets(nom, cells)
    assert t[["gene", "tissue"]].values.tolist() == [["ENSG1", "Brain_Cortex"]]
    assert t["chr"].iloc[0] == 4 and t["coloc_phenotype_id"].iloc[0] == SNCA_I


def test_probe_coloc_is_per_intron_and_abf_only_reaches_the_representative_intron():
    k = dict(analysis="aging__ad", trait="ad", LOCUS_ID="locus60_chr11", tissue="Brain_Cortex")
    pairs = pd.DataFrame([
        {**k, "gene": "G1", "modality": "sQTL", "phenotype_id": "i1", "p12": 1e-5,
         "cs_matches_gtex": True, "n_shared": 4000, "PP4": 0.81},
        {**k, "gene": "G1", "modality": "sQTL", "phenotype_id": "i1", "p12": 1e-6,
         "cs_matches_gtex": True, "n_shared": 4000, "PP4": 0.30},       # other prior
        {**k, "gene": "G1", "modality": "sQTL", "phenotype_id": "i2", "p12": 1e-5,
         "cs_matches_gtex": False, "n_shared": 4000, "PP4": 0.99},      # fails GTEx match
    ])
    hier = pd.DataFrame([
        {**k, "gene": "G2", "modality": "sQTL", "estimator": "abf", "PP4": 0.85},
        {**k, "gene": "G2", "modality": "eQTL", "estimator": "abf", "PP4": 0.20},
    ])
    srep = pd.DataFrame({"tissue": ["Brain_Cortex"], "gene": ["G2"], "phenotype_id": ["rep2"]})
    pc = probe_coloc(pairs, hier, srep).set_index("probe_key")
    assert pc.at["i1", "coloc_PP4_probe"] == pytest.approx(0.81)
    assert pc.at["i1", "coloc_estimator_probe"] == "susie"
    assert "i2" not in pc.index
    assert pc.at["rep2", "coloc_estimator_probe"] == "abf"
    assert pc.at["G2", "coloc_PP4_probe"] == pytest.approx(0.20)


def test_brainseq_refuses_the_mixed_ancestry_arm_and_is_not_yet_wired():
    with pytest.raises(SystemExit, match="ea_only"):
        check_qtl_source("brainseq", "all_samples")
    with pytest.raises(SystemExit, match="not wired"):
        check_qtl_source("brainseq", "ea_only")
    with pytest.raises(SystemExit, match="unknown"):
        check_qtl_source("somewhere")
    assert check_qtl_source("gtex") is None
