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
    LD_MULTI_SNP,
    PEQTL_SMR,
    SMR_MULTI_COL,
    SMR_MULTI_OUT_COLS,
    SMR_MULTI_STATUS,
    SMR_OUT_COLS,
    SMR_STATUS,
    SOURCE_MODALITIES,
    WEAK_F,
    collect_attrition,
    instrumented_family_sizes,
    parse_smr_log,
    multi_family_sizes,
    smr_command,
    smr_multi_status,
    smr_out_suffix,
    smr_status,
    tss_fallback_genes,
    build_brainseq_targets,
    build_targets,
    check_qtl_source,
    classify,
    esd_rows,
    esd_rows_brainseq,
    flist_rows,
    locus_chr,
    ma_from_loci,
    probe_coloc,
    probe_coloc_brainseq,
    probe_gene,
    read_smr,
    run_root,
    smr_command,
    smr_threshold,
    source_root,
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
# Multiple testing: two families, corrected apart
# --------------------------------------------------------------------------- #
def _probe_rows(primary_instrumented: int, secondary_instrumented: int,
                uninstrumented: int = 2) -> pd.DataFrame:
    rows = []
    for i in range(primary_instrumented):
        rows.append({"analysis": "aging__ad", "modality": "sQTL", "tissue": "Brain_Cortex",
                     "probeID": f"P{i}", "probe_family": "primary", "p_SMR": 1e-6})
    for i in range(secondary_instrumented):
        rows.append({"analysis": "aging__ad", "modality": "sQTL", "tissue": "Brain_Cortex",
                     "probeID": f"S{i}", "probe_family": "secondary", "p_SMR": 1e-6})
    for i in range(uninstrumented):
        rows.append({"analysis": "aging__ad", "modality": "sQTL", "tissue": "Brain_Cortex",
                     "probeID": f"U{i}", "probe_family": "primary", "p_SMR": np.nan})
    return pd.DataFrame(rows)


def test_the_primary_denominator_does_not_count_the_secondary_family():
    """The defect this replaced: the two families were pooled into one Bonferroni
    denominator, and then only the primary count was printed, so the reported threshold
    could not be reproduced from the reported n."""
    fam = instrumented_family_sizes(_probe_rows(3, 4)).set_index("probe_family")
    assert int(fam.at["primary", "n_instrumented_family"]) == 3
    assert int(fam.at["secondary", "n_instrumented_family"]) == 4
    assert 7 not in set(fam["n_instrumented_family"])


def test_an_uninstrumented_probe_is_not_in_any_denominator():
    """There is no p-value to correct, so counting it would only make the threshold stricter
    for reasons unrelated to how many tests were done."""
    fam = instrumented_family_sizes(_probe_rows(3, 0, uninstrumented=5))
    assert int(fam.loc[fam["probe_family"] == "primary",
                       "n_instrumented_family"].iloc[0]) == 3


def test_the_same_probe_in_two_tissues_counts_twice():
    d = _probe_rows(1, 0, uninstrumented=0)
    d = pd.concat([d, d.assign(tissue="Brain_Cortex_2")], ignore_index=True)
    fam = instrumented_family_sizes(d)
    assert int(fam["n_instrumented_family"].iloc[0]) == 2


# --------------------------------------------------------------------------- #
# The SMR-side status axis, and the split of `neither`
# --------------------------------------------------------------------------- #
def test_smr_status_separates_untested_from_tested_and_null():
    """`no_instrument` means the probe was never tested; folding it into a null was what
    made `neither` the largest and least readable cell in both reports."""
    assert smr_status(np.nan, np.nan, np.nan, 0.001) == "no_instrument"
    assert smr_status(0.5, 0.4, 20, 0.001) == "instrumented_tested_null"
    assert smr_status(1e-6, 0.4, 20, 0.001) == "smr_heidi_supported"
    assert smr_status(1e-6, 1e-4, 20, 0.001) == "smr_signal_heidi_rejects"
    assert smr_status(1e-6, np.nan, np.nan, 0.001) == "smr_signal_heidi_unavailable"
    assert smr_status(1e-6, 0.4, 2, 0.001) == "smr_signal_heidi_unavailable"


def test_every_status_the_classifier_can_emit_is_declared():
    emitted = {smr_status(*a) for a in [
        (np.nan, np.nan, np.nan, 0.001), (0.5, 0.4, 20, 0.001), (1e-6, 0.4, 20, 0.001),
        (1e-6, 1e-4, 20, 0.001), (1e-6, np.nan, np.nan, 0.001)]}
    assert emitted <= set(SMR_STATUS)
    assert emitted == set(SMR_STATUS)


def test_neither_is_split_by_whether_the_probe_was_ever_instrumented():
    no_coloc = 0.1
    assert classify(no_coloc, np.nan, np.nan, np.nan, 0.001) == "neither_no_instrument"
    assert classify(no_coloc, 0.5, 0.4, 20, 0.001) == "neither_tested_null"
    assert {"neither_no_instrument", "neither_tested_null"} <= set(AGREEMENT)
    assert "neither" not in AGREEMENT


def test_the_agreement_classes_still_track_smr_status_where_coloc_calls():
    coloc = 0.95
    assert classify(coloc, np.nan, np.nan, np.nan, 0.001) == "coloc_no_instrument"
    assert classify(coloc, 0.5, 0.4, 20, 0.001) == "coloc_smr_not_significant"
    assert classify(coloc, 1e-6, 0.4, 20, 0.001) == "coloc_and_smr_heidi_not_rejected"
    assert classify(coloc, 1e-6, 1e-4, 20, 0.001) == "coloc_and_smr_heidi_rejected"
    assert classify(coloc, 1e-6, 0.4, 2, 0.001) == "coloc_and_smr_heidi_untestable"


def test_a_weak_instrument_marker_exists_and_is_only_a_marker():
    """F is reported; nothing in the module may drop a probe for failing it."""
    assert WEAK_F == 10.0


# --------------------------------------------------------------------------- #
# SNP attrition into HEIDI
# --------------------------------------------------------------------------- #
_LOG = """smr --bfile x --gwas-summary y

16212 SNPs to be included from [/p/A_g.chr11.esi].
489 individuals to be included from [/p/1000G.EUR.QC.11.fam].
493922 SNPs to be included from [/p/1000G.EUR.QC.11.bim].
2677 SNPs are included after allele checking.
Genotype data for 489 individuals and 2677 SNPs to be included from [/p/1000G.EUR.QC.11.bed].
5 Probes to be included from [/p/A_g.chr11.epi].
eQTL summary data of 1 Probes to be included from [/p/A_g.chr11.besd].
GWAS summary data of 214144 SNPs to be included from [/p/aging__scz.ma].
"""


def test_the_log_parser_reads_each_harmonization_step(tmp_path):
    f = tmp_path / "A_g.chr11.log"
    f.write_text(_LOG)
    got = parse_smr_log(f)
    assert got["n_besd_snps"] == 16212
    assert got["n_panel_snps"] == 493922
    assert got["n_smr_shared"] == 2677
    assert got["n_gwas_snps"] == 214144
    assert got["n_probes_epi"] == 5
    assert got["n_probes_besd"] == 1


def test_the_smr_shared_count_is_not_named_after_alleles():
    """SMR prints it as "after allele checking", but it is the BESD n panel n GWAS count and
    the GWAS side does essentially all of the cutting. The old name invited reading a
    coverage number as an allele-QC failure, which is exactly what happened on 2026-09-11."""
    from isograph_benchmark.real_data.smr_heidi import _ATTRITION

    assert "n_smr_shared" in _ATTRITION
    assert "n_allele_harmonized" not in _ATTRITION


def test_a_step_the_log_never_printed_stays_nan_rather_than_zero(tmp_path):
    """Zero would assert that everything was dropped; absence means SMR did not report it."""
    f = tmp_path / "S_g.chr3.log"
    f.write_text("16212 SNPs to be included from [/p/S_g.chr3.esi].\n")
    got = parse_smr_log(f)
    assert got["n_besd_snps"] == 16212
    assert np.isnan(got["n_smr_shared"])
    assert np.isnan(got["n_freq_mismatch"])


def test_collect_attrition_keys_runs_by_path(tmp_path):
    """Without the source inputs it reports only what the log carries -- and notably NOT a
    BESD-retention fraction, whose denominator is imputation density rather than quality."""
    d = tmp_path / "smr" / "aging__scz" / "caudate"
    d.mkdir(parents=True)
    (d / "A_g.chr11.log").write_text(_LOG)
    out = collect_attrition(tmp_path)
    assert len(out) == 1
    r = out.iloc[0]
    assert (r["analysis"], r["tissue"], r["modality"], r["chr"]) == (
        "aging__scz", "caudate", "A_g", 11)
    assert r["n_smr_shared"] == 2677
    assert "frac_besd_retained" not in out.columns


def test_collect_attrition_is_empty_when_nothing_has_run(tmp_path):
    assert collect_attrition(tmp_path).empty


# --------------------------------------------------------------------------- #
# Probe-position fallback
# --------------------------------------------------------------------------- #
def test_a_gene_without_an_hg19_tss_is_flagged_as_a_position_fallback(monkeypatch):
    """`flist_rows` centres such a probe on the median ESD position, which moves the cis
    window and can change which SNPs are eligible instruments -- so it is QC, not a detail."""
    import isograph_benchmark.real_data.smr_heidi as smr

    targets = pd.DataFrame({"gene": ["ENSG1", "ENSG2"], "symbol": ["AAA", "BBB"]})
    monkeypatch.setattr(smr, "gene_tss_hg19", lambda g: pd.DataFrame(
        {"gene": ["ENSG1", "ENSG2"], "symbol": ["AAA", "BBB"],
         "strand": ["+", "+"], "tss": [1_000_000, np.nan]}))
    assert tss_fallback_genes(targets) == {"ENSG2"}


def test_no_targets_means_no_fallbacks():
    assert tss_fallback_genes(pd.DataFrame(columns=["gene", "symbol"])) == set()


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
    (np.nan, 0.20, 0.30, 12, "neither_tested_null"),
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


def test_brainseq_refuses_the_mixed_ancestry_arm_and_any_arm_whose_checks_failed():
    with pytest.raises(SystemExit, match="ea_only"):
        check_qtl_source("brainseq", "all_samples", failures=[])
    with pytest.raises(SystemExit, match="checks"):
        check_qtl_source("brainseq", "ea_only", failures=["caudate: sign pin not run"])
    assert check_qtl_source("brainseq", "ea_only", failures=[]) is None
    with pytest.raises(SystemExit, match="unknown"):
        check_qtl_source("somewhere")
    assert check_qtl_source("gtex") is None


def test_source_roots_keep_the_brainseq_arm_apart_from_gtex():
    g, b = source_root("gtex"), source_root("brainseq", "ea_only")
    assert g.name == "gtex" and b.parts[-2:] == ("brainseq", "ea_only")
    assert run_root("brainseq", 1e-5, "ea_only").parent.parent == b
    assert run_root("brainseq", PEQTL_SMR, "ea_only") == b
    with pytest.raises(SystemExit, match="--arm"):
        source_root("brainseq")
    assert SOURCE_MODALITIES == {"gtex": ("eQTL", "sQTL"), "brainseq": ("A_g", "S_g")}


# --------------------------------------------------------------------------- #
# BrainSEQ ESD and targets
# --------------------------------------------------------------------------- #
def _bs_q(rows):
    return pd.DataFrame(rows, columns=["phenotype_id", "variant_id", "af", "slope"]).assign(
        pval_nominal=1e-9, slope_se=0.05)


def _alleles(rows):
    return pd.DataFrame(rows, columns=["variant_id", "bs_ref", "bs_alt"])


def test_brainseq_esd_effect_allele_is_alt_and_the_pin_turns_only_pinnable_genes():
    """tensorQTL's slope is per ALT whichever way the panel codes the SNP, and an S_g slope
    enters the ESD on the discovery axis only when the pin can orient that gene."""
    q = _bs_q([("ENSG1.2", "rs1", 0.3, 0.5), ("ENSG1.2", "rs2", 0.6, -0.2),
               ("ENSG2.1", "rs1", 0.3, 0.4)])
    alleles = _alleles([("rs1", "A", "G"), ("rs2", "C", "T")])
    bim = _bim([(4, "rs1", 0, 89_000_000, "G", "A"), (4, "rs2", 0, 89_000_500, "C", "T")])
    pin = pd.DataFrame({"gene": ["ENSG1", "ENSG2"], "sign": [-1.0, -1.0],
                        "pinnable": [True, False]})
    e = esd_rows_brainseq(q, alleles, bim, pin).set_index(["phenotype_id", "SNP"])
    assert set(ESD_COLS) - {"SNP"} <= set(e.columns)
    assert (e.at[("ENSG1.2", "rs1"), "A1"], e.at[("ENSG1.2", "rs1"), "A2"]) == ("G", "A")
    assert (e.at[("ENSG1.2", "rs2"), "A1"], e.at[("ENSG1.2", "rs2"), "A2"]) == ("T", "C")
    assert e.at[("ENSG1.2", "rs1"), "Beta"] == pytest.approx(-0.5)     # pinned: flipped
    assert e.at[("ENSG1.2", "rs2"), "Beta"] == pytest.approx(0.2)
    assert e.at[("ENSG2.1", "rs1"), "Beta"] == pytest.approx(0.4)      # unpinnable: as mapped
    assert e.at[("ENSG1.2", "rs1"), "Bp"] == 89_000_000                # hg19 panel position
    raw = esd_rows_brainseq(q, alleles, bim).set_index(["phenotype_id", "SNP"])
    assert raw.at[("ENSG1.2", "rs1"), "Beta"] == pytest.approx(0.5)    # A_g: never pinned


def test_brainseq_esd_drops_ambiguous_mismatched_alleleless_and_monomorphic_variants():
    q = _bs_q([("ENSG1.2", "rs1", 0.3, 0.5), ("ENSG1.2", "rs2", 0.3, 0.5),
               ("ENSG1.2", "rs3", 0.3, 0.5), ("ENSG1.2", "rs4", 0.3, 0.5),
               ("ENSG1.2", "rs5", 1.0, 0.5)])
    alleles = _alleles([("rs1", "A", "T"), ("rs2", "A", "G"), ("rs3", "A", "G"),
                        ("rs5", "A", "G")])                             # rs4 not in the pvar
    bim = _bim([(1, "rs1", 0, 10, "A", "T"), (1, "rs2", 0, 20, "C", "T"),
                (1, "rs3", 0, 30, "A", "G"), (1, "rs4", 0, 40, "A", "G"),
                (1, "rs5", 0, 50, "A", "G")])
    assert esd_rows_brainseq(q, alleles, bim)["SNP"].tolist() == ["rs3"]


def test_brainseq_targets_cross_every_nominated_locus_with_every_region():
    k = dict(analysis="aging__lbd", trait="lbd", LOCUS_ID="locus08_chr4")
    gnom = pd.DataFrame([{**k, "gene": "G1", "symbol": "SNCA"}])
    gcells = pd.DataFrame([{**k, "gene": "G1", "tissue": "Brain_Frontal_Cortex_BA9",
                            "PP4_sQTL": 0.9}])
    bnom = pd.DataFrame([{**k, "gene": "G1", "symbol": "SNCA"},
                         {**k, "gene": "G2", "symbol": "MMRN1"}])
    hier = pd.DataFrame([{**k, "gene": "G1", "tissue": "dlpfc", "modality": "S_g",
                          "PP4": 0.88, "estimator": "susie"}])
    t = build_brainseq_targets(
        gnom, gcells, bnom, hier, ("caudate", "dlpfc"),
        {"caudate": "Brain_Caudate_basal_ganglia", "dlpfc": "Brain_Frontal_Cortex_BA9"})
    assert len(t) == 4 and (t["chr"] == 4).all()
    s = t.set_index(["gene", "tissue"])
    assert s.at[("G1", "dlpfc"), "target_source"] == "both"
    assert s.at[("G2", "caudate"), "target_source"] == "brainseq"
    assert s.at[("G2", "caudate"), "symbol"] == "MMRN1"
    assert s.at[("G1", "dlpfc"), "gtex_call_in_matched_tissue"]
    assert not s.at[("G1", "caudate"), "gtex_call_in_matched_tissue"]
    assert s.at[("G1", "dlpfc"), "coloc_PP4_S_g"] == pytest.approx(0.88)
    assert pd.isna(s.at[("G1", "dlpfc"), "coloc_PP4_A_g"])


def test_brainseq_probe_coloc_is_the_hierarchy_cell_for_that_region_and_axis():
    k = dict(analysis="aging__lbd", trait="lbd", LOCUS_ID="locus08_chr4", gene="G1",
             tissue="dlpfc")
    hier = pd.DataFrame([{**k, "modality": "S_g", "PP4": 0.88, "estimator": "susie"},
                         {**k, "modality": "A_g", "PP4": 0.10, "estimator": "abf"}])
    pc = probe_coloc_brainseq(hier).set_index(["tissue", "modality", "probe_key"])
    assert pc.at[("dlpfc", "S_g", "G1"), "coloc_PP4_probe"] == pytest.approx(0.88)
    assert pc.at[("dlpfc", "A_g", "G1"), "coloc_estimator_probe"] == "abf"


# --------------------------------------------------------------------------- #
# Multi-SNP SMR sensitivity arm (--smr-multi)
# --------------------------------------------------------------------------- #
def test_msmr_header_is_the_smr_header_with_p_smr_multi_before_heidi():
    """Pinned to the header SMR 1.4.2 writes (src/SMR_data_p1.cpp:2632), not to memory."""
    assert SMR_MULTI_OUT_COLS == [
        "probeID", "ProbeChr", "Gene", "Probe_bp", "topSNP", "topSNP_chr", "topSNP_bp",
        "A1", "A2", "Freq", "b_GWAS", "se_GWAS", "p_GWAS", "b_eQTL", "se_eQTL", "p_eQTL",
        "b_SMR", "se_SMR", "p_SMR", "p_SMR_multi", "p_HEIDI", "nsnp_HEIDI"]
    # the single-SNP columns survive unchanged and in order
    assert [c for c in SMR_MULTI_OUT_COLS if c != SMR_MULTI_COL] == SMR_OUT_COLS


def test_multi_arm_writes_its_own_directory_and_composes_with_the_threshold_arm():
    primary = run_root("gtex")
    assert run_root("gtex", multi=True) == primary / "sensitivity" / "smr_multi"
    # the pre-specified relaxed arm keeps the exact directory the wrapper documents
    assert run_root("gtex", peqtl_smr=1e-6) == primary / "sensitivity" / "peqtl_smr_1e-06"
    # both at once must not collapse onto either single-arm directory
    both = run_root("gtex", peqtl_smr=1e-6, multi=True)
    assert both == primary / "sensitivity" / "peqtl_smr_1e-06__smr_multi"
    assert both != run_root("gtex", multi=True) and both != run_root("gtex", peqtl_smr=1e-6)
    # a threshold equal to the default is the primary arm, not a sensitivity directory
    assert run_root("gtex", peqtl_smr=PEQTL_SMR) == primary


def test_smr_command_adds_the_multi_flags_only_when_asked():
    args = dict(bfile=Path("b"), ma=Path("m.ma"), besd=Path("x"), probes=Path("p"),
                out=Path("o"))
    single = smr_command(**args)
    assert "--smr-multi" not in single and "--ld-multi-snp" not in single
    multi = smr_command(**args, multi=True)
    assert "--smr-multi" in multi
    assert multi[multi.index("--ld-multi-snp") + 1] == f"{LD_MULTI_SNP:g}"
    # the multi arm changes nothing else about the call
    assert [a for a in multi if a not in ("--smr-multi", "--ld-multi-snp",
                                          f"{LD_MULTI_SNP:g}")] == single
    assert smr_out_suffix(False) == ".smr" and smr_out_suffix(True) == ".msmr"


def test_read_smr_checks_the_header_that_matches_the_extension(tmp_path):
    row = {c: 1 for c in SMR_MULTI_OUT_COLS}
    m = tmp_path / "x.msmr"
    pd.DataFrame([row]).to_csv(m, sep="\t", index=False)
    assert SMR_MULTI_COL in read_smr(m).columns

    # a .smr missing the single-SNP columns must fail loudly rather than parse
    bad = tmp_path / "x.smr"
    pd.DataFrame([{"probeID": "P1"}]).to_csv(bad, sep="\t", index=False)
    with pytest.raises(SystemExit):
        read_smr(bad)


def test_multi_unavailable_is_not_a_null_result():
    thr = 0.05 / 10
    # never instrumented: neither test ran
    assert smr_multi_status(np.nan, np.nan, 0.5, 10, thr) == "no_instrument"
    # instrumented, but too few cis SNPs survived pruning for the multi test
    assert smr_multi_status(1e-9, np.nan, 0.5, 10, thr) == "multi_unavailable"
    # tested and null
    assert smr_multi_status(1e-9, 0.5, 0.5, 10, thr) == "multi_tested_null"
    # significant, and HEIDI decides the rest
    assert smr_multi_status(1e-9, 1e-9, 0.5, 10, thr) == "multi_heidi_supported"
    assert smr_multi_status(1e-9, 1e-9, 1e-4, 10, thr) == "multi_signal_heidi_rejects"
    assert smr_multi_status(1e-9, 1e-9, 0.5, 2, thr) == "multi_signal_heidi_unavailable"
    assert set(SMR_MULTI_STATUS) >= {
        smr_multi_status(1e-9, p, h, n, thr)
        for p, h, n in [(np.nan, 0.5, 10), (0.5, 0.5, 10), (1e-9, 0.5, 10),
                        (1e-9, 1e-4, 10), (1e-9, 0.5, 2)]}


def test_the_multi_family_denominator_counts_only_probes_the_multi_test_ran_on():
    base = dict(analysis="ad", modality="sQTL", probe_family="primary", tissue="cortex")
    d = pd.DataFrame([
        {**base, "probeID": "P1", "p_SMR": 1e-9, SMR_MULTI_COL: 1e-9},
        {**base, "probeID": "P2", "p_SMR": 1e-9, SMR_MULTI_COL: np.nan},  # unavailable
        {**base, "probeID": "P3", "p_SMR": np.nan, SMR_MULTI_COL: np.nan},  # no instrument
    ])
    assert int(instrumented_family_sizes(d)["n_instrumented_family"].iloc[0]) == 2
    assert int(multi_family_sizes(d)["n_multi_family"].iloc[0]) == 1

    # the same probe in two tissues is two tests, in both denominators
    d2 = pd.concat([d, d.assign(tissue="cerebellum")], ignore_index=True)
    assert int(multi_family_sizes(d2)["n_multi_family"].iloc[0]) == 2

    # nothing tested -> an empty frame with the right columns, never a divide by zero
    empty = multi_family_sizes(d.assign(**{SMR_MULTI_COL: np.nan}))
    assert empty.empty and "n_multi_family" in empty.columns
