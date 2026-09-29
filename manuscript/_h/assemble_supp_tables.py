"""Assemble real-data supplementary tables from the analysis parquet ledgers.

Reads the committed result parquets under 03_module_characterization/_m and 04_module_trust/_m/stability,
emits one clean CSV per supplementary table under manuscript/_m/supp_tables/, and a
machine-checkable manifest. Presentation only; every number is copied verbatim from
the source ledgers -- regenerate, do not hand-edit.

Run: python manuscript/_h/assemble_supp_tables.py   (login node, no SLURM)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from isograph_benchmark.paths import ensure_dir, region_store, stage_out  # noqa: E402

OUT = ensure_dir(stage_out("manuscript", "supp_tables"))
QTL = stage_out("anchoring", "qtl_anchoring_meta")
BASE = stage_out("trust", "baseline_comparison")
TRUST = stage_out("trust.stability", "module_trust")
TRUST_TABLES = stage_out("trust.stability", "module_trust_tables")
DTU = stage_out("trust", "dtu_added_value")           # module context, held-out DTU
GATE = region_store("brainseq", "caudate_sczd")
COMP = stage_out("characterize")                      # composition adjustment (BrainSEQ)
CHAR = stage_out("characterize")                      # axis separation rollup
COMP_GTEX = stage_out("modules", "gtex", "_m", "composition")
MECH = stage_out("mechanism")                         # switch-mechanism stage
RBP = stage_out("regulation", "rbp")
ANCHOR = stage_out("anchoring", "module_genetic_anchoring_meta")
SCZ = stage_out("integration", "scz_age_projection")


def write(df: pd.DataFrame, name: str) -> None:
    df.to_csv(OUT / name, index=False)
    print(f"  wrote {name:32s} {df.shape[0]:>4d} x {df.shape[1]}")


def round_num(df: pd.DataFrame, n: int = 4) -> pd.DataFrame:
    """Round numeric columns to n decimals, but p-values to 3 significant figures
    so small p's (e.g. 8e-6) survive instead of collapsing to 0.0."""
    df = df.copy()
    p_like = ("p_value", "p_fe", "p_re", "pheno_fdr", "mwu_p", "adjusted_p")
    for col in df.select_dtypes("number").columns:
        if any(tok in col for tok in p_like):
            df[col] = df[col].apply(
                lambda x: float(f"{x:.3g}") if pd.notna(x) else x)
        else:
            df[col] = df[col].round(n)
    return df


# --- S1 / S2  three-baseline comparison ------------------------------------
def baseline_tables() -> None:
    pooled = pd.read_parquet(BASE / "baseline_comparison_pooled.parquet")
    pooled = pooled[
        ["method", "features", "n_regions", "median_n_modules", "median_module_size",
         "mean_frac_pheno_sig", "mean_frac_both", "mean_frac_go_enriched",
         "total_pheno_sig", "total_both"]
    ]
    write(round_num(pooled), "tableS1_baseline_pooled.csv")

    per = pd.read_parquet(BASE / "baseline_comparison.parquet")
    per = per.sort_values(["cohort", "region", "method"]).reset_index(drop=True)
    write(round_num(per), "tableS2_baseline_per_region.csv")


# --- S3 / S4  QTL splicing-specificity contrast ----------------------------
def _tidy_contrast(df: pd.DataFrame) -> pd.DataFrame:
    df = df[["graph_method", "module_set", "k", "ratio_fe", "ratio_fe_low",
             "ratio_fe_high", "p_fe", "I2"]].copy()
    df.columns = ["method", "module_set", "n_analyses", "sqtl_eqtl_ratio",
                  "ratio_ci_low", "ratio_ci_high", "p_value", "I2"]
    order = {"all_modules": 0, "pheno_sig_modules": 1, "go_invisible_modules": 2,
             "go_visible_modules": 3}
    mord = {"isograph": 0, "wgcna_switch_only": 1, "wgcna_multiplex": 2}
    df["_o"] = df.module_set.map(order)
    df["_m"] = df.method.map(mord)
    df = df.sort_values(["_o", "_m"]).drop(columns=["_o", "_m"]).reset_index(drop=True)
    return df


def qtl_tables() -> None:
    contrast = _tidy_contrast(pd.read_parquet(QTL / "qtl_anchoring_meta_contrast.parquet"))
    write(round_num(contrast), "tableS3_qtl_specificity_contrast.csv")

    common = _tidy_contrast(pd.read_parquet(QTL / "qtl_anchoring_meta_contrast_common.parquet"))
    write(round_num(common), "tableS4_qtl_specificity_matched_baseline.csv")

    raw = pd.read_parquet(QTL / "qtl_anchoring_meta.parquet")
    raw = raw[raw.k > 0][["graph_method", "xqtl_kind", "module_set", "k",
                          "n_fg_total", "or_fe", "or_fe_low", "or_fe_high",
                          "p_fe", "I2"]].copy()
    raw.columns = ["method", "xqtl_kind", "module_set", "n_analyses", "n_fg_total",
                   "or", "or_ci_low", "or_ci_high", "p_value", "I2"]
    write(round_num(raw), "tableS5_qtl_raw_enrichment_or.csv")


# --- S6  GO-invisible disease switch modules -------------------------------
def gate_table() -> None:
    import json
    g = pd.read_parquet(GATE / "go_invisible_gate.parquet").copy()
    bg = json.load(open(GATE / "go_invisible_gate_background.json"))
    cols = ["module", "n_genes", "pheno_fdr", "go_invisible", "n_go_terms",
            "genes_with_real_switch", "max_switch_strength", "n_sig_switch_tx",
            "frac_cds_changed", "frac_coding_status_change", "frac_biotype_switch",
            "frac_utr_changed", "top_switch_genes"]
    g = g[cols]
    bg_row = {c: bg.get(c, "") for c in cols}
    bg_row["module"] = "_background (all switch tx)"
    g = pd.concat([g, pd.DataFrame([bg_row])], ignore_index=True)
    write(round_num(g), "tableS6_go_invisible_gate.csv")


# --- S7  module trust ledgers -----------------------------------------------
# All five S7 tables are copied from the CSVs written by
# `isograph_benchmark/real_data/module_trust_tables.py` (stage 04 `_h/04e`) rather than
# re-derived from the parquets here. One assembler owning the arithmetic is what keeps the
# figure, the tables and the Results text quoting the same numbers.
def trust_table() -> None:
    """S7 -- per-region funnel, and S7a-S7d, the per-module ledgers behind it.

    S7 carries BOTH methods. The reproducibility subsection is a matched contrast
    throughout, and an IsoGraph-only table cannot be checked against the baseline the text
    quotes beside every claim.
    """
    if not TRUST_TABLES.exists():
        print("  skip tableS7*: run 04_module_trust/_h/04e.module_trust_tables.sh first")
        return

    def _read(name: str) -> pd.DataFrame:
        return pd.read_csv(TRUST_TABLES / f"{name}.csv")

    write(round_num(_read("region_funnel")), "tableS7_module_trust_funnel.csv")

    stab = _read("split_half_modules")
    cols = ["cohort", "region", "method", "module_id", "n_genes", "n_genes_assigned",
            "coassign_density", "null_mean", "best_match_jaccard", "perm_p", "fdr",
            "trusted"]
    write(round_num(stab[[c for c in cols if c in stab.columns]]),
          "tableS7a_split_half_module_ledger.csv")

    proj = _read("projection_modules")
    # raw_* stay in: they are the evidence that the raw projected age correlation is a
    # property of the target cohort, which is why the reported statistic is standardised.
    pcols = ["pair", "method", "direction", "module_id", "status", "n_features",
             "n_switch_features", "signed_kme", "kme_null_mean", "kme_perm_p",
             "kme_perm_q", "preservation_r_native_pc1", "age_r_source", "age_r_target",
             "age_z_source", "age_z_target", "sign_match", "both_sig", "raw_sign_match",
             "raw_both_sig"]
    write(round_num(proj[[c for c in pcols if c in proj.columns]]),
          "tableS7b_projection_module_ledger.csv")

    perm = _read("crosscohort_permutation")
    write(round_num(perm), "tableS7c_crosscohort_permutation.csv")

    func = _read("functional_preservation")
    write(round_num(func), "tableS7d_functional_preservation.csv")

    # S7e/S7f answer the two questions the subsection's claims invite and the other
    # tables cannot: whether split-half agreement is a property of the chosen Leiden
    # resolution, and why the projected-age statistic is standardised rather than raw.
    res_f = TRUST_TABLES / "resolution_sensitivity.csv"
    if res_f.exists():
        write(round_num(_read("resolution_sensitivity")),
              "tableS7e_resolution_sensitivity.csv")
    else:
        print("  skip tableS7e: no resolution_sensitivity.csv (sweep not run)")

    write(round_num(_read("projection_sign_scale")),
          "tableS7f_projection_sign_scale.csv")


# --- S13a  switch vs abundance axis separation ------------------------------
def separation_table() -> None:
    """Per-analysis separation of the switch and abundance coordinates.

    The ledger behind figSeparation panel a and the median |r| range the manuscript
    quotes. All 17 analyses are here, including the five GTEx regions that the
    composition table (S13) cannot cover for want of a matched snRNA reference: the
    separation test needs no deconvolution, so this is the one per-analysis table that
    spans the whole set.

    Quartiles and the three tail fractions are kept beside the median, because the claim
    is about the bulk of the distribution ("largely distinct"), not about a central value:
    a median |r| near 0.12 with 40% of genes under 0.1 and 4% over 0.5 says something a
    median alone does not.
    """
    f = CHAR / "axis_orthogonality_summary.parquet"
    if not f.exists():
        print("  skip tableS13a: no axis_orthogonality_summary.parquet")
        return
    df = pd.read_parquet(f)
    write(round_num(df), "tableS13a_axis_orthogonality.csv")


# --- S13  cell-type composition adjustment ----------------------------------
def composition_table() -> None:
    """Base vs composition-adjusted switch-unique counts, both cohorts.

    This is the ledger behind figCompositionRobustness. Both arms are reported
    together, and the five GTEx regions without a defensibly matched snRNA reference
    are absent by design (not deconvolved rather than forced against a mismatched
    panel) -- the figure names them; this table covers only what was tested.

    `n_overlap`/`n_new` are the gene-level turnover the counts alone hide: adjustment
    both drops and adds genes, so the adjusted set is not nested in the unadjusted one
    and no ratio of the two counts is reported.
    """
    bs = pd.read_parquet(COMP / "composition_adjustment.parquet")
    bs.insert(0, "cohort", "BrainSEQ")
    gt = pd.read_parquet(COMP_GTEX / "composition_adjustment_gtex.parquet")
    gt.insert(0, "cohort", "GTEx")
    df = pd.concat([bs, gt], ignore_index=True)
    # The stage writes these columns as `comp_unique_*`, from before the criterion was
    # renamed: the manuscript calls these genes switch-unique, because "composition-unique"
    # reads as a statement about cell composition rather than about isoform usage. Renamed
    # here, in the display layer, so the published table matches the text without touching
    # the canonical column names the analysis stages share.
    df = df.rename(columns={"comp_unique_base": "switch_unique_base",
                            "comp_unique_adj": "switch_unique_adj"})
    df = _with_marker_denominator(df)
    write(round_num(df), "tableS13_composition_adjustment.csv")


# Stage directory -> the `region` label the composition ledgers use. GTEx labels are
# "GTEx <region dir>"; BrainSEQ's are hand-named, so they are spelled out.
_BRAINSEQ_REGION = {"caudate": "aging caudate", "caudate_sczd": "SCZD (caudate)",
                    "dlpfc": "aging DLPFC", "hippocampus": "aging hippocampus"}


def _with_marker_denominator(df: pd.DataFrame) -> pd.DataFrame:
    """Add `n_marker_tests`, the module x cell-type tests behind the marker counts.

    The ledgers carry `marker_enriched` / `marker_depleted` as bare counts. A count with no
    denominator cannot be read -- the number of tests differs ~3x between regions with the
    module count and the reference panel -- so the denominator is counted here from the
    same per-region files, and the counts are re-derived from them as a consistency check
    rather than trusted.
    """
    rows = []
    for f in sorted(stage_out("modules").glob(
            "*/*/_m/isograph_vae/celltype_composition/marker_enrichment.parquet")):
        cohort_dir, region_dir = f.parts[-6], f.parts[-5]
        region = (_BRAINSEQ_REGION[region_dir] if cohort_dir == "brainseq"
                  else f"GTEx {region_dir}")
        mk = pd.read_parquet(f)
        rows.append(dict(region=region, n_marker_tests=len(mk),
                         _enr=int((mk["fdr_enrich"] < 0.05).sum()),
                         _dep=int((mk["fdr_deplete"] < 0.05).sum())))
    if not rows:
        print("  note: no marker_enrichment.parquet found; n_marker_tests left empty")
        return df.assign(n_marker_tests=pd.NA)
    mk = pd.DataFrame(rows)
    out = df.merge(mk, on="region", how="left", validate="one_to_one")
    have = out["n_marker_tests"].notna()
    bad = have & ((out["_enr"] != out["marker_enriched"])
                  | (out["_dep"] != out["marker_depleted"]))
    if bad.any():
        raise SystemExit("marker counts disagree with marker_enrichment.parquet for: "
                         + ", ".join(out.loc[bad, "region"]))
    out = out.drop(columns=["_enr", "_dep"])
    # Denominator directly before the two counts it belongs to.
    cols = [c for c in out.columns if c != "n_marker_tests"]
    i = cols.index("marker_depleted")
    return out[cols[:i] + ["n_marker_tests"] + cols[i:]]


# --- S21  module context and held-out DTU evidence --------------------------
def dtu_context_table() -> None:
    """One row per split-half analysis: the held-out DTU test, before and after the
    co-expression comparator, the sub-threshold reading and the giant-module check.

    The two permutation q values are BH across the six analyses, as in the source
    summary. The giant-module columns are post hoc (means over the ten directions).
    """
    s = pd.read_parquet(DTU / "region_summary.parquet")
    g = pd.read_parquet(DTU / "giant_module_sensitivity.parquet")
    giant = g.groupby(["cohort", "region"], as_index=False).agg(
        r_without_giant_modules=("r_no_giant", "mean"),
        r_without_largest_module=("r_no_largest", "mean"),
        n_giant_modules_median=("n_giant_modules", "median"),
        frac_tested_in_largest_median=("frac_tested_in_largest", "median"))
    cols = {
        "cohort": "cohort", "region": "region", "n_replicates": "n_split_directions",
        "n_genes_median": "n_genes_median", "coverage": "coverage",
        "n_modules_median": "n_modules_median",
        "r_module": "r_module_mean", "r_module_median": "r_module_median",
        "r_module_min": "r_module_min", "r_module_max": "r_module_max",
        "r_module_n_positive": "r_module_n_positive",
        "p_strat": "r_module_perm_p_value", "q_strat": "r_module_perm_q_value",
        "r_module_given_abund": "r_given_coexpr_mean",
        "r_module_given_abund_median": "r_given_coexpr_median",
        "r_module_given_abund_min": "r_given_coexpr_min",
        "r_module_given_abund_max": "r_given_coexpr_max",
        "r_module_given_abund_n_positive": "r_given_coexpr_n_positive",
        "p_given_abund_strat": "r_given_coexpr_perm_p_value",
        "q_given_abund_strat": "r_given_coexpr_perm_q_value",
        "r_abund": "r_coexpr_context_alone_mean",
        "region_class": "pattern",
        "n_subthreshold": "n_subthreshold_genes",
        "subthr_rate_top": "subthreshold_replication_top_tertile",
        "subthr_rate_bottom": "subthreshold_replication_bottom_tertile",
        "subthr_or_median": "subthreshold_adjusted_or_median",
        "subthr_or_n_above_1": "subthreshold_or_n_above_1",
    }
    d = s[list(cols)].rename(columns=cols).merge(giant, on=["cohort", "region"], how="left")
    write(round_num(d), "tableS21_module_context_heldout_dtu.csv")


# --- S14  long-read orthogonal confirmation ---------------------------------
def orthogonal_table() -> None:
    """Per-gene long-read confirmation of the genetically anchored switch pairs.

    `max_anchored_if` is the qualification that decides whether a gene's usage
    correlation is interpretable at all; genes below the 0.05 bar are retained in the
    table with confirmed_at_usable_abundance = False rather than dropped.
    """
    d = pd.read_parquet(
        MECH / "switch_orthogonal_confirm" / "anchored_gene_confirmation.parquet")
    write(round_num(d), "tableS14_longread_orthogonal_confirmation.csv")


# --- S15  independent-caller (satuRn) concordance ---------------------------
def isa_table() -> None:
    """satuRn / IsoformSwitchAnalyzeR concordance across all 17 analyses.

    Analyses with zero switch genes are VACUOUS by construction, not null results;
    they are kept as rows with n_switch_genes = 0 and NaN statistics so a reader
    cannot mistake the testable denominator.
    """
    import json

    rows = []
    root = MECH / "isa_concordance"
    for d in sorted(x for x in root.iterdir() if x.is_dir()):
        f = d / "summary.json"
        if not f.exists():
            continue
        # The CLI writes bare NaN, which json.loads accepts but strict parsers do not.
        j = json.loads(f.read_text())
        rows.append({
            "analysis": d.name, "cohort": j.get("cohort"), "region": j.get("region"),
            "trait": j.get("trait"), "caller": j.get("caller"),
            "n_genes": j.get("n_genes"), "n_switch_genes": j.get("n_switch_genes"),
            "n_background_genes": j.get("n_background_genes"),
            "median_evidence_switch": j.get("median_evidence_switch"),
            "median_evidence_background": j.get("median_evidence_background"),
            "rank_biserial": j.get("rank_biserial"), "mwu_p": j.get("mwu_p"),
            "adjusted_or": j.get("adjusted_or"), "adjusted_p": j.get("adjusted_p"),
        })
    df = pd.DataFrame(rows).sort_values(["cohort", "region", "trait"])
    write(round_num(df.reset_index(drop=True)), "tableS15_isa_concordance.csv")


# --- S16  module-level genetic anchoring (deliberately a table, not a figure) --
def module_anchoring_table() -> None:
    """Per-module splicing-specificity contrast vs a size-matched permutation null.

    Reported as a table ON PURPOSE. The result is that splicing anchoring is a
    POOLED-gene property, not a per-module one, and it is NOT concentrated in
    GO-invisible modules -- plotting it would invite exactly the reading the analysis
    rules out. The table names the exemplar anchored programs and nothing more.
    """
    f = ANCHOR / "module_genetic_anchoring_all.parquet"
    if not f.exists():
        print(f"  skip tableS16: {f.name} absent")
        return
    d = pd.read_parquet(f)
    if "perm_p" in d.columns:
        d = d.sort_values("perm_p")
    write(round_num(d.reset_index(drop=True)), "tableS16_module_genetic_anchoring.csv")


# --- S17  per-RBP eCLIP binding support -------------------------------------
def rbp_binding_table() -> None:
    """Per-RBP within-gene eCLIP contrast: switched vs constitutive exons.

    One two-sided exact McNemar per RBP over its unique nominated regulon genes
    (deduped across regions), BH across the testable RBPs -- NOT one test per
    region x module, which would count the same region-invariant binding fact once per
    region. This is binding CAPACITY at alternative-exon sequence in HepG2/K562, not
    neuronal occupancy of these regulons.
    """
    d = pd.read_parquet(RBP / "rbp_binding_support.parquet")
    d = d.sort_values("rate_diff", ascending=False).reset_index(drop=True)
    write(round_num(d), "tableS17_rbp_eclip_binding_support.csv")


# --- S18  SCZ-risk convergence on age-sensitive switch modules ---------------
def scz_convergence_table() -> None:
    """Per-module SCZ colocalized-locus counts behind Fig 4E.

    Carries the candidate RBP regulators that the folded half-width Fig 4E panel
    cannot show legibly.
    """
    f = SCZ / "convergence.parquet"
    if not f.exists():
        print(f"  skip tableS18: {f.name} absent")
        return
    d = pd.read_parquet(f).sort_values("n_coloc_loci", ascending=False)
    write(round_num(d.reset_index(drop=True)), "tableS18_scz_convergence.csv")


def qtl_sensitivity_table() -> None:
    """S19 — the three pre-specified sensitivity arms for the anchoring contrast.

    The reviewer question this answers is not "is the contrast significant" but "is it
    an artefact of (a) selective constraint, (b) the sGene/eGene threshold, or (c) a
    binary call that cannot count independent signals". One row per (arm, graph_method,
    covariate_set, module_set) with the pooled ratio, so the primary arm and each
    sensitivity sit in one table and can be read against each other.

    The constraint arm carries BOTH covariate sets fitted on the identical
    constraint-complete gene subset — the adjusted and unadjusted rows there differ
    only by the covariates, which is what makes it a nested-model test rather than a
    change of universe.
    """
    base = stage_out("anchoring", "qtl_anchoring_meta")
    # Arm labels stay comma-free and match the sensitivity subdirectory naming, so the
    # CSV needs no quoting and a reader can map a row straight back to its output dir.
    arms = [("binary_standard_PRIMARY", base / "qtl_anchoring_meta_contrast.parquet")]
    sens = base / "sensitivity"
    if sens.exists():
        for d in sorted(p for p in sens.iterdir() if p.is_dir()):
            arms.append((d.name, d / "qtl_anchoring_meta_contrast.parquet"))
    parts = []
    for label, f in arms:
        if not f.exists():
            print(f"  tableS19: {label} absent ({f.name}) — omitted")
            continue
        d = pd.read_parquet(f)
        d.insert(0, "arm", label)
        if "covariate_set" not in d.columns:
            d["covariate_set"] = "standard"
        parts.append(d)
    if not parts:
        print("  skip tableS19: no anchoring contrast tables found")
        return
    out = pd.concat(parts, ignore_index=True)
    keep = [c for c in ("arm", "graph_method", "covariate_set", "module_set", "k",
                        "ratio_fe", "ratio_fe_low", "ratio_fe_high", "p_fe",
                        "ratio_re", "I2") if c in out.columns]
    write(round_num(out[keep]), "tableS19_qtl_anchoring_sensitivity.csv")


def coloc_convergence_table() -> None:
    """S20 — module-level coloc convergence for all five traits.

    Two sheets' worth of content in one table: the per-(trait, source) global test and
    the per-module counts. The `testable` flag and the two background columns are the
    load-bearing parts — a reader who takes only `n_coloc_genes` will over-read it, which
    is exactly what the count-only SCZ panel invited.
    """
    d = stage_out("anchoring", "module_coloc_convergence")
    gf, cf = d / "global.parquet", d / "convergence.parquet"
    if not (gf.exists() and cf.exists()):
        print("  skip tableS20: module_coloc_convergence outputs absent")
        return
    g = pd.read_parquet(gf)
    gcols = [c for c in ("trait", "source", "n_pool", "n_coloc_genes",
                         "n_modules_with_coloc", "concentration_obs",
                         "concentration_null_mean", "concentration_p",
                         "n_anchored_modules", "frac_coloc_anchored",
                         "frac_pool_anchored", "anchored_hyperg_p",
                         "frac_allgenes_anchored",
                         "anchored_hyperg_p_allgenes_denom") if c in g.columns]
    write(round_num(g[gcols]), "tableS20a_coloc_convergence_global.csv")
    c = pd.read_parquet(cf)
    ccols = [c_ for c_ in ("trait", "source", "module_id", "testable",
                           "n_module_genes_in_pool", "n_coloc_genes", "n_coloc_loci",
                           "n_go_invisible", "hyperg_p", "hyperg_fdr", "loo_worst_p",
                           "magma_p", "coloc_genes") if c_ in c.columns]
    write(round_num(c[ccols]), "tableS20b_coloc_convergence_per_module.csv")


def main() -> None:
    print(f"Writing supplementary tables to {OUT}")
    baseline_tables()
    qtl_tables()
    qtl_sensitivity_table()
    coloc_convergence_table()
    gate_table()
    trust_table()
    separation_table()
    composition_table()
    dtu_context_table()
    orthogonal_table()
    isa_table()
    module_anchoring_table()
    rbp_binding_table()
    scz_convergence_table()
    print("done.")


if __name__ == "__main__":
    main()
