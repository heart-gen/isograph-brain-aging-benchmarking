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
from isograph_benchmark.paths import ensure_dir, region_store, rel, stage_out  # noqa: E402

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
ASE = stage_out("mechanism", "ase_junction_switch")    # within-donor allelic test
ASE_REGIONS = ("caudate", "hippocampus", "dlpfc")
COLOC = stage_out("anchoring", "coloc")                # CLPP isoform events
SUSIE_ALL = stage_out("anchoring", "coloc_signal_susie", "all_introns")
LDSC = stage_out("anchoring", "ldsc")


def write(df: pd.DataFrame, name: str) -> None:
    df.to_csv(OUT / name, index=False)
    print(f"  wrote {name:32s} {df.shape[0]:>4d} x {df.shape[1]}")


def round_num(df: pd.DataFrame, n: int = 4, sig_cols: tuple[str, ...] = ()) -> pd.DataFrame:
    """Round numeric columns to n decimals, but p-values to 3 significant figures
    so small p's (e.g. 8e-6) survive instead of collapsing to 0.0. `sig_cols` names
    extra columns to keep at 3 significant figures (p-values the default tokens miss,
    and tiny coefficients such as LDSC tau)."""
    df = df.copy()
    p_like = ("p_value", "p_fe", "p_re", "pheno_fdr", "mwu_p", "adjusted_p")
    for col in df.select_dtypes("number").columns:
        if col in sig_cols or any(tok in col for tok in p_like):
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


def allelic_tables() -> None:
    """S22 / S23 — within-donor allelic test of isoform choice (BrainSEQ recount).

    S22 is one row per region for the BH-tested gate family, with the calibration the
    claim rests on (homozygous-at-lead null, lambda_GC), the agreement checks (score test,
    between-donor sign) and the module-level cis-control correlation. S23 is the
    per-module ledger behind that correlation. Read straight from the stage summaries.
    """
    import json

    cis = json.loads((ASE / "module_cis_control_summary.json").read_text())
    rows = []
    for region in ASE_REGIONS:
        s = json.loads((ASE / region / "allelic_summary.json").read_text())
        g, c = s["families"]["gate_family"], s["families"]["coloc_nominated"]
        m, md = cis["regions"].get(region, {}), s["module_direction"]
        rows.append({
            "region": region,
            "pairs": g["n_pairs"], "pairs_fitted": g["n_fitted"],
            "genes_fitted": g["n_genes_fitted"], "pairs_nominal_p05": g["n_nominal_p05"],
            "pairs_q05": g["n_q05"], "genes_q05": g["n_genes_q05"],
            "sign_agree_between_all": g["sign_agreement_between"],
            "sign_agree_between_q05": g["sign_agreement_between_q05"],
            "not_converged": g["n_not_converged"], "at_bound": g["n_at_bound"],
            "at_bound_q05": g["n_at_bound_q05"],
            "hom_null_fits": s["hom_null"]["n_fitted"],
            "hom_null_frac_p05": s["hom_null"]["frac_p05"],
            "hom_null_lambda_gc": s["hom_null"]["lambda_gc"],
            "het_lambda_gc": s["lambda_gc_het"],
            "score_vs_glmm_spearman": s["score_vs_glmm_spearman"],
            "median_lead_site_kb": s["median_lead_site_kb"],
            "module_axis_genes": md["n_genes"],
            "module_axis_median_abs_rho": md["median_abs_rho"],
            "module_axis_null_median": md["null_median_abs_rho"],
            "module_axis_perm_p": md["perm_p"],
            "module_cis_modules_ranked": m.get("n_modules_ranked"),
            "module_cis_rho": m.get("rho"),
            "module_cis_perm_p": m.get("perm_p_two_sided"),
            "coloc_nominated_fitted": c["n_fitted"],
            "coloc_nominated_nominal_p05": c["n_nominal_p05"],
        })
    pooled = cis["pooled"]
    rows.append({"region": "pooled (" + ", ".join(pooled["regions"]) + ")",
                 "module_cis_rho": pooled["mean_rho"],
                 "module_cis_perm_p": pooled["perm_p_two_sided"]})
    out = pd.DataFrame(rows)
    counts = [c for c in out.columns if c.startswith(("pairs", "genes", "hom_null_fits",
                                                      "not_converged", "at_bound",
                                                      "module_axis_genes",
                                                      "module_cis_modules",
                                                      "coloc_nominated"))]
    out[counts] = out[counts].astype("Int64")
    write(round_num(out, sig_cols=("module_axis_perm_p", "module_cis_perm_p")),
          "tableS22_allelic_imbalance_regions.csv")

    mod = pd.read_parquet(ASE / "module_cis_control.parquet")
    mod = mod.rename(columns={"genes_cis": "genes_cis_q05", "rate": "cis_rate",
                              "aging_rank": "age_rank_within_trait",
                              "trait_n": "modules_in_trait"})
    write(round_num(mod), "tableS23_module_cis_control.csv")


def clpp_event_table() -> None:
    """S24 — every CLPP-nominated sQTL/eQTL isoform event, with whether its junction
    maps into the tissue-matched IsoGraph switch pair (the 76 events / 30 genes).

    `structural_consequence` is left out on purpose: it flags each transcript against a
    gene reference, so almost every mapped event carries every flag.
    """
    e = pd.read_parquet(COLOC / "coloc_isoform_events_combined.parquet")
    cols = ["analysis", "trait", "case", "gene_name", "gene", "kind", "tissue",
            "best_rsid", "risk_allele", "risk_qtl_effect", "junction", "clpp",
            "go_invisible", "junction_in_switch_pair", "switch_pair", "n_switch_pairs",
            "junction_transcripts", "junction_transcript_polarity_r", "brainseq_region",
            "brainseq_switch_pair", "brainseq_replicates_switch"]
    e = (e[cols].sort_values(["junction_in_switch_pair", "clpp"], ascending=[False, False])
         .reset_index(drop=True))
    write(round_num(e), "tableS24_clpp_isoform_events.csv")


def signal_coloc_tables() -> None:
    """S25 / S26 — signal-level colocalization, all-introns arm.

    S25: every analysis x gene with PP4_sQTL >= 0.8, with its eQTL PP4, tissue
    consistency, and the estimator and prior robustness of the headline cell
    (coloc.susie where both traits could be fine-mapped, coloc.abf otherwise). The aging
    and SCZD schizophrenia gene sets overlap, so rows are analysis x gene; collapse on
    (gene, trait) for gene-by-trait cells. `sqtl_preferential` (PP4_eQTL < 0.5) is a
    selection on the sQTL call, not a splicing-specificity estimate -- S26 is the
    unbiased paired contrast.
    """
    g = pd.read_parquet(SUSIE_ALL / "genes.parquet")
    h = pd.read_parquet(SUSIE_ALL / "cells_hierarchy.parquet",
                        columns=["analysis", "gene", "modality", "tissue", "phenotype_id",
                                 "PP4", "estimator", "prior_robustness"])
    best = (h[h.modality == "sQTL"].sort_values("PP4", ascending=False)
            .drop_duplicates(["analysis", "gene"]))
    s = g[g.PP4_sQTL >= 0.8].merge(best, on=["analysis", "gene"], how="left")
    assert (s.tissue == s.max_tissue_sQTL).all(), "headline cell is not the max tissue"
    s["sqtl_preferential"] = s.PP4_eQTL < 0.5
    s = s.rename(columns={"phenotype_id": "headline_intron",
                          "estimator": "headline_estimator"})
    s["headline_estimator"] = s.headline_estimator.map(
        {"susie": "coloc.susie", "abf": "coloc.abf"})
    cols = ["analysis", "trait", "symbol", "gene", "module_id", "go_invisible",
            "PP4_sQTL", "PP4_eQTL", "sqtl_preferential", "n_tissue",
            "n_tissue_sQTL_coloc", "n_tissue_eQTL_coloc", "max_tissue_sQTL",
            "headline_intron", "headline_estimator", "prior_robustness"]
    s = s[cols].sort_values("PP4_sQTL", ascending=False).reset_index(drop=True)
    write(round_num(s), "tableS25_signal_coloc_nominations.csv")

    c = pd.read_parquet(SUSIE_ALL / "contrast.parquet")
    write(round_num(c, sig_cols=("mcnemar_p", "wilcoxon_cond_p", "wilcoxon_pp4_p")),
          "tableS26_coloc_modality_contrast.csv")


def coloc_event_table() -> None:
    """S28 — every colocalization event (nomination x tissue) of the all-introns arm,
    mapped to transcripts and to the tissue-matched IsoGraph switch pair. Fig 5b,c.

    `signal_level_nomination` is True when coloc.susie scored the gene x trait in at least
    one tissue (33 of 119); the rest rest on the coloc.abf fallback, whose rows name GTEx's
    representative intron only. `clpp_resolves_switch` marks gene x trait cells that the
    CLPP layer (S24) also resolves to a switch pair. No risk-allele direction here: that is
    resolved only in the CLPP layer. `structural_consequence` is left out for the reason
    given under S24.
    """
    e = pd.read_parquet(SUSIE_ALL / "coloc_isoform_events.parquet")
    e["junction_in_switch_pair"] = e.junction_in_switch_pair.fillna(False).astype(bool)
    key = ["gene_name", "trait"]
    nom = e.groupby(key).estimator.agg(lambda s: (s == "susie").any())
    e = e.join(nom.rename("signal_level_nomination"), on=key)
    c = pd.read_parquet(COLOC / "coloc_isoform_events_combined.parquet")
    c = c[c.junction_in_switch_pair.fillna(False).astype(bool)]
    ck = set(zip(c.gene_name, c.trait))
    e["clpp_resolves_switch"] = [k in ck for k in zip(e.gene_name, e.trait)]
    assert e.groupby(key).ngroups == 119 and nom.sum() == 33
    e["estimator"] = e.estimator.map({"susie": "coloc.susie", "abf": "coloc.abf"})
    cols = ["analysis", "trait", "gene_name", "gene", "tissue", "signal_level_nomination",
            "estimator", "PP4_sQTL", "PP4_eQTL", "prior_robustness", "fallback_reason",
            "phenotype_id", "junction", "junction_transcripts", "go_invisible",
            "junction_in_switch_pair", "switch_pair", "n_switch_pairs",
            "junction_transcript_polarity_r", "brainseq_region", "brainseq_switch_pair",
            "brainseq_replicates_switch", "switch_pair_scope", "clpp_resolves_switch"]
    e = (e[cols].sort_values(["junction_in_switch_pair", "signal_level_nomination",
                              "PP4_sQTL"], ascending=[False, False, False])
         .reset_index(drop=True))
    write(round_num(e), "tableS28_coloc_isoform_events.csv")


def longread_coloc_table() -> None:
    """S29 — per-gene long-read confirmation of the colocalization-anchored switch pairs.

    The coloc-layer counterpart of S14 (CLPP layer), behind Fig 4a,b: the 55 genes whose
    colocalizing junction is carried by a tissue-matched switch-pair isoform. Genes whose
    anchored isoform never reaches 0.05 mean isoform fraction stay in the table with
    confirmed_at_usable_abundance = False.
    """
    import json

    oc = MECH / "switch_orthogonal_confirm" / "signal_coloc"
    d = pd.read_parquet(oc / "anchored_gene_confirmation.parquet")
    assert len(d) == 55, len(d)
    # Distinct transcript pairs are the primary unit (a pair anchored from both sides is two
    # orientations with one outcome); the orientation counts stay as the sensitivity columns.
    dp = pd.read_parquet(oc / "distinct_pair_confirmation.parquet")
    det = dp[dp["pair_detected"]]
    per = pd.DataFrame({
        "n_distinct_pairs": dp.groupby("gene").size(),
        "n_distinct_detected": det.groupby("gene").size(),
        "n_distinct_switch_like": det.groupby("gene")["switch_like"].sum(),
        "n_distinct_detected_usable": det[det["anchored_usable"]].groupby("gene").size(),
        "n_distinct_switch_like_usable":
            det[det["anchored_usable"]].groupby("gene")["switch_like"].sum(),
    }).fillna(0).astype(int).reset_index()
    d = d.merge(per, on="gene", how="left")
    assert d["n_distinct_pairs"].notna().all()
    summ = json.loads((oc / "anchored_summary.json").read_text())["distinct_pairs"]
    assert d["n_distinct_detected"].sum() == summ["n_detected"]
    assert d["n_distinct_switch_like"].sum() == summ["n_switch_like"]
    assert d["n_distinct_switch_like_usable"].sum() == summ["n_switch_like_usable"]
    d = d.rename(columns={
        "n_anchored_pairs": "n_orientations",
        "n_pairs_detected": "n_orientations_detected",
        "n_switch_like": "n_orientations_switch_like",
        "n_pairs_usable": "n_orientations_usable",
        "n_switch_like_usable": "n_orientations_switch_like_usable",
    })
    lead = ["gene", "gene_name", "traits", "max_PP4_sQTL"] + list(per.columns[1:])
    d = d[lead + [c for c in d.columns if c not in lead]]
    d = d.sort_values(["confirmed_at_usable_abundance", "n_distinct_switch_like"],
                      ascending=[False, False]).reset_index(drop=True)
    write(round_num(d), "tableS29_longread_coloc_confirmation.csv")


def smr_table() -> None:
    """S30 — SMR/HEIDI on the colocalization nominations (GTEx QTLs).

    One row per probe in the primary confirmatory family (one pre-designated probe per
    gene and QTL class), with its coloc PP4, SMR estimate and P, HEIDI P, instrument
    strength, status and agreement with coloc. The secondary event-localization family
    (a gene's other introns) is corrected apart and left in the analysis ledger.
    `no_instrument` means untested, not negative.
    """
    s = pd.read_parquet(stage_out("anchoring", "smr_heidi", "gtex", "smr_results.parquet"))
    s = s[s.probe_family == "primary"]
    cols = ["analysis", "trait", "symbol", "gene", "modality", "probeID", "tissue",
            "topSNP", "coloc_PP4_probe", "coloc_estimator_probe", "b_SMR", "se_SMR",
            "p_SMR", "smr_threshold", "p_HEIDI", "nsnp_HEIDI", "F_instrument",
            "weak_instrument", "smr_status", "agreement"]
    s = s[cols].sort_values(["modality", "trait", "p_SMR"]).reset_index(drop=True)
    write(round_num(s, sig_cols=("b_SMR", "se_SMR", "p_SMR", "p_HEIDI")),
          "tableS30_smr_heidi.csv")


def brainseq_coloc_table() -> None:
    """S31 — BrainSEQ switch- (S_g) versus abundance-axis (A_g) colocalization.

    Same-cohort genetic anchoring, not replication (BrainSEQ is the discovery cohort).
    One row per analysis plus a pooled row: genes testable on both axes, calls per axis,
    the discordant counts and the exact McNemar P, and paired Wilcoxon tests of the
    conditional and raw PP4.
    """
    c = pd.read_parquet(stage_out("anchoring", "coloc_brainseq", "ea_only",
                                  "contrast.parquet"))
    write(round_num(c, sig_cols=("mcnemar_p", "wilcoxon_cond_p", "wilcoxon_pp4_p")),
          "tableS31_brainseq_axis_coloc_contrast.csv")


def rbp_tables() -> None:
    """S32 / S33 — candidate RBP regulons of the switch modules.

    S32: every module x RBP cell (mature-transcript scope, per RBP) with the
    hypergeometric over-representation test and the opportunity-adjusted binomial GLM
    (length, GC, UTR/CDS composition, transcript count). `supported` = hypergeometric
    q < 0.05 and adjusted q < 0.05 with OR > 1. Member-gene lists are left in the ledger.
    S33: ENCODE eCLIP (HepG2/K562) binding at switched versus constitutive exon intervals,
    one exact McNemar per RBP over its unique nominated genes, BH across testable RBPs.
    """
    r = pd.read_parquet(RBP / "rbp_regulon.parquet")
    r["supported"] = (r.q < 0.05) & (r.qvalue_glm < 0.05) & (r.odds_ratio > 1)
    assert len(r) == 11040 and (r.q < 0.05).sum() == 384 and r.supported.sum() == 89
    r = r.drop(columns=["module_genes", "universe_genes"], errors="ignore")
    r = r.sort_values(["supported", "q"], ascending=[False, True]).reset_index(drop=True)
    write(round_num(r, sig_cols=("p", "q", "pvalue_glm", "qvalue_glm")),
          "tableS32_rbp_regulons.csv")
    b = pd.read_parquet(RBP / "rbp_binding_support.parquet")
    assert len(b) == 35 and b.binding_supported.sum() == 21
    b = b.sort_values("mcnemar_fdr").reset_index(drop=True)
    write(round_num(b, sig_cols=("mcnemar_p", "mcnemar_fdr")), "tableS33_rbp_eclip_binding.csv")


def ldsc_table() -> None:
    """S27 — partitioned heritability of the switch-derived QTL annotations.

    `coef_p` is the one-sided test that the annotation's per-SNP coefficient exceeds 0,
    conditional on baselineLD (single-annotation models) or on baselineLD plus the other
    QTL layer (`joint`). `coef_q_bh` is its Benjamini-Hochberg q over the one testing
    family: the 15 single-annotation models of the aging switch layer (5 traits x
    cis/sQTL/eQTL). Joint and SCZD-module rows stay descriptive and have no q. The text
    and Fig 5e quote the q.
    """
    ld = pd.read_parquet(LDSC / "ldsc_partitioned.parquet")
    order = {"cis_only": 0, "sqtl_only": 1, "eqtl_only": 2, "joint": 3}
    ld = (ld.assign(_o=ld.model.map(order))
          .sort_values(["annotation", "trait", "_o", "annot"]).drop(columns="_o")
          .reset_index(drop=True))
    fam = (ld.annotation == "aging") & ld.model.isin(["cis_only", "sqtl_only", "eqtl_only"])
    if fam.sum() != 15:
        raise ValueError(f"expected 15 single-annotation aging tests, found {fam.sum()}")
    p = ld.loc[fam, "coef_p"].sort_values(ascending=False)
    rank = pd.Series(range(len(p), 0, -1), index=p.index)
    ld["coef_q_bh"] = (p * len(p) / rank).cummin().clip(upper=1).reindex(ld.index)
    write(round_num(ld, sig_cols=("enrichment_p", "coef", "coef_se", "coef_p",
                                  "coef_q_bh")), "tableS27_ldsc_partitioned.csv")


# --- S34  cohort description (Methods) ---------------------------------------
# One row per discovery analysis, read from the bundle sample tables and manifests that
# every IsoGraph fit consumed. GTEx exact age is released top-coded at 70.
BUNDLES = rel("inputs", "bundles")
BUNDLE_ANALYSES = (  # (bundle suite, region, cohort, analysis)
    ("brainseq_v1", "caudate", "BrainSEQ", "Aging"),
    ("brainseq_v1", "dlpfc", "BrainSEQ", "Aging"),
    ("brainseq_v1", "hippocampus", "BrainSEQ", "Aging"),
    ("brainseq_sczd", "caudate", "BrainSEQ", "Schizophrenia"),
)
GTEX_RACE = {1: "Asian", 2: "Black", 3: "White", 4: "American Indian", 98: "Unknown",
             99: "Unknown"}
BSEQ_RACE = {"AA": "Black", "CAUC": "White", "HISP": "Hispanic", "AS": "Asian",
             "Multi-Racial": "Multiracial"}
HARDY = {0: "ventilator", 1: "violent/fast", 2: "fast natural", 3: "intermediate",
         4: "slow"}
MUSIC_REF = {  # Tran/LIBD reference per analysis (celltype_composition.GTEX_REF + BrainSEQ)
    ("BrainSEQ", "caudate"): "Caudate (direct)", ("BrainSEQ", "dlpfc"): "DLPFC (direct)",
    ("BrainSEQ", "hippocampus"): "Hippocampus (direct)",
    ("GTEx", "amygdala"): "Amygdala (direct)",
    ("GTEx", "anterior_cingulate_cortex_ba24"): "sACC (direct)",
    ("GTEx", "frontal_cortex_ba9"): "DLPFC (cortical proxy)",
    ("GTEx", "cortex"): "DLPFC (cortical proxy)",
    ("GTEx", "hippocampus"): "Hippocampus (direct)",
    ("GTEx", "caudate_basal_ganglia"): "NAc (striatal proxy)",
    ("GTEx", "putamen_basal_ganglia"): "NAc (striatal proxy)",
    ("GTEx", "nucleus_accumbens_basal_ganglia"): "NAc (direct)",
}


def _counts(s: pd.Series, labels: dict) -> str:
    vc = s.map(lambda v: labels.get(v, "Unknown") if pd.notna(v) else "Unknown").value_counts()
    return "; ".join(f"{k} {v}" for k, v in vc.items())


def _med_range(s: pd.Series, n: int = 1) -> str:
    s = s.dropna()
    return f"{s.median():.{n}f} [{s.min():.{n}f}-{s.max():.{n}f}]"


def _med_iqr(s: pd.Series, n: int = 1) -> str:
    s = s.dropna()
    return f"{s.median():.{n}f} [{s.quantile(.25):.{n}f}-{s.quantile(.75):.{n}f}]"


def cohort_table() -> None:
    import json
    specs = list(BUNDLE_ANALYSES) + [
        ("gtex_v11_brain", p.name, "GTEx", "Aging")
        for p in sorted((BUNDLES / "gtex_v11_brain").iterdir()) if p.is_dir()]
    rows = []
    for suite, region, cohort, analysis in specs:
        d = BUNDLES / suite / region
        s = pd.read_parquet(d / "samples.parquet")
        prov = json.loads((d / "manifest.json").read_text())["provenance"]
        bseq = cohort == "BrainSEQ"
        donor = s["BrNum"] if bseq else s["SUBJID"]
        female = (s["Sex"] == "F") if bseq else (s["SEX"] == 2)
        mod = {m: m for m in s["MoD"].dropna().unique()} if bseq else HARDY
        rows.append({
            "cohort": cohort, "region": region, "analysis": analysis,
            "library": (
                "; ".join(f"{k} {v}" for k, v in s["Protocol"].value_counts().items())
                if bseq else "poly(A) (GTEx v11)"),
            "n_samples": len(s), "n_donors": donor.nunique(),
            "n_control": int((s["Dx"] == "Control").sum()) if bseq else len(s),
            "n_schizophrenia": int((s["Dx"] == "SCZD").sum()) if bseq else 0,
            "n_female": int(female.sum()), "pct_female": round(100 * female.mean(), 1),
            "age_median_range": _med_range(s["Age" if bseq else "AGE"]),
            "ancestry": _counts(s["Race" if bseq else "RACE"], BSEQ_RACE if bseq else GTEX_RACE),
            "rin_median_iqr": _med_iqr(s["RIN" if bseq else "SMRIN"]),
            "pmi_h_median_iqr": _med_iqr(s["PMI"]) if bseq else "",
            "ischemic_time_min_median_iqr": "" if bseq else _med_iqr(s["SMTSISCH"], 0),
            "death_classification": (("Manner " if bseq else "Hardy ")
                                     + _counts(s["MoD" if bseq else "DTHHRDY"], mod)),
            "genes_before_filter": int(prov["n_genes_before_expression_filter"]),
            "genes_after_filter": int(prov["n_genes_after_expression_filter"]),
            "transcripts_before_filter": int(prov["n_transcripts_before_expression_filter"]),
            "transcripts_after_filter": int(prov["n_transcripts_after_expression_filter"]),
            "expression_filter": prov["expression_filter"],
            "music_reference": MUSIC_REF.get((cohort, region), "none (no matched reference)"),
        })
    df = pd.DataFrame(rows)
    assert len(df) == 17 and (df.n_samples == df.n_donors).all(), "one sample per donor"
    write(df, "tableS34_cohort_description.csv")


# --- S35  synthetic scenario grid (Methods) --------------------------------
# Seed rule mirrors isograph_benchmark/benchmark/run_synthetic.py.
SCALE_SCENARIOS = {"scale", "scale_realistic"}


def synthetic_grid_table() -> None:
    import yaml
    cfg = yaml.safe_load(rel("configs", "synthetic_grid.yaml").read_text())
    rows = []
    for name, sc in cfg["scenarios"].items():
        grid = {k: v for k, v in sc.items() if k != "seed_count"}
        n_cells = 1
        for v in grid.values():
            n_cells *= len(v)
        seeds = sc.get("seed_count", cfg["seed_count_scale"] if name in SCALE_SCENARIOS
                       else cfg["seed_count_core"])
        other = {k: v for k, v in grid.items() if k not in ("n_genes", "n_samples")}
        rows.append({
            "scenario": name,
            "n_genes": "; ".join(map(str, sc["n_genes"])),
            "n_samples": "; ".join(map(str, sc["n_samples"])),
            "swept_parameters": "; ".join(f"{k} = {', '.join(map(str, v))}"
                                          for k, v in other.items() if len(v) > 1),
            "fixed_parameters": "; ".join(f"{k} = {v[0]}"
                                          for k, v in other.items() if len(v) == 1),
            "grid_cells": n_cells, "seeds_per_cell": seeds, "datasets": n_cells * seeds,
        })
    write(pd.DataFrame(rows), "tableS35_synthetic_scenarios.csv")


# --- S36  MAGMA competitive module gene-set tests -----------------------------
def magma_table() -> None:
    """S36 -- every MAGMA module x trait competitive test, both IsoGraph resolutions and the
    gene-level WGCNA baseline (manuscript Table S33; backs Fig. S24).

    The load-bearing column is `giant_module` (>= 900 genes): the manuscript reports that
    most FDR-significant sets are giant in every pipeline, so the table must carry the size
    beside the P value rather than only the hits.
    """
    gwas = stage_out("anchoring", "gwas")
    prod, r5 = gwas / "magma_results_combined.parquet", gwas / "magma_results_combined_res5.parquet"
    if not (prod.exists() and r5.exists()):
        print("  skip tableS36: MAGMA combined outputs absent")
        return
    p = pd.read_parquet(prod).assign(resolution="2.0")
    s = pd.read_parquet(r5).assign(resolution="5.0")
    # The resolution-5.0 file re-lists the resolution-independent WGCNA rows; keep them once.
    s = s[s["backend"] != "wgcna_gene"]
    d = pd.concat([p, s], ignore_index=True)
    d.loc[d["backend"] == "wgcna_gene", "resolution"] = "n/a"
    parts = d["VARIABLE"].str.split("__", n=2, expand=True)
    d["cohort"], d["region"], d["module_id"] = parts[0], parts[1], parts[2]
    d["giant_module"] = d["NGENES"] >= 900
    d = d.rename(columns={"NGENES": "n_genes", "BETA": "beta", "BETA_STD": "beta_std",
                          "SE": "se", "P": "p_value", "FDR": "fdr"})
    d = d[["backend", "resolution", "trait", "cohort", "region", "module_id", "n_genes",
           "giant_module", "beta", "beta_std", "se", "p_value", "fdr"]]
    d = d.sort_values(["backend", "resolution", "trait", "fdr"]).reset_index(drop=True)
    write(round_num(d, sig_cols=("fdr",)), "tableS36_magma_module_gwas.csv")


# --- S37  gnomAD constraint and ClinVar density of switched exons ------------
def clinical_table() -> None:
    """S37 -- cross-analysis rollup of switch-gene LOEUF and switched-exon ClinVar density
    (manuscript Table S34; backs Fig. S25a,b). One row per (exon scope, module stratum).
    """
    f = MECH / "clinical_consequence_meta.parquet"
    if not f.exists():
        print("  skip tableS37: clinical_consequence_meta.parquet absent")
        return
    d = pd.read_parquet(f)
    d = d.rename(columns={"n_regions": "n_analyses",
                          "n_ratio_gt1_p05": "n_analyses_enriched_p05",
                          "n_ratio_lt1_p05": "n_analyses_depleted_p05",
                          "median_ratio": "median_density_ratio",
                          "fisher_p": "density_fisher_p"})
    write(round_num(d, sig_cols=("density_fisher_p", "loeuf_fisher_p")),
          "tableS37_clinical_consequence.csv")


# --- S38 / S39  short-read PSI corroboration ---------------------------------
def psi_tables() -> None:
    """S38 -- gene-level PSI corroboration of the switch layer, one row per analysis and
    switch set (manuscript Table S35); S39 -- per-event short-read junction confirmation of
    the CLPP-anchored switch pairs (manuscript Table S36).

    S38 is rolled up from the committed `summary_<set>.json` files, which are the source of
    record for `validate_switch_splicing.py`; the older SWITCH_VALIDATION_SUMMARY.md predates
    the switching-filter re-run and is not read.
    """
    import json

    sv = MECH / "switch_validation"
    rows = []
    for f in sorted(sv.glob("*/summary_*.json")):
        j = json.loads(f.read_text())
        g = j["gene_corroboration"]
        e = j.get("event_resolved") or {}
        rows.append({
            "cohort": j["cohort"], "region": j["region"], "trait": j["trait"],
            "switch_set": j["switch_set"], "n_genes_tested": g["n_genes"],
            "n_switch_positive": g["n_switch_pos"], "n_psi_positive": g["n_psi_pos"],
            "n_both": g["contingency"]["switch_psi"],
            "psi_positive_rate_switch": g["psi_pos_rate_in_switch_genes"],
            "psi_positive_rate_background": g["psi_pos_rate_in_background"],
            "fisher_or": g["fisher_or"], "fisher_p": g["fisher_p"],
            "adjusted_or": g["adjusted_or"], "adjusted_p": g["adjusted_p"],
            "n_resolved_junctions": e.get("n_resolved_junctions"),
            "n_matched_to_psi": e.get("n_matched_to_psi"),
            "n_both_significant": e.get("n_both_significant"),
            "n_concordant": e.get("n_concordant"),
        })
    if rows:
        d = pd.DataFrame(rows).sort_values(["switch_set", "cohort", "region", "trait"])
        write(round_num(d.reset_index(drop=True), sig_cols=("fisher_p",)),
              "tableS38_psi_gene_corroboration.csv")
    else:
        print("  skip tableS38: no switch_validation summaries")

    jc = MECH / "junction_coloc_confirm" / "junction_confirm.parquet"
    if not jc.exists():
        print("  skip tableS39: junction_confirm.parquet absent")
        return
    d = pd.read_parquet(jc)
    cols = [c for c in ("gene_name", "ens", "trait", "tissue", "region", "match",
                        "anchored_junction", "clpp", "risk_allele", "mode", "event_id",
                        "event_type", "is_reference_contrast", "n_samples", "median_psi",
                        "minor_form_usage", "psi_iqr", "frac_samples_minor_ge_thresh",
                        "verdict", "verdict_reason") if c in d.columns]
    d = d[cols].sort_values(["gene_name", "trait", "region", "event_id"], na_position="last")
    write(round_num(d.reset_index(drop=True)), "tableS39_psi_junction_confirmation.csv")


def main() -> None:
    print(f"Writing supplementary tables to {OUT}")
    cohort_table()
    synthetic_grid_table()
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
    allelic_tables()
    clpp_event_table()
    signal_coloc_tables()
    coloc_event_table()
    longread_coloc_table()
    smr_table()
    brainseq_coloc_table()
    rbp_tables()
    ldsc_table()
    magma_table()
    clinical_table()
    psi_tables()
    print("done.")


if __name__ == "__main__":
    main()
