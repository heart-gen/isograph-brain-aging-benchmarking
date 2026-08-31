"""Genetic anchoring of IsoGraph co-switch modules via GTEx brain xQTL.

Tests whether genes in IsoGraph co-switch modules are enriched for GTEx brain
*sQTL* (splicing) and, as a specificity contrast, *eQTL* (expression). Co-switch
modules are defined on isoform usage, so the on-thesis prediction is sQTL
enrichment that is stronger than (or present without) eQTL enrichment — genetic
evidence for a regulated DTU-without-DGE layer.

Enrichment is a power-matched gene-set test, not QTL re-discovery: within each
xQTL's tested-gene universe (sGenes/eGenes, one row per gene with qval) intersected
with IsoGraph's tested genes, a logistic model regresses sGene/eGene status on
module membership while adjusting for the standard QTL-detectability confounds
(cis-variant count, gene length, isoform/intron multiplicity). The module-membership
odds ratio is the matched enrichment.

Modules are defined in the discovery cohort (e.g. BrainSEQ for SCZD) and tested in
independent GTEx donors, so there is no discovery/validation circularity. Genes are
matched on unversioned Ensembl id. Reads only saved artifacts + the copied xQTL
catalog; deterministic.

`--method` selects the graph whose modules are anchored: `isograph` (default) or the
matched WGCNA baselines `wgcna_switch_only` / `wgcna_multiplex`. Because those baselines
consume the same switch/multiplex features as IsoGraph, anchoring them tests whether the
splicing-QTL specificity is a property of the switch features or of IsoGraph's inference.

READ THE TWO ARMS TOGETHER. Co-switch module genes are cis-QTL *depleted* for BOTH
modalities (OR < 1 in every tissue tested, both kinds), which is the expected baseline
for coordinated/network genes. The result is that splicing-QTL is spared RELATIVE to
expression-QTL — a ratio of two depletions, never an sQTL enrichment. Describing it as
enrichment misstates the finding.

Sensitivity arms (`--outcome` / `--covariate-set`) address the two standing objections
to that reading. The PRIMARY arm is `binary` + `standard` and keeps the bare filename,
so a sensitivity run can never overwrite the numbers the figures and tables use:

  --covariate-set constraint  Adds gnomAD v4.1 LOEUF, missense z and log mean expression
      in this cohort/region's own count matrix. Tests the leading alternative
      explanation — that the shared depletion, and hence the contrast, is selective
      constraint on network-central genes. LOEUF is undefined for ~20% of tested genes,
      so the universe is restricted to constraint-complete genes FIRST and BOTH covariate
      sets are then fitted on that identical subset; otherwise "survives adjustment"
      would be confounded with "survives subsetting".
  --outcome continuous  Rank-INT of -log10(pval_beta), the permutation statistic the
      sGene/eGene call thresholds away. Threshold-free, and adds precision without
      adding data. NOT a substitute for the primary: it sharpens the eQTL arm far more
      than the sQTL arm, so leading with it would improve the ratio for a reason that
      is really about the denominator.
  --outcome dose  Poisson on the number of independent SuSiE credible sets — how many
      distinct cis signals a gene carries, LD-resolved rather than variant-counted.

Deliberately NOT used: GTEx all-pairs nominal summary statistics. They reintroduce LD
and the multiple-variants-per-gene problem the permutation pass already solved, cost
TB across 13 tissues x 2 modalities, and collapse to a per-gene statistic weaker than
`pval_beta` (min-p is not calibrated for the number of variants tested).

Per analysis it writes, under <artifact-parent>/_m/ (non-isograph methods and every
sensitivity arm get a filename suffix so the primary IsoGraph outputs are untouched):
  qtl_anchoring[_<method>][_<outcome>][_<covariate_set>].parquet — one row per
      (xqtl_kind, module_set, covariate_set): graph_method, matched OR, CI, p-value,
      foreground size, qtl rate in foreground vs background, outcome, covariate_set.
  QTL_ANCHORING[...].md — the sQTL-vs-eQTL contrast writeup.
  qtl_anchoring[...].json — run parameters (tissue, fdr, universe sizes, covariate
      coverage for the constraint arm).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.interpret_modules import DEFAULT_GTF_CACHE
from isograph_benchmark.real_data.partition_provenance import load_enrichment
from isograph_benchmark.real_data.project_tiers import _bundle_path
from isograph_benchmark.real_data.sweep_leiden import _artifact_dir

DEFAULT_XQTL_DIR = rel("inputs", "raw", "gtex_v11", "xqtl")
DEFAULT_CONSTRAINT = rel("inputs", "raw", "clinical", "gnomad.v4.1.constraint_metrics.tsv")
_QVAL = 0.05

# graph method -> (artifact subdir under <_m>, module_enrichment table basename).
# The matched WGCNA baselines consume the same switch/multiplex feature matrix as
# IsoGraph, so anchoring them tests whether the splicing-QTL specificity is a property
# of the switch features or of IsoGraph's VAE + Leiden inference.
_METHODS = {
    "isograph": ("isograph_vae", "isograph_modules.parquet"),
    "wgcna_switch_only": ("wgcna_switch_only", "wgcna_switch_modules.parquet"),
    "wgcna_multiplex": ("wgcna_multiplex", "wgcna_multiplex_modules.parquet"),
}

# IsoGraph region naming -> GTEx v11 brain tissue file prefix.
_GTEX_TISSUE = {
    "amygdala": "Brain_Amygdala",
    "anterior_cingulate_cortex_ba24": "Brain_Anterior_cingulate_cortex_BA24",
    "caudate_basal_ganglia": "Brain_Caudate_basal_ganglia",
    "cerebellar_hemisphere": "Brain_Cerebellar_Hemisphere",
    "cerebellum": "Brain_Cerebellum",
    "cortex": "Brain_Cortex",
    "frontal_cortex_ba9": "Brain_Frontal_Cortex_BA9",
    "hippocampus": "Brain_Hippocampus",
    "hypothalamus": "Brain_Hypothalamus",
    "nucleus_accumbens_basal_ganglia": "Brain_Nucleus_accumbens_basal_ganglia",
    "putamen_basal_ganglia": "Brain_Putamen_basal_ganglia",
    "spinal_cord_cervical_c_1": "Brain_Spinal_cord_cervical_c-1",
    "substantia_nigra": "Brain_Substantia_nigra",
}
# BrainSEQ regions -> closest GTEx brain tissue.
_BRAINSEQ_TISSUE = {
    "caudate": "Brain_Caudate_basal_ganglia",
    "caudate_sczd": "Brain_Caudate_basal_ganglia",
    "hippocampus": "Brain_Hippocampus",
    "dlpfc": "Brain_Frontal_Cortex_BA9",
}


def resolve_tissue(analysis: str, region: str | None) -> str:
    if analysis == "brainseq-sczd":
        return _BRAINSEQ_TISSUE["caudate_sczd"]
    if analysis == "brainseq-aging":
        return _BRAINSEQ_TISSUE[region]
    if analysis == "gtex-aging":
        return _GTEX_TISSUE[region]
    raise ValueError(f"no tissue mapping for analysis={analysis!r} region={region!r}")


def _bare(series: pd.Series) -> pd.Series:
    return series.astype(str).str.split(".", n=1).str[0]


def load_qtl_genes(xqtl_dir: Path, tissue: str, kind: str, fdr: float) -> pd.DataFrame:
    """One row per tested gene: gene (bare), is_qtl, evidence, num_var, gene_len, group_size.

    ``is_qtl`` is the thresholded sGene/eGene call (the primary outcome). ``evidence``
    is the permutation-calibrated gene-level p-value on a -log10 scale, kept because
    thresholding at qval<=fdr discards most of the information the permutation pass
    already computed — see the ``continuous`` outcome.
    """
    suffix = "sGenes" if kind == "sQTL" else "eGenes"
    path = xqtl_dir / f"{tissue}.v11.{suffix}.txt.gz"
    d = pd.read_csv(path, sep="\t")
    d["gene"] = _bare(d["gene_id"])
    d["gene_len"] = (d["gene_end"] - d["gene_start"]).abs() + 1
    d["is_qtl"] = (d["qval"] <= fdr).astype(int)
    # pval_beta is the beta-approximated permutation p: already LD-aware and already
    # calibrated for the number of variants tested, which min-p over all pairs is not.
    d["evidence"] = -np.log10(pd.to_numeric(d["pval_beta"], errors="coerce").clip(lower=1e-300))
    cols = {"is_qtl": "max", "evidence": "max", "num_var": "max", "gene_len": "max"}
    if "group_size" in d.columns:
        cols["group_size"] = "max"
    return d.groupby("gene", as_index=False).agg(cols)


def load_credible_set_counts(xqtl_dir: Path, tissue: str, kind: str) -> pd.Series:
    """gene -> number of DISTINCT SuSiE credible sets (independent cis signals).

    The dose outcome. A binary sGene/eGene call cannot distinguish a gene with one
    weak signal from a gene with four independent ones; the credible-set count can,
    and it is LD-resolved rather than variant-counted. Genes tested but with no
    credible set are absent here and are filled with 0 by the caller.
    """
    path = xqtl_dir / f"{tissue}.v11.{'sQTLs' if kind == 'sQTL' else 'eQTLs'}.SuSiE_summary.parquet"
    if not path.exists():
        return pd.Series(dtype=float)
    # The two catalogs are shaped differently: sQTL phenotypes are intron clusters and
    # carry a separate gene_id, while for eQTL the phenotype IS the gene and no gene_id
    # column exists. Resolve that before grouping rather than assuming one schema.
    import pyarrow.parquet as pq

    names = set(pq.read_schema(path).names)
    gene_col = "gene_id" if "gene_id" in names else "phenotype_id"
    d = pd.read_parquet(path, columns=sorted({gene_col, "phenotype_id", "cs_id"}))
    d["gene"] = _bare(d[gene_col])
    # A credible set is identified by (phenotype, cs_id): one gene can carry several
    # splicing phenotypes, each with its own independent sets.
    return (d[["gene", "phenotype_id", "cs_id"]].drop_duplicates()
            .groupby("gene").size().astype(float))


def gene_constraint(constraint_path: Path) -> pd.DataFrame:
    """gene -> LOEUF and missense z from gnomAD v4.1, one row per gene.

    LOEUF (``lof.oe_ci.upper``) is the covariate that tests the leading alternative
    explanation for the anchoring result: co-switch module genes are cis-QTL DEPLETED
    for both modalities, which is what selective constraint on network-central genes
    would produce. If constraint drives it, adjusting for LOEUF should collapse the
    splicing-specificity contrast.
    """
    if not constraint_path.exists():
        return pd.DataFrame(columns=["gene", "loeuf", "mis_z"])
    d = pd.read_csv(constraint_path, sep="\t", low_memory=False,
                    usecols=lambda c: c in {"gene_id", "mane_select", "canonical",
                                            "lof.oe_ci.upper", "mis.z_score"})
    # One transcript per gene: MANE Select where it exists, else the canonical.
    for flag in ("mane_select", "canonical"):
        if flag in d.columns:
            pick = d[d[flag].astype(str).str.lower().isin({"true", "t", "1"})]
            if not pick.empty:
                d = pick
                break
    d["gene"] = _bare(d["gene_id"])
    out = (d.rename(columns={"lof.oe_ci.upper": "loeuf", "mis.z_score": "mis_z"})
            .groupby("gene", as_index=False)[["loeuf", "mis_z"]].median())
    return out


def gene_expression(bundle_path: Path) -> pd.Series:
    """gene -> mean log1p expression in this cohort/region's own count matrix.

    Tissue-matched, and the same matrix the modules were discovered on. Expression
    governs xQTL discovery power asymmetrically — eQTL power tracks expression more
    directly than sQTL power tracks intron coverage — so it does NOT cancel in the
    paired sQTL/eQTL contrast the way a symmetric covariate would.
    """
    counts_f, genes_f = bundle_path / "gene_counts.npz", bundle_path / "genes.parquet"
    if not (counts_f.exists() and genes_f.exists()):
        return pd.Series(dtype=float)
    mat = np.load(counts_f, allow_pickle=True)["data"]
    genes = pd.read_parquet(genes_f)["gene_id"]
    if mat.shape[0] != len(genes):
        return pd.Series(dtype=float)
    mean_log = np.log1p(mat).mean(axis=1)
    return pd.Series(mean_log, index=_bare(genes).values).groupby(level=0).max()


def isoform_counts(gtf_cache: Path) -> pd.Series:
    c = pd.read_parquet(gtf_cache, columns=["transcript_id", "gene_id", "feature"])
    tx = c[c["feature"] == "transcript"]
    return tx.groupby(_bare(tx["gene_id"]))["transcript_id"].nunique()


def build_gene_sets(mod_dir: Path, enrich_path: Path, fdr: float,
                    context: str | None = None) -> dict[str, set[str]]:
    modules = pd.read_parquet(mod_dir / "modules.parquet")
    modules["gene"] = _bare(modules["gene_id"])
    modules["module_id"] = modules["module_id"].astype(str)
    sets = {"all_modules": set(modules["gene"])}

    # Leiden module ids are re-assigned on every fit, so joining a stale enrichment
    # table here silently scrambles pheno_fdr / n_go_terms across modules rather
    # than failing. That is exactly how the 2026-06-29 GO-invisible partition went
    # wrong; load_enrichment turns it into a hard error.
    enrich = load_enrichment(
        enrich_path, modules, context=context or f"build_gene_sets {mod_dir}"
    )
    if enrich is not None:
        enrich["module_id"] = enrich["module_id"].astype(str)
        sig = enrich[enrich["pheno_fdr"] <= fdr]
        sig_ids = set(sig["module_id"])
        inv_ids = set(sig.loc[sig["n_go_terms"] == 0, "module_id"])
        vis_ids = set(sig.loc[sig["n_go_terms"] > 0, "module_id"])
        member = lambda ids: set(modules.loc[modules["module_id"].isin(ids), "gene"])
        sets["pheno_sig_modules"] = member(sig_ids)
        sets["go_invisible_modules"] = member(inv_ids)
        sets["go_visible_modules"] = member(vis_ids)
    return {k: v for k, v in sets.items() if v}


#: Outcome -> (column, model). All three report the module-membership effect on a LOG
#: scale, so the paired sQTL-vs-eQTL contrast (a difference of logs) is coherent within
#: an outcome. Effects are NOT comparable ACROSS outcomes and must not be pooled.
OUTCOMES = {
    "binary": "is_qtl",       # thresholded sGene/eGene call — logistic; the primary
    "continuous": "evidence",  # rank-INT of -log10(pval_beta) — OLS; threshold-free
    "dose": "n_credible_sets",  # independent SuSiE signals — Poisson; LD-resolved
}

#: Covariate sets. ``standard`` is the published QTL-detectability adjustment;
#: ``constraint`` adds selective constraint and tissue-matched expression, which
#: together test whether the contrast is a constraint/expression artefact.
COVARIATE_SETS = {
    "standard": ("log_num_var", "log_len", "log_iso", "log_group"),
    "constraint": ("log_num_var", "log_len", "log_iso", "log_group",
                   "loeuf", "mis_z", "log_expr"),
}


def _rank_int(x: pd.Series) -> pd.Series:
    """Rank-based inverse-normal transform (Blom). Makes the continuous outcome
    scale-free and robust to the heavy right tail of -log10 p."""
    from scipy.special import ndtri

    r = x.rank(method="average")
    return pd.Series(ndtri((r - 0.375) / (len(r) + 0.25)), index=x.index)


def matched_enrichment(universe: pd.DataFrame, in_set: set[str],
                       outcome: str = "binary",
                       covariate_set: str = "standard") -> dict:
    ycol = OUTCOMES[outcome]
    df = universe.copy()
    df["in_set"] = df["gene"].isin(in_set).astype(int)
    covars = [c for c in COVARIATE_SETS[covariate_set] if c in df.columns]
    df = df.dropna(subset=[ycol, "in_set", *covars])
    n_fg = int(df["in_set"].sum())
    # rate_fg / rate_bg stay on the BINARY call in every outcome: they are the raw
    # unmatched contrast the report quotes, and re-defining them per outcome would
    # make the summary tables silently incomparable.
    has_bin = "is_qtl" in df.columns
    rate_fg = float(df.loc[df["in_set"] == 1, "is_qtl"].mean()) if n_fg and has_bin else np.nan
    rate_bg = float(df.loc[df["in_set"] == 0, "is_qtl"].mean()) if has_bin else np.nan
    out = {"n_universe": int(len(df)), "n_foreground": n_fg,
           "n_fg_qtl": int(df.loc[df["in_set"] == 1, "is_qtl"].sum()) if has_bin else 0,
           "rate_fg": round(rate_fg, 4) if np.isfinite(rate_fg) else np.nan,
           "rate_bg": round(rate_bg, 4) if np.isfinite(rate_bg) else np.nan,
           "outcome": outcome, "covariate_set": covariate_set,
           "n_covariates": len(covars)}
    if n_fg < 5 or df["in_set"].nunique() < 2 or df[ycol].nunique() < 2:
        out.update(odds_ratio=np.nan, or_ci_low=np.nan, or_ci_high=np.nan,
                   pvalue=np.nan, method="insufficient")
        return out

    X = sm.add_constant(df[["in_set", *covars]].astype(float))
    y = df[ycol].astype(float)
    try:
        if outcome == "binary":
            res = sm.Logit(y, X).fit(disp=0, maxiter=200)
            label = "logit_matched"
        elif outcome == "continuous":
            res = sm.OLS(_rank_int(y), X).fit()
            label = "ols_rankint_matched"
        else:
            res = sm.GLM(y, X, family=sm.families.Poisson()).fit()
            label = "poisson_matched"
        ci = res.conf_int().loc["in_set"]
        # Exponentiated effect: an OR (logit), a rate ratio (Poisson), or exp of a
        # standardized mean difference (OLS). Kept in the same column so the meta
        # machinery is shared; `outcome` records which it is.
        out.update(odds_ratio=float(np.exp(res.params["in_set"])),
                   or_ci_low=float(np.exp(ci[0])), or_ci_high=float(np.exp(ci[1])),
                   pvalue=float(res.pvalues["in_set"]),
                   method=f"{label}_{covariate_set}")
    except Exception as exc:  # separation / non-convergence -> unmatched Fisher
        if outcome != "binary":
            out.update(odds_ratio=np.nan, or_ci_low=np.nan, or_ci_high=np.nan,
                       pvalue=np.nan, method=f"failed ({type(exc).__name__})")
            return out
        from scipy.stats import fisher_exact
        tab = pd.crosstab(df["in_set"], df["is_qtl"])
        orr, p = fisher_exact(tab.reindex(index=[1, 0], columns=[1, 0]).values)
        out.update(odds_ratio=float(orr), or_ci_low=np.nan, or_ci_high=np.nan,
                   pvalue=float(p), method=f"fisher_unmatched ({type(exc).__name__})")
    return out


def variant_suffix(method: str, outcome: str, covariate_set: str) -> str:
    """Filename suffix. The PRIMARY analysis (isograph, binary, standard covariates)
    keeps the bare filename, so sensitivity runs can never overwrite the numbers the
    figures and tables are built from."""
    parts = [] if method == "isograph" else [method]
    if outcome != "binary":
        parts.append(outcome)
    if covariate_set != "standard":
        parts.append(covariate_set)
    return "".join(f"_{p}" for p in parts)


def run_anchoring(analysis: str, region: str | None, variant: str, fdr: float,
                  xqtl_dir: Path, tissue: str | None, gtf_cache: Path,
                  method: str = "isograph", outcome: str = "binary",
                  covariate_set: str = "standard",
                  constraint_path: Path | None = None) -> pd.DataFrame:
    constraint_path = constraint_path or DEFAULT_CONSTRAINT
    iso_dir = _artifact_dir(analysis, region, variant)
    subdir, enrich_name = _METHODS[method]
    mod_dir = iso_dir.parent / subdir
    tissue = tissue or resolve_tissue(analysis, region)
    # universe = the method's own tested-gene set (node_diagnostics for IsoGraph; the
    # clustered genes for WGCNA, which has no separate node table).
    nd = mod_dir / "node_diagnostics.parquet"
    universe_src = nd if nd.exists() else mod_dir / "modules.parquet"
    iso_genes = set(_bare(pd.read_parquet(universe_src, columns=["gene_id"])["gene_id"]))
    iso_per_gene = isoform_counts(gtf_cache)
    enrich_path = iso_dir.parent / "module_enrichment" / enrich_name
    gene_sets = build_gene_sets(
        mod_dir, enrich_path, fdr,
        context=f"qtl_anchoring {analysis}/{region or ''} [{method}]")

    # Constraint / expression covariates are loaded once and reused across both xQTL
    # arms so the two models are adjusted identically — an asymmetric adjustment would
    # itself manufacture a splicing-specificity contrast.
    constraint = (gene_constraint(constraint_path) if covariate_set == "constraint"
                  else pd.DataFrame(columns=["gene", "loeuf", "mis_z"]))
    expr = (gene_expression(_bundle_path(analysis, region))
            if covariate_set == "constraint" else pd.Series(dtype=float))

    rows = []
    universe_sizes = {}
    covariate_coverage = {}
    for kind in ("sQTL", "eQTL"):
        qtl = load_qtl_genes(xqtl_dir, tissue, kind, fdr)
        qtl = qtl[qtl["gene"].isin(iso_genes)].copy()
        qtl["log_num_var"] = np.log1p(qtl["num_var"])
        qtl["log_len"] = np.log(qtl["gene_len"])
        qtl["log_iso"] = np.log(qtl["gene"].map(iso_per_gene).fillna(1.0).clip(lower=1))
        if "group_size" in qtl.columns:
            qtl["log_group"] = np.log(qtl["group_size"].clip(lower=1))
        if outcome == "dose":
            cs = load_credible_set_counts(xqtl_dir, tissue, kind)
            # Tested but with no credible set is a real 0, not missing data.
            qtl["n_credible_sets"] = qtl["gene"].map(cs).fillna(0.0)
        fit_sets = [covariate_set]
        if covariate_set == "constraint":
            qtl = qtl.merge(constraint, on="gene", how="left")
            qtl["log_expr"] = qtl["gene"].map(expr)
            n_before = len(qtl)
            # LOEUF is undefined for ~20% of tested genes (non-coding, and genes with
            # too few expected LoF). Dropping them inside the model would confound
            # "the contrast survives constraint adjustment" with "the contrast
            # survives restricting to LOEUF-covered genes". Restrict FIRST, then fit
            # BOTH covariate sets on the identical subset, so the two rows differ only
            # by the covariates and the comparison is a genuine nested-model test.
            qtl = qtl.dropna(subset=["loeuf", "mis_z", "log_expr"]).copy()
            covariate_coverage[kind] = {
                "n_tested": int(n_before),
                "n_constraint_complete": int(len(qtl)),
                "frac_retained": round(len(qtl) / n_before, 4) if n_before else np.nan,
            }
            fit_sets = ["standard", "constraint"]
        universe_sizes[kind] = int(len(qtl))
        for cset in fit_sets:
            for set_name, genes in gene_sets.items():
                res = matched_enrichment(qtl, genes, outcome=outcome, covariate_set=cset)
                rows.append({"analysis": analysis, "region": region or "", "tissue": tissue,
                             "graph_method": method, "xqtl_kind": kind,
                             "module_set": set_name, "fdr": fdr, **res})
    gate = pd.DataFrame(rows)

    out_dir = ensure_dir(iso_dir.parent)
    suffix = variant_suffix(method, outcome, covariate_set)
    gate.to_parquet(out_dir / f"qtl_anchoring{suffix}.parquet", index=False, compression="zstd")
    (out_dir / f"qtl_anchoring{suffix}.json").write_text(json.dumps(
        {"analysis": analysis, "region": region, "tissue": tissue, "variant": variant,
         "method": method, "outcome": outcome, "covariate_set": covariate_set,
         "fdr": fdr, "n_universe_genes": len(iso_genes),
         "universe_sizes": universe_sizes,
         "covariate_coverage": covariate_coverage,
         "module_set_sizes": {k: len(v) for k, v in gene_sets.items()}}, indent=2))
    _write_report(out_dir, analysis, region, tissue, gate, suffix, method,
                  outcome, covariate_set)
    return gate


_MODEL_BLURB = {
    "binary": ("Power-matched enrichment (logistic: sGene/eGene status ~ module "
               "membership + covariates) within each xQTL's tested-gene universe "
               "intersected with the method's tested genes."),
    "continuous": ("Threshold-free sensitivity (OLS: rank-inverse-normal of "
                   "-log10 permutation p ~ module membership + covariates). Uses the "
                   "same gene-level permutation statistic the sGene/eGene call "
                   "thresholds, so it adds precision without adding data."),
    "dose": ("Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ "
             "module membership + covariates). Distinguishes a gene with one weak cis "
             "signal from one with several independent ones; LD-resolved."),
}


def _markdown_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    rows = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    rows += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(rows)


def _write_report(out_dir: Path, analysis: str, region: str | None, tissue: str,
                  gate: pd.DataFrame, suffix: str = "", method: str = "isograph",
                  outcome: str = "binary", covariate_set: str = "standard") -> None:
    show = gate.copy()
    for c in ("odds_ratio", "or_ci_low", "or_ci_high"):
        show[c] = show[c].round(2)
    show["pvalue"] = show["pvalue"].apply(lambda p: f"{p:.2e}" if pd.notna(p) else "NA")
    show = show.rename(columns={"method": "fit_method"})
    table = show[["xqtl_kind", "module_set", "n_foreground", "rate_fg", "rate_bg",
                  "odds_ratio", "or_ci_low", "or_ci_high", "pvalue", "fit_method"]]
    lines = [
        f"# Genetic anchoring — {method} co-switch modules vs GTEx {tissue} xQTL "
        f"({analysis}" + (f"/{region}" if region else "") + ")",
        "",
        _MODEL_BLURB[outcome],
        "",
        ("Covariates: log cis-variant count, log gene length, log isoform count "
         "[+ log intron group size for sQTL]."
         if covariate_set == "standard" else
         "Covariates: log cis-variant count, log gene length, log isoform count "
         "[+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, "
         "and log mean expression in this cohort/region's own count matrix**. The "
         "constraint set tests whether the splicing-specificity contrast survives "
         "adjustment for selective constraint and expression level — the leading "
         "alternative explanation for module genes being cis-QTL depleted in BOTH "
         "modalities."),
        "",
        "Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring "
        f"--analysis {analysis}" + (f" --region {region}" if region else "")
        + (f" --method {method}" if method != "isograph" else "") + "`",
        "",
        "## Matched odds ratios",
        "",
        _markdown_table(table),
        "",
        "## Reading",
        "",
        "- **Read the two arms together, not the ratio alone.** Co-switch module genes "
        "are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); "
        "the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. "
        "sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, "
        "and must never be described as sQTL enrichment.",
        "- Matching on cis-variant count / gene length / isoform multiplicity controls "
        "the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw "
        "(unmatched) contrast for reference.",
        "- Scope: cis-sQTL enrichment shows module *members* undergo genetically "
        "regulated splicing; it does not by itself prove the *co-switching* is genetic "
        "(a shared trans regulator / cell composition could coordinate it).",
    ]
    (out_dir / f"QTL_ANCHORING{suffix}.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="GTEx sQTL/eQTL anchoring of co-switch modules.")
    p.add_argument("--analysis", required=True, help="brainseq-sczd, brainseq-aging, gtex-aging")
    p.add_argument("--region", default=None)
    p.add_argument("--variant", default="standard")
    p.add_argument("--method", default="isograph", choices=list(_METHODS),
                   help="graph method whose modules to anchor (default isograph)")
    p.add_argument("--fdr", type=float, default=_QVAL)
    p.add_argument("--xqtl-dir", default=str(DEFAULT_XQTL_DIR))
    p.add_argument("--tissue", default=None, help="override GTEx tissue prefix")
    p.add_argument("--gtf-cache", default=str(DEFAULT_GTF_CACHE))
    p.add_argument("--outcome", default="binary", choices=list(OUTCOMES),
                   help="binary = thresholded sGene/eGene (PRIMARY); continuous = "
                        "rank-INT of -log10 permutation p; dose = SuSiE credible-set "
                        "count. Sensitivity outcomes write to their own files.")
    p.add_argument("--covariate-set", default="standard", choices=list(COVARIATE_SETS),
                   help="standard = published QTL-detectability covariates (PRIMARY); "
                        "constraint = adds gnomAD LOEUF, missense z and log expression.")
    p.add_argument("--constraint-metrics", default=str(DEFAULT_CONSTRAINT))
    args = p.parse_args()
    gate = run_anchoring(args.analysis, args.region, args.variant, args.fdr,
                         Path(args.xqtl_dir), args.tissue, Path(args.gtf_cache),
                         method=args.method, outcome=args.outcome,
                         covariate_set=args.covariate_set,
                         constraint_path=Path(args.constraint_metrics))
    # --covariate-set constraint emits BOTH covariate sets on the same subset, so the
    # summary is keyed on (module_set, covariate_set), not module_set alone.
    piv = gate.pivot_table(index=["module_set", "covariate_set"], columns="xqtl_kind",
                           values="odds_ratio")
    print(f"[{args.method} / {args.outcome}] module_set      covariates   "
          f"sQTL     eQTL   sQTL/eQTL")
    for (mset, cset), r in piv.iterrows():
        s_or, e_or = r.get("sQTL", np.nan), r.get("eQTL", np.nan)
        ratio = s_or / e_or if np.isfinite(s_or) and np.isfinite(e_or) and e_or else np.nan
        print(f"{mset:22}{cset:12}{s_or:7.2f}  {e_or:7.2f}     {ratio:7.3f}")


if __name__ == "__main__":
    main()
