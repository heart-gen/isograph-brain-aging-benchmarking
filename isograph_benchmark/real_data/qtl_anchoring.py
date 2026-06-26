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

Per analysis it writes, under <artifact-parent>/_m/:
  qtl_anchoring.parquet — one row per (xqtl_kind, module_set): matched OR, CI,
      p-value, foreground size, qtl rate in foreground vs background, method.
  QTL_ANCHORING.md — the sQTL-vs-eQTL contrast writeup.
  qtl_anchoring.json — run parameters (tissue, fdr, universe sizes).
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
from isograph_benchmark.real_data.sweep_leiden import _artifact_dir

DEFAULT_XQTL_DIR = rel("inputs", "raw", "gtex_v11", "xqtl")
_QVAL = 0.05

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
    """One row per tested gene: gene (bare), is_qtl, num_var, gene_len, group_size."""
    suffix = "sGenes" if kind == "sQTL" else "eGenes"
    path = xqtl_dir / f"{tissue}.v11.{suffix}.txt.gz"
    d = pd.read_csv(path, sep="\t")
    d["gene"] = _bare(d["gene_id"])
    d["gene_len"] = (d["gene_end"] - d["gene_start"]).abs() + 1
    d["is_qtl"] = (d["qval"] <= fdr).astype(int)
    cols = {"gene": "first", "is_qtl": "max", "num_var": "max", "gene_len": "max"}
    if "group_size" in d.columns:
        cols["group_size"] = "max"
    agg = d.groupby("gene", as_index=False).agg({k: v for k, v in cols.items() if k != "gene"})
    return agg


def isoform_counts(gtf_cache: Path) -> pd.Series:
    c = pd.read_parquet(gtf_cache, columns=["transcript_id", "gene_id", "feature"])
    tx = c[c["feature"] == "transcript"]
    return tx.groupby(_bare(tx["gene_id"]))["transcript_id"].nunique()


def build_gene_sets(iso_dir: Path, fdr: float) -> dict[str, set[str]]:
    modules = pd.read_parquet(iso_dir / "modules.parquet")
    modules["gene"] = _bare(modules["gene_id"])
    modules["module_id"] = modules["module_id"].astype(str)
    sets = {"all_modules": set(modules["gene"])}

    enrich_path = iso_dir.parent / "module_enrichment" / "isograph_modules.parquet"
    if enrich_path.exists():
        enrich = pd.read_parquet(enrich_path)
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


def matched_enrichment(universe: pd.DataFrame, in_set: set[str]) -> dict:
    df = universe.copy()
    df["in_set"] = df["gene"].isin(in_set).astype(int)
    covars = ["log_num_var", "log_len", "log_iso"]
    if "log_group" in df.columns:
        covars.append("log_group")
    df = df.dropna(subset=["is_qtl", "in_set", *covars])
    n_fg = int(df["in_set"].sum())
    rate_fg = float(df.loc[df["in_set"] == 1, "is_qtl"].mean()) if n_fg else np.nan
    rate_bg = float(df.loc[df["in_set"] == 0, "is_qtl"].mean())
    out = {"n_universe": int(len(df)), "n_foreground": n_fg,
           "n_fg_qtl": int(df.loc[df["in_set"] == 1, "is_qtl"].sum()),
           "rate_fg": round(rate_fg, 4), "rate_bg": round(rate_bg, 4)}
    if n_fg < 5 or df["in_set"].nunique() < 2 or df["is_qtl"].nunique() < 2:
        out.update(odds_ratio=np.nan, or_ci_low=np.nan, or_ci_high=np.nan,
                   pvalue=np.nan, method="insufficient")
        return out
    X = sm.add_constant(df[["in_set", *covars]].astype(float))
    try:
        res = sm.Logit(df["is_qtl"].astype(float), X).fit(disp=0, maxiter=200)
        ci = res.conf_int().loc["in_set"]
        out.update(odds_ratio=float(np.exp(res.params["in_set"])),
                   or_ci_low=float(np.exp(ci[0])), or_ci_high=float(np.exp(ci[1])),
                   pvalue=float(res.pvalues["in_set"]), method="logit_matched")
    except Exception as exc:  # separation / non-convergence -> unmatched Fisher
        from scipy.stats import fisher_exact
        tab = pd.crosstab(df["in_set"], df["is_qtl"])
        orr, p = fisher_exact(tab.reindex(index=[1, 0], columns=[1, 0]).values)
        out.update(odds_ratio=float(orr), or_ci_low=np.nan, or_ci_high=np.nan,
                   pvalue=float(p), method=f"fisher_unmatched ({type(exc).__name__})")
    return out


def run_anchoring(analysis: str, region: str | None, variant: str, fdr: float,
                  xqtl_dir: Path, tissue: str | None, gtf_cache: Path) -> pd.DataFrame:
    iso_dir = _artifact_dir(analysis, region, variant)
    tissue = tissue or resolve_tissue(analysis, region)
    iso_genes = set(_bare(pd.read_parquet(iso_dir / "node_diagnostics.parquet")["gene_id"]))
    iso_per_gene = isoform_counts(gtf_cache)
    gene_sets = build_gene_sets(iso_dir, fdr)

    rows = []
    universe_sizes = {}
    for kind in ("sQTL", "eQTL"):
        qtl = load_qtl_genes(xqtl_dir, tissue, kind, fdr)
        qtl = qtl[qtl["gene"].isin(iso_genes)].copy()
        qtl["log_num_var"] = np.log1p(qtl["num_var"])
        qtl["log_len"] = np.log(qtl["gene_len"])
        qtl["log_iso"] = np.log(qtl["gene"].map(iso_per_gene).fillna(1.0).clip(lower=1))
        if "group_size" in qtl.columns:
            qtl["log_group"] = np.log(qtl["group_size"].clip(lower=1))
        universe_sizes[kind] = int(len(qtl))
        for set_name, genes in gene_sets.items():
            res = matched_enrichment(qtl, genes)
            rows.append({"analysis": analysis, "region": region or "", "tissue": tissue,
                         "xqtl_kind": kind, "module_set": set_name, "fdr": fdr, **res})
    gate = pd.DataFrame(rows)

    out_dir = ensure_dir(iso_dir.parent)
    gate.to_parquet(out_dir / "qtl_anchoring.parquet", index=False, compression="zstd")
    (out_dir / "qtl_anchoring.json").write_text(json.dumps(
        {"analysis": analysis, "region": region, "tissue": tissue, "variant": variant,
         "fdr": fdr, "n_isograph_genes": len(iso_genes),
         "universe_sizes": universe_sizes,
         "module_set_sizes": {k: len(v) for k, v in gene_sets.items()}}, indent=2))
    _write_report(out_dir, analysis, region, tissue, gate)
    return gate


def _markdown_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    rows = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    rows += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(rows)


def _write_report(out_dir: Path, analysis: str, region: str | None, tissue: str,
                  gate: pd.DataFrame) -> None:
    show = gate.copy()
    for c in ("odds_ratio", "or_ci_low", "or_ci_high"):
        show[c] = show[c].round(2)
    show["pvalue"] = show["pvalue"].apply(lambda p: f"{p:.2e}" if pd.notna(p) else "NA")
    table = show[["xqtl_kind", "module_set", "n_foreground", "rate_fg", "rate_bg",
                  "odds_ratio", "or_ci_low", "or_ci_high", "pvalue", "method"]]
    lines = [
        f"# Genetic anchoring — IsoGraph co-switch modules vs GTEx {tissue} xQTL "
        f"({analysis}" + (f"/{region}" if region else "") + ")",
        "",
        "Power-matched enrichment (logistic: qtl status ~ module membership + "
        "log cis-variant count + log gene length + log isoform count [+ log intron "
        "group size for sQTL]) within each xQTL's tested-gene universe intersected "
        "with IsoGraph's tested genes.",
        "",
        "Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring "
        f"--analysis {analysis}" + (f" --region {region}" if region else "") + "`",
        "",
        "## Matched odds ratios",
        "",
        _markdown_table(table),
        "",
        "## Reading",
        "",
        "- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL "
        "OR is near 1** for the same module set => the genetic signal on co-switch "
        "modules is splicing-specific (DTU-without-DGE) rather than expression-level.",
        "- Matching on cis-variant count / gene length / isoform multiplicity controls "
        "the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw "
        "(unmatched) contrast for reference.",
        "- Scope: cis-sQTL enrichment shows module *members* undergo genetically "
        "regulated splicing; it does not by itself prove the *co-switching* is genetic "
        "(a shared trans regulator / cell composition could coordinate it).",
    ]
    (out_dir / "QTL_ANCHORING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="GTEx sQTL/eQTL anchoring of co-switch modules.")
    p.add_argument("--analysis", required=True, help="brainseq-sczd, brainseq-aging, gtex-aging")
    p.add_argument("--region", default=None)
    p.add_argument("--variant", default="standard")
    p.add_argument("--fdr", type=float, default=_QVAL)
    p.add_argument("--xqtl-dir", default=str(DEFAULT_XQTL_DIR))
    p.add_argument("--tissue", default=None, help="override GTEx tissue prefix")
    p.add_argument("--gtf-cache", default=str(DEFAULT_GTF_CACHE))
    args = p.parse_args()
    gate = run_anchoring(args.analysis, args.region, args.variant, args.fdr,
                         Path(args.xqtl_dir), args.tissue, Path(args.gtf_cache))
    s = gate[gate.xqtl_kind == "sQTL"].set_index("module_set")["odds_ratio"]
    e = gate[gate.xqtl_kind == "eQTL"].set_index("module_set")["odds_ratio"]
    print("module_set            sQTL_OR  eQTL_OR")
    for k in s.index:
        print(f"{k:22}{s[k]:7.2f}  {e.get(k, float('nan')):7.2f}")


if __name__ == "__main__":
    main()
