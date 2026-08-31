"""Assemble the real-data MAIN-text tables from the analysis parquet ledgers.

Mirror of assemble_supp_tables.py for the main display items. Presentation only:
every number is copied verbatim from the committed source ledgers under
04_module_characterization/_m -- regenerate, do not hand-edit.

Table 2 (Colocalized splicing-led genes) is the promoted main biology table: the
12 genes where a disease-GWAS-colocalizing sQTL resolves onto an IsoGraph switch
pair (the DTU-without-DGE class), assembled from the per-gene deep-dive panel and
its curated literature layer. It backs the biology payoff (Fig 4) and gives the
main text an at-a-glance biology table (all other biology tables are supplementary).

NB the synthetic-benchmark summary is recommended to move to the supplement
(see manuscript/MANUSCRIPT_PLAN.md Sec 15), so Table 2 here is intended as the biology main table.

Run: python manuscript/_h/assemble_main_tables.py   (login node, no SLURM)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from isograph_benchmark.paths import ensure_dir, stage_out  # noqa: E402

OUT = ensure_dir(stage_out("manuscript", "main_tables"))
DEEP = stage_out("anchoring", "deep_dive")
QTL = stage_out("anchoring", "qtl_anchoring_meta")

# eCAVIAR CLPP thresholds used by the coloc pipeline (coloc_summary.py):
# >= 0.01 "colocalized" (standard permissive eCAVIAR bar); >= 0.05 "strong".
# Manuscript confidence tiers (stars): * suggestive >0.01, ** moderate >0.05,
# *** high-confidence >0.10.
CLPP_MIN = 0.01
CLPP_STRONG = 0.05
CLPP_HIGH = 0.10


def _clpp_stars(x: float) -> str:
    """Confidence tier for a CLPP posterior: * >0.01, ** >0.05, *** >0.10."""
    if pd.isna(x):
        return ""
    if x > CLPP_HIGH:
        return "***"
    if x > CLPP_STRONG:
        return "**"
    if x > CLPP_MIN:
        return "*"
    return ""

# Human labels for the module-set contrast rows.
_MODSET = {
    "all_modules": "All modules",
    "pheno_sig_modules": "Phenotype-associated",
    "go_invisible_modules": "GO-invisible (DTU-without-DGE)",
    "go_visible_modules": "GO-visible (immune/abundance)",
}


def _clean_tissue(t: str) -> str:
    """GTEx tissue id -> human label: Brain_Frontal_Cortex_BA9 -> Frontal Cortex BA9."""
    if not isinstance(t, str) or not t:
        return ""
    return t.replace("Brain_", "").replace("_", " ").strip()


def _fmt_clpp(x: float) -> str:
    return "" if pd.isna(x) else f"{x:.3f}"


def _status(curation: str) -> str:
    """Literature curation -> compact status label for the table."""
    return {"documented": "known isoform biology",
            "novel_candidate": "novel candidate"}.get(str(curation), str(curation))


def table2_splicing_led() -> pd.DataFrame:
    """Table 2 -- colocalized splicing-led genes (IsoGraph-resolved switch)."""
    panel = pd.read_parquet(DEEP / "deep_dive_panel.parquet")
    lit = pd.read_parquet(DEEP / "deep_dive_literature.parquet")

    sl = panel[panel["verdict"].str.startswith("splicing-led")].copy()
    sl = sl.merge(
        lit[["gene_name", "curation"]].rename(columns={"gene_name": "gene"}),
        on="gene", how="left")

    tbl = pd.DataFrame({
        "Gene": sl["gene"],
        "Ensembl": sl["ens"],
        "Trait(s)": sl["traits"],
        "Concordant traits": sl["concordant_traits"].fillna(""),
        "Lead variant": sl["top_rsid"],
        "Risk allele": sl["top_risk_allele"],
        "Top tissue": sl["top_tissue"].map(_clean_tissue),
        "Max CLPP": sl["max_clpp"].map(
            lambda x: f"{x:.3f} {_clpp_stars(x)}".strip()),
        "n resolved events": sl["n_resolved_events"],
        "GO-invisible": sl["go_invisible"].map({True: "yes", False: "no"}),
        "LOEUF": sl["loeuf"].round(3),
        "Status": sl["curation"].map(_status),
    })

    # Deterministic order: documented cases first, then by descending CLPP.
    tbl["_doc"] = (sl["curation"].values == "documented").astype(int)
    tbl["_clpp"] = sl["max_clpp"].values
    tbl = (tbl.sort_values(["_doc", "_clpp"], ascending=[False, False])
              .drop(columns=["_doc", "_clpp"]).reset_index(drop=True))
    return tbl


def table_qtl_specificity() -> pd.DataFrame:
    """Statistical anchor: sQTL/eQTL splicing-specificity contrast per module set,
    with the matched-baseline method comparison (the p-value-bearing biology table)."""
    c = pd.read_parquet(QTL / "qtl_anchoring_meta_contrast.parquet")

    def _p(x: float) -> str:
        return "" if pd.isna(x) else f"{x:.2g}"

    rows = []
    # IsoGraph across the four module sets (the primary contrast).
    for _, r in c[c.graph_method == "isograph"].iterrows():
        rows.append({
            "Method": "IsoGraph",
            "Module set": _MODSET.get(r.module_set, r.module_set),
            "Analyses (k)": int(r.k),
            "sQTL/eQTL ratio": round(r.ratio_fe, 3),
            "95% CI": f"{r.ratio_fe_low:.2f}-{r.ratio_fe_high:.2f}",
            "p (FE)": _p(r.p_fe),
            "I2": round(r.I2, 2),
        })
    # Matched-baseline controls on the phenotype-associated + GO-invisible sets.
    for meth in ("wgcna_switch_only", "wgcna_multiplex"):
        for ms in ("pheno_sig_modules", "go_invisible_modules"):
            sub = c[(c.graph_method == meth) & (c.module_set == ms)]
            if sub.empty:
                continue
            r = sub.iloc[0]
            rows.append({
                "Method": meth,
                "Module set": _MODSET.get(ms, ms),
                "Analyses (k)": int(r.k),
                "sQTL/eQTL ratio": round(r.ratio_fe, 3),
                "95% CI": f"{r.ratio_fe_low:.2f}-{r.ratio_fe_high:.2f}",
                "p (FE)": _p(r.p_fe),
                "I2": round(r.I2, 2),
            })
    return pd.DataFrame(rows)


def _to_markdown(df: pd.DataFrame) -> str:
    """Render a GitHub-flavoured pipe table (no tabulate dependency)."""
    cols = [str(c) for c in df.columns]
    head = "| " + " | ".join(cols) + " |"
    rule = "| " + " | ".join("---" for _ in cols) + " |"
    rows = [
        "| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
        for row in df.itertuples(index=False, name=None)
    ]
    return "\n".join([head, rule, *rows])


def write_csv_and_md(df: pd.DataFrame, stem: str, title: str, legend: str) -> None:
    df.to_csv(OUT / f"{stem}.csv", index=False)
    md = [f"# {title}", "", legend, "", _to_markdown(df), ""]
    (OUT / f"{stem}.md").write_text("\n".join(md))
    print(f"  wrote {stem}.csv / {stem}.md   {df.shape[0]} x {df.shape[1]}")


def main() -> None:
    print("Assembling main-text tables ->", OUT)

    # Table 2 -- statistical anchor (p-value-bearing): sQTL/eQTL specificity contrast.
    t2 = table_qtl_specificity()
    # NB the wording below was corrected on 2026-08-29 after the stale downstream
    # outputs were re-run. The previous legend asserted that the effect "concentrates
    # in the GO-invisible modules" at 1.172 / I2 = 0.00; on inputs that agree with
    # their own sources that arm is 1.068, p = 0.077, and GO-invisible no longer
    # separates from GO-visible. Do not restore the old wording.
    legend2 = (
        "**Table 2. Splicing-QTL are spared relative to expression-QTL in IsoGraph's "
        "phenotype-associated co-switch modules -- an IsoGraph-only method effect.** "
        "Paired within-analysis sQTL-odds-ratio / eQTL-odds-ratio contrast (>1 = "
        "splicing genetics spared over expression genetics), inverse-variance "
        "fixed-effect meta-analysis across brain xQTL analyses (k), with 95% CI, p, "
        "and I2 heterogeneity. The contrast removes the shared cis-QTL depletion "
        "baseline of constrained network genes. The effect is carried by the "
        "phenotype-associated set (1.111, p = 3.6e-4); on the 8-tissue set common to "
        "all methods it is 1.108 (p = 0.001) at I2 = 0.00. **It does NOT localise to "
        "the GO-invisible modules:** GO-invisible (1.068, p = 0.077) and GO-visible "
        "(1.084, p = 0.050) are indistinguishable, so this table does not support a "
        "GO-invisible-specific genetic claim -- the DTU-without-DGE content claim "
        "rests on the GO-invisible gate (Fig S-real-3) instead. The primary internal "
        "control is the matched WGCNA baselines, which consume identical switch "
        "features and are null in every module set (p >= 0.41). This is the "
        "statistical anchor of the genetic-anchoring result (cf. per-gene resolution "
        "in Table 3, whose colocalization posteriors are individually modest). "
        "Verbatim from "
        "05_genetic_anchoring/_m/qtl_anchoring_meta/qtl_anchoring_meta_contrast.parquet "
        "(regenerated 2026-08-29)."
    )
    write_csv_and_md(
        t2, "table2_qtl_specificity_contrast",
        "Table 2 -- sQTL/eQTL splicing-specificity contrast", legend2)

    # Table 3 -- per-gene resolution companion (coherence, NOT per-locus significance).
    t3 = table2_splicing_led()
    mc = pd.read_parquet(DEEP / "deep_dive_panel.parquet")
    mc = mc.loc[mc["verdict"].str.startswith("splicing-led"), "max_clpp"]
    n_hi = int((mc > CLPP_HIGH).sum())
    n_mod = int(((mc > CLPP_STRONG) & (mc <= CLPP_HIGH)).sum())
    legend3 = (
        "**Table 3. Genes where a disease-GWAS-colocalizing sQTL resolves onto an "
        "IsoGraph switch pair (splicing-led genes).** The 12 genes at which a brain "
        "sQTL colocalizing with a GWAS credible set (eCAVIAR CLPP) maps onto an "
        "IsoGraph switch pair -- the DTU-without-DGE class; all 12 lie in GO-invisible "
        "modules. Max CLPP confidence tiers: * suggestive (>0.01, the coloc inclusion "
        "bar), ** moderate (>0.05), *** high-confidence (>0.10). **Colocalization "
        "posteriors are individually modest** "
        f"({n_hi} high-confidence, {n_mod} moderate; CTSH is the sole *** case at "
        "0.39): the defensible claim is the *set-level* coherence -- splicing-led "
        "equivalent to GO-invisible, cross-disease concordance (SNCA in LBD+PD) -- not "
        "any single locus. The set-level statistical support is Table 2 (contrast) and "
        "partitioned heritability (S-LDSC), not per-gene CLPP. Risk allele aligned to "
        "the GWAS trait; LOEUF is gnomAD constraint; Status = established disease "
        "isoform biology vs novel candidate. Verbatim from 05_genetic_anchoring/_m/deep_dive/ "
        "(deep_dive_panel.parquet, deep_dive_literature.parquet); see Tables S8-S12 "
        "for per-event/RBP/clinical/literature layers."
    )
    write_csv_and_md(
        t3, "table3_splicing_led_genes",
        "Table 3 -- Splicing-led colocalized genes (per-gene resolution)", legend3)


if __name__ == "__main__":
    main()
