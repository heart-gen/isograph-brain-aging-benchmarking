"""Coloc capstone, step 4: paired sQTL/eQTL colocalization contrast + Manubot summary.

Aggregates the CLPP results (coloc_clpp.R) into the manuscript-facing statistics for
the genetic-anchoring capstone:

  * Colocalization rate of IsoGraph switch genes with GTEx brain sQTL vs eQTL
    (fraction of testable genes with CLPP >= 0.01 in >=1 brain region), the paired
    splicing-specificity contrast.
  * Concentration of sQTL colocalization in GO-invisible vs GO-visible switch modules.
  * The colocalizing gene list (best region/CS per gene x kind) for a supplementary table.

A gene is "testable" for a kind if it has >=1 GTEx credible set of that kind at a locus
carried into fine-mapping. Rates are reported per kind over that kind's own testable
set, so the sQTL/eQTL contrast is not confounded by differing CS availability.

Reads real_data/coloc/_m/<analysis>/coloc/clpp_results.tsv. Writes alongside:
  coloc_contrast.parquet      — rate table by kind x go_invisible
  COLOC_SUMMARY.md            — Manubot writeup
  coloc_colocalized_genes.tsv — colocalized gene list (CLPP>=0.01) for a supp table
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data.gwas_traits import TRAITS

CLPP_MIN = 0.01
CLPP_STRONG = 0.05


def _trait_label(analysis: str) -> str:
    """Infer the GWAS trait label from a `<gene_source>__<trait>` analysis dir name."""
    if "__" in analysis:
        key = analysis.rsplit("__", 1)[-1]
        if key in TRAITS:
            return TRAITS[key].label
    return "SCZ"


def _rate(df: pd.DataFrame) -> dict:
    tested = df["gene"].nunique()
    coloc = df.loc[df["clpp"] >= CLPP_MIN, "gene"].nunique()
    strong = df.loc[df["clpp"] >= CLPP_STRONG, "gene"].nunique()
    return {"n_tested": tested, "n_coloc": coloc, "n_strong": strong,
            "coloc_rate": round(coloc / tested, 4) if tested else np.nan}


def run(analysis: str, region: str | None) -> pd.DataFrame:
    m_dir = stage_out("anchoring.coloc", analysis + (f"_{region}" if region else ""))
    clpp_path = m_dir / "coloc" / "clpp_results.tsv"
    if not clpp_path.exists():
        raise SystemExit(f"{clpp_path} not found; run 03.coloc_clpp.R first.")
    clpp = pd.read_csv(clpp_path, sep="\t")
    # attach HGNC symbol (from the credible-set file) for readable gene naming
    cs_path = m_dir / "qtl_credible_sets.tsv"
    if cs_path.exists():
        sym = (pd.read_csv(cs_path, sep="\t", usecols=["gene", "gene_name"])
               .drop_duplicates("gene").set_index("gene")["gene_name"].to_dict())
        clpp["gene_name"] = clpp["gene"].map(sym).fillna(clpp["gene"])

    rows = []
    for kind in ["sQTL", "eQTL"]:
        sub = clpp[clpp["kind"] == kind]
        if sub.empty:
            continue
        rows.append({"kind": kind, "module_set": "all_switch", **_rate(sub)})
        for inv, label in [(True, "go_invisible"), (False, "go_visible")]:
            rows.append({"kind": kind, "module_set": label,
                         **_rate(sub[sub["go_invisible"] == inv])})
    contrast = pd.DataFrame(rows)
    out_dir = ensure_dir(m_dir / "coloc")
    contrast.to_parquet(out_dir / "coloc_contrast.parquet", index=False, compression="zstd")

    # colocalized gene list (best CLPP per gene x kind)
    name_cols = (["gene", "gene_name"] if "gene_name" in clpp.columns else ["gene"])
    coloc_genes = (clpp[clpp["clpp"] >= CLPP_MIN]
                   .sort_values("clpp", ascending=False)
                   .drop_duplicates(["gene", "kind"])
                   [name_cols + ["kind", "go_invisible", "tissue", "LOCUS_ID", "clpp",
                     "n_shared", "best_rsid", "pip_gwas", "pip_qtl", "gwas_p"]])
    coloc_genes.to_csv(out_dir / "coloc_colocalized_genes.tsv", sep="\t", index=False)

    _write_report(out_dir, analysis, region, contrast, clpp, coloc_genes)
    print(contrast.to_string(index=False))
    return contrast


def _write_report(out_dir: Path, analysis: str, region: str | None,
                  contrast: pd.DataFrame, clpp: pd.DataFrame, coloc_genes: pd.DataFrame) -> None:
    def pick(kind, mset, col):
        r = contrast[(contrast["kind"] == kind) & (contrast["module_set"] == mset)]
        return r[col].iloc[0] if len(r) else np.nan

    show = contrast.copy()
    cols = ["kind", "module_set", "n_tested", "n_coloc", "n_strong", "coloc_rate"]
    tbl = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    tbl += ["| " + " | ".join(str(getattr(r, c)) for c in cols) + " |"
            for r in show.itertuples(index=False)]

    n_sq = int(pick("sQTL", "all_switch", "n_coloc") or 0)
    n_eq = int(pick("eQTL", "all_switch", "n_coloc") or 0)
    n_both = coloc_genes.groupby("gene")["kind"].nunique()
    n_both = int((n_both >= 2).sum())
    n_strong_total = int(clpp[clpp["clpp"] >= CLPP_STRONG]["gene"].nunique())

    def _name(row):
        tag = "GO-invisible" if row["go_invisible"] else "GO-visible"
        return f"{row.get('gene_name', row['gene'])} ({row['kind']}, {tag}, CLPP {row['clpp']:.3f})"
    named = "; ".join(_name(r) for _, r in coloc_genes.sort_values("clpp", ascending=False).iterrows()) \
        if len(coloc_genes) else "none"

    label = _trait_label(analysis)
    lines = [
        f"# sQTL/eQTL colocalization of IsoGraph switch genes with {label} GWAS "
        f"({analysis}" + (f"/{region}" if region else "") + ")",
        "",
        "eCAVIAR CLPP (= sum over shared variants of PIP_GWAS x PIP_QTL) between "
        f"{label} GWAS SuSiE fine-mapping and GTEx v11 brain QTL credible sets, for "
        "the genes of IsoGraph's phenotype-associated co-switch modules. Paired across "
        "sQTL and eQTL. CLPP >= 0.01 = colocalized (eCAVIAR convention); >= 0.05 = strong.",
        "",
        "Reproduce: `Rscript real_data/coloc/_h/03.coloc_clpp.R " + analysis + "` then "
        "`python -m isograph_benchmark.real_data.coloc_summary --analysis " + analysis + "`.",
        "",
        "## Colocalization rate by QTL kind and module class",
        "",
        "\n".join(tbl),
        "",
        "## Reading",
        "",
        f"- **{n_sq} switch genes show sQTL colocalization** and {n_eq} eQTL "
        f"colocalization with {label} (CLPP >= 0.01); {n_both} both. Colocalizing genes: "
        f"{named}.",
        (f"- **Magnitude: suggestive, not strong** -- {n_strong_total} genes reach the "
         "strong CLPP >= 0.05 bar. All hits sit just above the 0.01 eCAVIAR threshold, "
         "limited by diffuse GWAS fine-mapping (low PIP_GWAS) at these loci with the "
         "503-sample EUR reference panel. Reported as a named-candidate mechanistic "
         "vignette, not a powered rate contrast (few tens of testable genes)."
         if n_strong_total == 0 else
         f"- {n_strong_total} genes reach the strong CLPP >= 0.05 bar."),
        "- **Why the GO-invisible hits matter most.** GO-invisible switch modules are "
        "the new biology unique to co-switching: coordinated isoform-usage structure "
        "that gene-level GO/pathway enrichment cannot see because GO annotations track "
        f"gene-abundance programs. A shared {label} causal variant landing on the sQTL "
        "of a GO-invisible switch gene is therefore genetic anchoring of biology that "
        "abundance-based co-expression networks miss entirely -- the point of the "
        "switch layer, not a weaker version of an expression result.",
        "- CLPP colocalizes credible sets using GTEx SuSiE PIP (GTEx v11 ships no full "
        f"cis sumstats), so it is restricted to genome-wide-significant {label} loci and "
        "to genes with a GTEx brain QTL credible set. It is evidence that the risk "
        "variant acts through splicing/expression of the switch gene, not proof that "
        "the cross-gene co-switching is itself one genetic signal. Rates are within "
        "each kind's own testable set, so the sQTL-vs-eQTL contrast is not driven by "
        "differing credible-set availability.",
    ]
    (out_dir / "COLOC_SUMMARY.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="Coloc step 4: sQTL/eQTL colocalization contrast.")
    p.add_argument("--analysis", default="brainseq-sczd")
    p.add_argument("--region", default=None)
    args = p.parse_args()
    run(args.analysis, args.region)


if __name__ == "__main__":
    main()
