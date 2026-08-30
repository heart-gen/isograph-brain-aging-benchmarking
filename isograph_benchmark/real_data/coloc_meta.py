"""Cross-trait rollup of the coloc capstone: one table over disease + aging traits.

Each `05_genetic_anchoring/_m/coloc/<gene_source>__<trait>/coloc/` holds a per-trait CLPP contrast
(coloc_summary.py). This aggregates them into the manuscript-facing cross-trait view:

  * per-trait colocalization counts (sQTL / eQTL, total / strong, GO-invisible share);
  * the pooled colocalized-gene table (best CLPP per gene x trait x kind), which is the
    named-candidate evidence — canonical neurodegeneration genes (e.g. SNCA for LBD,
    TPP1/SCFD1 for ALS) anchored to switch modules through splicing.

Writes 05_genetic_anchoring/_m/coloc/COLOC_META.md and coloc_meta.parquet (+ colocalized_genes.tsv).
"""
from __future__ import annotations

import argparse

import pandas as pd

from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data.gwas_traits import TRAITS

CLPP_STRONG = 0.05
# gene-source -> case (mirrors switch_bundles / the S-LDSC annotation split)
_CASE = {"brainseq-sczd": "disease", "aging": "aging"}


def _trait_of(dirname: str) -> str:
    return dirname.rsplit("__", 1)[-1]


def _source_of(dirname: str) -> str:
    return dirname.rsplit("__", 1)[0]


def collect() -> tuple[pd.DataFrame, pd.DataFrame]:
    base = stage_out("anchoring.coloc")
    counts, genes = [], []
    for d in sorted(base.glob("*__*")):
        cdir = d / "coloc"
        contrast_f = cdir / "coloc_contrast.parquet"
        genes_f = cdir / "coloc_colocalized_genes.tsv"
        if not contrast_f.exists():
            continue
        trait = _trait_of(d.name)
        case = _CASE.get(_source_of(d.name), "other")
        label = TRAITS[trait].label if trait in TRAITS else trait.upper()
        contrast = pd.read_parquet(contrast_f)
        for kind in ["sQTL", "eQTL"]:
            r = contrast[(contrast.kind == kind) & (contrast.module_set == "all_switch")]
            inv = contrast[(contrast.kind == kind) & (contrast.module_set == "go_invisible")]
            counts.append({
                "trait": label, "case": case, "kind": kind,
                "n_tested": int(r.n_tested.iloc[0]) if len(r) else 0,
                "n_coloc": int(r.n_coloc.iloc[0]) if len(r) else 0,
                "n_strong": int(r.n_strong.iloc[0]) if len(r) else 0,
                "n_coloc_go_invisible": int(inv.n_coloc.iloc[0]) if len(inv) else 0,
            })
        if genes_f.exists():
            g = pd.read_csv(genes_f, sep="\t")
            if len(g):
                g.insert(0, "trait", label)
                g.insert(1, "case", case)
                genes.append(g)
    cnt = pd.DataFrame(counts)
    gdf = (pd.concat(genes, ignore_index=True).sort_values(["case", "trait", "clpp"],
           ascending=[True, True, False]) if genes else pd.DataFrame())
    return cnt, gdf


def _write_report(cnt: pd.DataFrame, gdf: pd.DataFrame, out_dir) -> None:
    lines = [
        "# Cross-trait colocalization of the IsoGraph switch layer",
        "",
        "eCAVIAR CLPP of switch-module genes against GTEx v11 brain sQTL / eQTL "
        "credible sets, for a disease trait (SCZ, on the SCZD switch layer) and the "
        "neurodegenerative aging traits (AD, PD, LBD, ALS, on the pooled aging switch "
        "layer). CLPP >= 0.01 colocalized; >= 0.05 strong.",
        "",
        "## Colocalization counts by trait and QTL kind",
        "",
    ]
    cols = ["case", "trait", "kind", "n_tested", "n_coloc", "n_strong",
            "n_coloc_go_invisible"]
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join("---" for _ in cols) + " |")
    for r in cnt.itertuples(index=False):
        lines.append("| " + " | ".join(str(getattr(r, c)) for c in cols) + " |")
    lines += ["", "## Colocalized genes (best CLPP per gene x trait x kind)", ""]
    if len(gdf):
        gc = ["case", "trait", "gene_name", "kind", "go_invisible", "tissue", "clpp"]
        gc = [c for c in gc if c in gdf.columns]
        lines.append("| " + " | ".join(gc) + " |")
        lines.append("| " + " | ".join("---" for _ in gc) + " |")
        for r in gdf.iterrows():
            row = r[1]
            vals = [f"{row[c]:.3f}" if c == "clpp" else str(row[c]) for c in gc]
            lines.append("| " + " | ".join(vals) + " |")
    n_strong = int((gdf["clpp"] >= CLPP_STRONG).sum()) if len(gdf) else 0
    n_inv = int(gdf["go_invisible"].sum()) if len(gdf) else 0
    lines += [
        "", "## Reading", "",
        f"- Across traits, {len(gdf)} switch genes colocalize (CLPP >= 0.01); "
        f"{n_strong} are strong (>= 0.05); {n_inv}/{len(gdf)} sit in GO-invisible "
        "switch modules.",
        "- The aging case is the stronger one: colocalizing switch genes include "
        "canonical neurodegeneration loci reached through **splicing** (sQTL) of a "
        "GO-invisible switch module -- e.g. SNCA (LBD), TPP1 and SCFD1 (ALS). This is "
        "genetic anchoring of biology that gene-abundance co-expression networks miss, "
        "because GO/pathway enrichment tracks abundance programs, not co-switching.",
        "- Honest scope: CLPP uses GTEx SuSiE credible sets (no full cis sumstats), so "
        "it is restricted to genome-wide-significant trait loci with a brain QTL "
        "credible set, and GTEx brain QTLs are bulk-tissue (cell-type-specific splicing "
        "under-sampled). It shows the risk variant acts through the switch gene's "
        "splicing/expression, not that the co-switching is itself one genetic signal.",
    ]
    (out_dir / "COLOC_META.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    argparse.ArgumentParser(description="Cross-trait coloc rollup.").parse_args()
    out_dir = ensure_dir(stage_out("anchoring.coloc"))
    cnt, gdf = collect()
    if cnt.empty:
        raise SystemExit("No coloc contrasts found.")
    cnt.to_parquet(out_dir / "coloc_meta.parquet", index=False, compression="zstd")
    if len(gdf):
        gdf.to_csv(out_dir / "colocalized_genes.tsv", sep="\t", index=False)
    _write_report(cnt, gdf, out_dir)
    print(cnt.to_string(index=False))
    print(f"\nWrote {out_dir/'COLOC_META.md'}")


if __name__ == "__main__":
    main()
