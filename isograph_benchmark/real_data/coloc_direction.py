"""Signed direction + isoform event for each colocalized switch gene.

A colocalization (coloc_clpp.R) says a GWAS credible set and a GTEx brain QTL credible
set share the same causal variant for a switch gene. That is unsigned: it does not say
which way the risk allele pushes splicing/expression. This module resolves the sign and
names the concrete isoform event, answering the two follow-ups:

  #11 signed direction — for the colocalizing variant, align the GWAS risk-increasing
      allele to the GTEx effect (ALT) allele and report:
        risk allele -> intron usage up/down (sQTL)  /  gene expression up/down (eQTL).
  #12 isoform event    — for sQTL, the GTEx `phenotype_id` IS the LeafCutter intron
      junction (chr:start:end:clu_N_strand:ENSG); we parse it into concrete junction
      coordinates so the "switch" is a named splicing change, not an abstract score.

Inputs (all already on disk, no GWAS/QTL recompute):
  * <analysis>/coloc/coloc_colocalized_genes.tsv  — the colocalized genes + best_rsid.
  * <analysis>/variant_rsid_map.tsv               — hg38 variant_id <-> rsID.
  * GWAS sumstats via gwas_traits.stream_snps_in_set — signed beta + effect allele.
  * GTEx v11 <tissue>.<kind>.signif_pairs.parquet — signed `slope` (per ALT allele).

GTEx slope convention: `slope` is the effect of the ALT allele of `variant_id`
(chrN_pos_REF_ALT_b38) on the phenotype (intron usage for sQTL, expression for eQTL).
The GWAS beta is the effect of the trait registry's A1 (effect) allele. We map both onto
the trait-increasing (risk) allele, so a positive `risk_qtl_effect` means the risk allele
raises intron usage / expression.

Output (<analysis>/coloc/): coloc_direction.parquet + COLOC_DIRECTION.md.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel, stage_out
from isograph_benchmark.real_data import gwas_traits as gt
from isograph_benchmark.real_data.qtl_anchoring import _bare

_COLOC_ROOT = stage_out("anchoring.coloc")
_XQTL_DIR = rel("inputs", "raw", "gtex_v11", "xqtl")
_KIND_SUFFIX = {"sQTL": "sQTLs", "eQTL": "eQTLs"}


def _trait_key(analysis: str) -> str | None:
    """`aging__ad` -> `ad`; `brainseq-sczd__scz` -> `scz`. Returns None if the resolved
    token is not a registered trait (e.g. the legacy `brainseq-sczd` dir with no
    `__<trait>` suffix, superseded by `brainseq-sczd__scz`)."""
    key = analysis.rsplit("__", 1)[-1]
    return key if key in gt.TRAITS else None


def _load_variant_map(analysis_dir: Path) -> dict[str, str]:
    """rsID -> hg38 variant_id (chrN_pos_REF_ALT_b38)."""
    m = pd.read_csv(analysis_dir / "variant_rsid_map.tsv", sep="\t", header=None,
                    names=["variant_id", "rsid"], dtype=str)
    return dict(zip(m["rsid"], m["variant_id"]))


def _gwas_risk(trait: str, rsids: set[str], tmp_dir: Path) -> pd.DataFrame:
    """Per rsID: the trait-increasing (risk) allele + |beta| from the GWAS sumstats.

    stream_snps_in_set returns a1 (effect allele), a2, and signed beta. The risk allele
    is a1 when beta>0 else a2; risk_beta is |beta| (effect of the risk allele).
    """
    spec = gt.get(trait)
    g = gt.stream_snps_in_set(spec, rsids, tmp_dir)
    if g.empty:
        return pd.DataFrame(columns=["rsid", "risk_allele", "gwas_other_allele",
                                     "risk_beta", "gwas_p"])
    g = g.dropna(subset=["beta"]).copy()
    up = g["beta"] > 0
    g["risk_allele"] = g["a1"].where(up, g["a2"]).astype(str).str.upper()
    g["gwas_other_allele"] = g["a2"].where(up, g["a1"]).astype(str).str.upper()
    g["risk_beta"] = g["beta"].abs()
    g["gwas_p"] = g["p"]
    return g[["rsid", "risk_allele", "gwas_other_allele", "risk_beta", "gwas_p"]]


def _load_slopes(kind: str, tissue: str, genes: set[str]) -> pd.DataFrame:
    """GTEx signed slopes for `genes` in one tissue+kind.

    sQTL keys the gene in `group_id` and the intron in `phenotype_id`; eQTL keys the
    gene in `phenotype_id`. slope is the ALT-allele effect on the phenotype.
    """
    path = _XQTL_DIR / f"{tissue}.v11.{_KIND_SUFFIX[kind]}.signif_pairs.parquet"
    if not path.exists():
        return pd.DataFrame()
    gene_col = "group_id" if kind == "sQTL" else "phenotype_id"
    cols = ["phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal"]
    if kind == "sQTL":
        cols = ["phenotype_id", "group_id", "variant_id", "slope", "slope_se",
                "pval_nominal"]
    d = pd.read_parquet(path, columns=cols)
    d["gene"] = _bare(d[gene_col])
    d = d[d["gene"].isin(genes)].copy()
    if d.empty:
        return d
    d["tissue"] = tissue
    d["kind"] = kind
    return d[["gene", "phenotype_id", "variant_id", "slope", "slope_se",
              "pval_nominal", "tissue", "kind"]]


def _parse_intron(phenotype_id: str, kind: str) -> str:
    """sQTL phenotype_id `chr1:999613:999692:clu_47_-:ENSG...` -> `chr1:999613-999692(-)`.

    The intron junction is the concrete splicing event; for eQTL there is no intron
    (the phenotype is the whole gene) so we return an empty string.
    """
    if kind != "sQTL" or not isinstance(phenotype_id, str):
        return ""
    parts = phenotype_id.split(":")
    if len(parts) < 4:
        return phenotype_id
    chrom, start, end, clu = parts[0], parts[1], parts[2], parts[3]
    strand = clu.split("_")[-1] if "_" in clu else ""
    return f"{chrom}:{start}-{end}({strand})"


def run(analysis: str) -> pd.DataFrame:
    analysis_dir = _COLOC_ROOT / analysis
    coloc_dir = analysis_dir / "coloc"
    cg_path = coloc_dir / "coloc_colocalized_genes.tsv"
    if not cg_path.exists():
        raise SystemExit(f"no colocalized genes for {analysis}: {cg_path} missing")
    cg = pd.read_csv(cg_path, sep="\t")
    if cg.empty:
        print(f"{analysis}: no colocalized genes; nothing to sign.")
        return cg

    trait = _trait_key(analysis)
    if trait is None:
        print(f"{analysis}: unrecognized trait suffix; skipping.")
        return pd.DataFrame()
    rsid2vid = _load_variant_map(analysis_dir)

    tmp_dir = ensure_dir(analysis_dir / "tmp")
    gwas = _gwas_risk(trait, set(cg["best_rsid"].dropna()), tmp_dir)

    # GTEx slopes for the exact (tissue, kind) pairs that appear among coloc hits.
    slopes = []
    for (tissue, kind), sub in cg.groupby(["tissue", "kind"]):
        s = _load_slopes(kind, tissue, set(sub["gene"]))
        if not s.empty:
            slopes.append(s)
    slopes = (pd.concat(slopes, ignore_index=True) if slopes
              else pd.DataFrame(columns=["gene", "phenotype_id", "variant_id", "slope",
                                         "slope_se", "pval_nominal", "tissue", "kind"]))

    cg = cg.copy()
    cg["variant_id"] = cg["best_rsid"].map(rsid2vid)
    cg = cg.merge(gwas, left_on="best_rsid", right_on="rsid", how="left")

    # Join slopes on (gene, tissue, kind, variant_id): the colocalizing variant's effect
    # on each phenotype (intron for sQTL -> one row per resolved isoform event).
    out = cg.merge(slopes, on=["gene", "tissue", "kind", "variant_id"], how="left")

    # Allele alignment: variant_id = chrN_pos_REF_ALT_b38; slope is per ALT.
    vparts = out["variant_id"].str.split("_", expand=True)
    out["ref"] = vparts[2].str.upper() if vparts.shape[1] > 2 else pd.NA
    out["alt"] = vparts[3].str.upper() if vparts.shape[1] > 3 else pd.NA

    risk = out["risk_allele"].astype("string")
    on_alt = risk == out["alt"]
    on_ref = risk == out["ref"]
    aligned = np.where(on_alt, 1.0, np.where(on_ref, -1.0, np.nan))
    out["risk_qtl_effect"] = out["slope"] * aligned
    out["allele_match"] = np.where(on_alt | on_ref, "matched", "ambiguous")
    out["intron_event"] = [
        _parse_intron(p, k) for p, k in zip(out["phenotype_id"], out["kind"])
    ]

    def _direction(row) -> str:
        e = row["risk_qtl_effect"]
        if pd.isna(e):
            return "unresolved (variant not in GTEx signif_pairs or allele mismatch)"
        phen = "intron usage" if row["kind"] == "sQTL" else "expression"
        arrow = "increases" if e > 0 else "decreases"
        tgt = row["intron_event"] or row["gene_name"]
        return f"risk allele {row['risk_allele']} {arrow} {phen} of {tgt}"

    out["direction"] = out.apply(_direction, axis=1)
    out["analysis"] = analysis
    out["trait"] = trait

    keep = ["analysis", "trait", "gene", "gene_name", "kind", "go_invisible", "tissue",
            "best_rsid", "variant_id", "ref", "alt", "risk_allele", "risk_beta",
            "gwas_p", "phenotype_id", "intron_event", "slope", "slope_se",
            "risk_qtl_effect", "allele_match", "clpp", "direction"]
    out = out[[c for c in keep if c in out.columns]].sort_values(
        ["kind", "clpp"], ascending=[True, False]).reset_index(drop=True)

    out.to_parquet(coloc_dir / "coloc_direction.parquet", index=False)
    _write_report(coloc_dir, analysis, out)
    print(f"{analysis}: signed {out['risk_qtl_effect'].notna().sum()}/{len(out)} "
          f"colocalized (gene,phenotype) rows -> {coloc_dir/'coloc_direction.parquet'}")
    return out


def _write_report(coloc_dir: Path, analysis: str, out: pd.DataFrame) -> None:
    resolved = out[out["risk_qtl_effect"].notna()]
    lines = [
        f"# Signed colocalization direction — {analysis}", "",
        "Each colocalized switch gene, with the GWAS risk (trait-increasing) allele "
        "aligned to the GTEx effect allele. `risk_qtl_effect` > 0 means the risk allele "
        "raises intron usage (sQTL) or expression (eQTL); for sQTL, `intron_event` is the "
        "concrete LeafCutter junction (the resolved isoform event).", "",
        f"- colocalized (gene, phenotype) rows: **{len(out)}**",
        f"- signed (allele-matched + slope found): **{len(resolved)}**",
        f"- GO-invisible among signed: **{int(resolved['go_invisible'].sum())}**", "",
    ]
    if not resolved.empty:
        show = resolved.head(20)
        lines += ["| gene | kind | tissue | risk | effect | isoform event | CLPP |",
                  "|------|------|--------|------|--------|---------------|------|"]
        for _, r in show.iterrows():
            sign = "↑" if r["risk_qtl_effect"] > 0 else "↓"
            ev = r["intron_event"] or "(gene-level)"
            lines.append(
                f"| {r['gene_name']} | {r['kind']} | {r['tissue']} | {r['risk_allele']} "
                f"| {sign} | {ev} | {r['clpp']:.3f} |")
    (coloc_dir / "COLOC_DIRECTION.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="Signed direction + isoform event per coloc hit.")
    p.add_argument("--analysis", nargs="*", default=None,
                   help="analysis dir(s) under real_data/coloc/_m; default: all with "
                        "a coloc/coloc_colocalized_genes.tsv")
    args = p.parse_args()
    analyses = args.analysis or sorted(
        d.name for d in _COLOC_ROOT.iterdir()
        if (d / "coloc" / "coloc_colocalized_genes.tsv").exists())
    combined = []
    for a in analyses:
        df = run(a)
        if not df.empty:
            combined.append(df)
    if combined:
        allout = pd.concat(combined, ignore_index=True)
        allout.to_parquet(_COLOC_ROOT / "coloc_direction_combined.parquet", index=False)
        print(f"combined -> {_COLOC_ROOT/'coloc_direction_combined.parquet'} "
              f"({len(allout)} rows across {len(combined)} analyses)")


if __name__ == "__main__":
    main()
