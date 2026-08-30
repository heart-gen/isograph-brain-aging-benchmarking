"""Coloc capstone, step 1: select IsoGraph switch genes under SCZ peaks + their QTL CS.

Adapts the organoid colocalization machinery (cerebral_organoid_schizophrenia
colocalization/_h/01-03) to the IsoGraph switch layer. Where the organoid pipeline
starts from velocity/WGCNA/rescue candidate genes, this starts from the genes of
IsoGraph's phenotype-associated co-switch modules, and it colocalizes them against
GTEx brain *sQTL and eQTL* credible sets (paired, for the splicing-specificity
contrast) rather than eQTL alone.

The coloc question this sets up: does an SCZ GWAS causal variant share a credible set
with the sQTL of an IsoGraph switch gene (splicing colocalization), and does that
happen more than for the matched eQTL of the same gene? A colocalizing sQTL is the
mechanistic capstone for the genetic-anchoring headline.

Only credible-set-based colocalization is possible: GTEx v11 ships SuSiE fine-mapping
(SuSiE_summary parquet: per phenotype, the 95% credible-set variants with PIP), not
full cis allpairs, so downstream (coloc_clpp.R) uses eCAVIAR CLPP = sum PIP_gwas *
PIP_qtl over shared variants, plus coloc.susie where LD permits. This restricts the
test to genome-wide-significant loci; disclosed as scope.

Matching: IsoGraph modules and GTEx SuSiE both carry Ensembl gene_id, matched bare
(unversioned). hg19 gene coordinates (for the PGC3 window + LD panel, both b37) come
from MAGMA's NCBI37.3.gene.loc keyed by the GTEx symbol. PGC3 EUR is b37, matching.

Writes under real_data/coloc/_m/<analysis>[/<region>]/:
  candidate_genes.tsv   — switch gene, module, go_invisible, symbol, hg19 coords, lead SNP/P, retained
  candidate_loci.tsv    — merged retained ±1Mb windows (LOCUS_ID, chr, start, stop, genes)
  qtl_credible_sets.tsv — gene x tissue x kind(sQTL/eQTL) x CS variant, PIP (b38 ids)
  susie/<locus>.gwas.tsv, susie/<locus>.snps.txt, susie/loci_testable.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel, stage_out
from isograph_benchmark.real_data import gwas_traits as gt
from isograph_benchmark.real_data.qtl_anchoring import (
    _GTEX_TISSUE,
    _bare,
    build_gene_sets,
    resolve_tissue,
)
from isograph_benchmark.real_data.sweep_leiden import _artifact_dir
from isograph_benchmark.real_data.switch_bundles import BUNDLES, get_bundle

GENE_LOC_HG19 = Path("/ocean/projects/bio250020p/shared/opt/magma-v1.10/NCBI37.3.gene.loc")
PANEL_DIR = Path("/ocean/projects/bio250020p/shared/resources/ldsc/1000G_EUR_Phase3_plink")
DEFAULT_XQTL_DIR = rel("inputs", "raw", "gtex_v11", "xqtl")
WINDOW = 1_000_000       # +/- bp around gene body
P_THRESH = 1e-5          # min GWAS signal required in a window
_QVAL = 0.05

# QTL kinds tested, and the SuSiE_summary file suffix for each.
_QTL_KINDS = {"sQTL": "sQTLs", "eQTL": "eQTLs"}
# GTEx brain tissues to pull QTL credible sets from (all 13, for discovery power;
# a switch gene may colocalize in any brain region). Mirrors qtl_anchoring tissue set.
_GTEX_BRAIN = sorted(set(_GTEX_TISSUE.values()))


def _load_gene_loc_hg19() -> pd.DataFrame:
    d = pd.read_csv(GENE_LOC_HG19, sep="\t", header=None,
                    names=["entrez", "chr", "start", "stop", "strand", "symbol"],
                    dtype={"chr": str})
    return d.drop_duplicates("symbol")[["symbol", "chr", "start", "stop"]]


def load_switch_genes(iso_dir: Path, fdr: float) -> pd.DataFrame:
    """IsoGraph switch genes with module + GO-invisible tags, one row per (gene, module).

    Uses the same phenotype-associated / GO-invisible module definitions as
    qtl_anchoring.build_gene_sets, but keeps the gene->module mapping so coloc hits
    can be attributed to GO-invisible vs GO-visible modules.
    """
    enrich_path = iso_dir.parent / "module_enrichment" / "isograph_modules.parquet"
    gene_sets = build_gene_sets(iso_dir, enrich_path, fdr)
    modules = pd.read_parquet(iso_dir / "modules.parquet")
    modules["gene"] = _bare(modules["gene_id"])
    modules["module_id"] = modules["module_id"].astype(str)

    sig = gene_sets.get("pheno_sig_modules", set())
    inv = gene_sets.get("go_invisible_modules", set())
    keep = modules[modules["gene"].isin(sig)].copy()
    keep["go_invisible"] = keep["gene"].isin(inv)
    return keep[["gene", "module_id", "go_invisible"]].drop_duplicates()


def load_qtl_credible_sets(xqtl_dir: Path, tissues: list[str], genes: set[str]) -> pd.DataFrame:
    """GTEx brain sQTL + eQTL credible-set variants for the given genes (bare Ensembl).

    Schemas differ: the eQTL SuSiE summary keys the gene in `phenotype_id` (Ensembl,
    versioned) and has no `gene_id`; the sQTL summary keys the intron in
    `phenotype_id` and the gene in `gene_id`. Both carry `gene_name`.
    """
    rows = []
    base = ["phenotype_id", "gene_name", "variant_id", "pip", "cs_id", "cs_size"]
    for kind, suffix in _QTL_KINDS.items():
        gene_col = "gene_id" if kind == "sQTL" else "phenotype_id"
        cols = base + (["gene_id"] if kind == "sQTL" else [])
        for tissue in tissues:
            path = xqtl_dir / f"{tissue}.v11.{suffix}.SuSiE_summary.parquet"
            if not path.exists():
                continue
            d = pd.read_parquet(path, columns=cols)
            d["gene"] = _bare(d[gene_col])
            d = d[d["gene"].isin(genes)]
            if d.empty:
                continue
            d["tissue"] = tissue
            d["kind"] = kind
            rows.append(d[["gene", "gene_name", "phenotype_id", "variant_id",
                           "pip", "cs_id", "cs_size", "tissue", "kind"]])
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    parts = out["variant_id"].str.split("_", expand=True)
    out["chr_hg38"] = parts[0].str.replace("^chr", "", regex=True)
    out["pos_hg38"] = pd.to_numeric(parts[1], errors="coerce").astype("Int64")
    out["ref"] = parts[2]
    out["alt"] = parts[3]
    return out


def _panel_bim(chrom: int) -> pd.DataFrame:
    """1000G EUR hg19 panel bim for one chromosome (rsid, pos, A1, A2)."""
    f = PANEL_DIR / f"1000G.EUR.QC.{chrom}.bim"
    d = pd.read_csv(f, sep="\t", header=None,
                    names=["chr", "rsid", "cm", "pos", "A1", "A2"],
                    dtype={"rsid": str, "A1": str, "A2": str})
    return d[["rsid", "pos", "A1", "A2"]]


def _panel_positions(chroms: list[int]) -> pd.DataFrame:
    """rsid -> hg19 (chr, pos) across the given panel chromosomes."""
    frames = []
    for c in chroms:
        d = _panel_bim(c)[["rsid", "pos"]].copy()
        d["chr"] = c
        frames.append(d)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["rsid", "pos", "chr"])


def _sig_snps(spec: gt.TraitSpec) -> pd.DataFrame:
    """Genome-wide-significant GWAS SNPs (P < P_THRESH) with hg19 chr/pos from the panel.

    Build-agnostic: the trait's rsIDs are joined to the hg19 LD panel bim for
    positions, so hg38 harmonised sumstats need no liftover (and only panel SNPs,
    the ones coloc can actually use, are kept).
    """
    sig = gt.stream_sig_snps(spec, P_THRESH)  # rsid,a1,a2,beta,se,p,n
    pos = _panel_positions(list(range(1, 23)))
    d = sig.merge(pos, on="rsid", how="inner")
    d = d.rename(columns={"rsid": "snp"})
    d["chr"] = d["chr"].astype(str)
    return d[["chr", "pos", "snp", "p"]].rename(columns={"pos": "bp"})


def _merge_loci(retained: pd.DataFrame) -> pd.DataFrame:
    """Merge overlapping ±window gene intervals on the same chromosome into loci."""
    r = retained.sort_values(["chr_int", "win_start"]).reset_index(drop=True)
    out, cur = [], None
    for _, row in r.iterrows():
        if cur is None:
            cur = dict(chr=row["chr"], chr_int=row["chr_int"], start=row["win_start"],
                       stop=row["win_stop"], genes=[row["gene"]],
                       lead_snp=row["lead_snp"], lead_p=row["lead_p"])
            continue
        if row["chr_int"] == cur["chr_int"] and row["win_start"] <= cur["stop"]:
            cur["stop"] = max(cur["stop"], row["win_stop"])
            cur["genes"].append(row["gene"])
            if pd.notna(row["lead_p"]) and (pd.isna(cur["lead_p"]) or row["lead_p"] < cur["lead_p"]):
                cur["lead_snp"], cur["lead_p"] = row["lead_snp"], row["lead_p"]
        else:
            out.append(cur)
            cur = dict(chr=row["chr"], chr_int=row["chr_int"], start=row["win_start"],
                       stop=row["win_stop"], genes=[row["gene"]],
                       lead_snp=row["lead_snp"], lead_p=row["lead_p"])
    if cur is not None:
        out.append(cur)
    loci = pd.DataFrame(out)
    loci["genes"] = loci["genes"].apply(lambda g: ",".join(sorted(set(g))))
    loci = loci.sort_values(["chr_int", "start"]).reset_index(drop=True)
    loci["LOCUS_ID"] = [f"locus{i+1:02d}_chr{c}" for i, c in enumerate(loci["chr"])]
    return loci


def _write_gwas_loci(loci: pd.DataFrame, out_susie: Path, spec: gt.TraitSpec,
                     tmp_dir: Path) -> pd.DataFrame:
    """Per-locus GWAS z + rsID lists, build-agnostic via the hg19 panel.

    For each locus the panel SNPs in its window define the candidate rsIDs; the trait
    sumstats are streamed once for that rsID set (by name, any build) and each SNP's
    hg19 position comes from the panel. z = beta/se; the R step re-signs z to the panel
    REF allele before susie_rss.
    """
    ensure_dir(out_susie)
    locus_pos: dict[str, dict[str, int]] = {}
    union: set[str] = set()
    for _, L in loci.iterrows():
        bim = _panel_bim(int(L["chr_int"]))
        inwin = bim[(bim["pos"] >= L["start"]) & (bim["pos"] <= L["stop"])]
        locus_pos[L["LOCUS_ID"]] = dict(zip(inwin["rsid"], inwin["pos"]))
        union |= set(inwin["rsid"])
    gw = gt.stream_snps_in_set(spec, union, tmp_dir)  # rsid,a1,a2,beta,se,p,n
    gw = gw[gw["rsid"].str.startswith("rs")]

    summ = []
    for _, L in loci.iterrows():
        lid = L["LOCUS_ID"]
        posmap = locus_pos[lid]
        g = gw[gw["rsid"].isin(posmap)].copy()
        g = g[(g["a1"].str.len() == 1) & (g["a2"].str.len() == 1) & (g["se"] > 0)
              & g["beta"].notna()].drop_duplicates("rsid")
        if g.empty:
            continue
        g["pos"] = g["rsid"].map(posmap).astype(int)
        g["z"] = g["beta"] / g["se"]
        g = g.rename(columns={"rsid": "rsid"})
        g[["rsid", "pos", "a1", "a2", "beta", "se", "p", "z"]].to_csv(
            out_susie / f"{lid}.gwas.tsv", sep="\t", index=False)
        (out_susie / f"{lid}.snps.txt").write_text("\n".join(g["rsid"]) + "\n")
        lead = g.loc[g["p"].idxmin()]
        summ.append({"LOCUS_ID": lid, "chr": int(L["chr_int"]), "start": int(L["start"]),
                     "stop": int(L["stop"]), "genes": L["genes"], "n_gwas_snp": len(g),
                     "true_lead_snp": lead["rsid"], "true_lead_bp": int(lead["pos"]),
                     "true_lead_p": float(lead["p"]), "min_neff": int(g["n"].median())})
    S = pd.DataFrame(summ)
    if not S.empty:
        S.to_csv(out_susie / "loci_testable.tsv", sep="\t", index=False)
    return S


def _write_exclude_regions(spec: gt.TraitSpec, out_dir: Path) -> None:
    """hg19 regions the R step drops from fine-mapping (MHC; APOE for AD/LBD)."""
    rows = [{"chr": c, "start": s, "stop": e} for (c, s, e) in spec.exclude_hg19]
    pd.DataFrame(rows).to_csv(out_dir / "exclude_regions.tsv", sep="\t", index=False)


def _resolve_switch_genes(gene_source: str, region: str | None, variant: str,
                          fdr: float, min_recurrence: int) -> pd.DataFrame:
    """Switch genes (gene, module_id, go_invisible) for a single analysis or a bundle.

    For a bundle, keep genes whose switch-module membership recurs in >= min_recurrence
    analyses (the size-controlled core switch layer, same definition as S-LDSC), and
    tag a gene GO-invisible if it is GO-invisible in any contributing analysis.
    """
    if gene_source in BUNDLES:
        pairs = get_bundle(gene_source)
        from collections import Counter
        counts: Counter = Counter()
        frames = []
        for a, r in pairs:
            sw = load_switch_genes(_artifact_dir(a, r, variant), fdr)
            counts.update(set(sw["gene"]))
            frames.append(sw)
        core = {g for g, n in counts.items() if n >= min_recurrence}
        allsw = pd.concat(frames, ignore_index=True)
        allsw = allsw[allsw["gene"].isin(core)]
        inv = allsw.groupby("gene")["go_invisible"].any()
        out = (allsw.drop_duplicates("gene")[["gene", "module_id"]]
               .assign(go_invisible=lambda d: d["gene"].map(inv)))
        print(f"Bundle '{gene_source}' core switch genes (recurrence>={min_recurrence}): "
              f"{len(out)}")
        return out
    iso_dir = _artifact_dir(gene_source, region, variant)
    return load_switch_genes(iso_dir, fdr)


def run(gene_source: str, trait: str, region: str | None, variant: str, fdr: float,
        xqtl_dir: Path, tissues: list[str], min_recurrence: int) -> None:
    spec = gt.get(trait)
    tag = gene_source + (f"_{region}" if region else "") + f"__{trait}"
    out_dir = ensure_dir(stage_out("anchoring.coloc", tag))
    tmp_dir = ensure_dir(stage_out("anchoring.coloc", "_tmp"))

    switch = _resolve_switch_genes(gene_source, region, variant, fdr, min_recurrence)
    genes = set(switch["gene"])
    print(f"IsoGraph switch genes ({gene_source}): {len(genes)}; trait={spec.label}")

    qtl = load_qtl_credible_sets(xqtl_dir, tissues, genes)
    if qtl.empty:
        raise SystemExit("No GTEx brain QTL credible sets for any switch gene.")
    testable = set(qtl["gene"])
    print(f"Switch genes with a GTEx brain QTL credible set: {len(testable)} "
          f"(sQTL {qtl[qtl.kind=='sQTL'].gene.nunique()}, eQTL {qtl[qtl.kind=='eQTL'].gene.nunique()})")
    qtl.to_csv(out_dir / "qtl_credible_sets.tsv", sep="\t", index=False)
    _write_exclude_regions(spec, out_dir)

    # hg19 coords for testable genes via GTEx symbol -> NCBI37.3.gene.loc. Prefer a
    # real HGNC symbol over an Ensembl-id fallback (some eQTL gene_names are Ensembl).
    sym = (qtl.assign(_is_sym=~qtl["gene_name"].str.startswith("ENSG"))
              .sort_values("_is_sym", ascending=False)
              .drop_duplicates("gene")[["gene", "gene_name"]])
    loc = _load_gene_loc_hg19()
    cand = (switch[switch["gene"].isin(testable)]
            .merge(sym, on="gene", how="left")
            .merge(loc, left_on="gene_name", right_on="symbol", how="left"))
    n_unmapped = cand["chr"].isna().sum()
    if n_unmapped:
        print(f"  {n_unmapped} testable gene rows unmapped to hg19 gene.loc (dropped)")
    cand = cand.dropna(subset=["chr", "start", "stop"]).copy()
    cand["chr_int"] = pd.to_numeric(cand["chr"], errors="coerce")
    cand = cand.dropna(subset=["chr_int"])
    cand["chr_int"] = cand["chr_int"].astype(int)
    cand["win_start"] = np.maximum(1, cand["start"] - WINDOW)
    cand["win_stop"] = cand["stop"] + WINDOW

    sig = _sig_snps(spec)
    print(f"{spec.label} SNPs with P < {P_THRESH:g} (panel-mapped): {len(sig):,}")

    def annotate(row):
        hits = sig[(sig["chr"] == str(row["chr_int"]))
                   & (sig["bp"] >= row["win_start"]) & (sig["bp"] <= row["win_stop"])]
        if hits.empty:
            return pd.Series({"n_sig": 0, "lead_snp": None, "lead_p": np.nan})
        lead = hits.loc[hits["p"].idxmin()]
        return pd.Series({"n_sig": len(hits), "lead_snp": lead["snp"], "lead_p": float(lead["p"])})

    cand = pd.concat([cand.reset_index(drop=True),
                      cand.apply(annotate, axis=1).reset_index(drop=True)], axis=1)
    cand["retained"] = cand["n_sig"] > 0
    cand.sort_values(["retained", "lead_p"], ascending=[False, True]).to_csv(
        out_dir / "candidate_genes.tsv", sep="\t", index=False)
    n_ret = int(cand["retained"].sum())
    print(f"Switch genes under a {spec.label} peak (P<{P_THRESH:g} in ±{WINDOW//10**6}Mb): "
          f"{n_ret}/{len(cand)}")

    retained = cand[cand["retained"]].drop_duplicates("gene")
    if retained.empty:
        print(f"No switch gene sits under a {spec.label} peak; no loci to colocalize.")
        (out_dir / "candidate_loci.tsv").write_text("")
        return
    loci = _merge_loci(retained)
    loci.to_csv(out_dir / "candidate_loci.tsv", sep="\t", index=False)
    S = _write_gwas_loci(loci, out_dir / "susie", spec, tmp_dir)
    print(f"Merged loci: {len(loci)}; loci with usable GWAS SNPs: {len(S)}")
    print(f"Wrote prep to {out_dir}")


def main() -> None:
    p = argparse.ArgumentParser(description="Coloc step 1: switch genes under trait peaks + QTL CS.")
    p.add_argument("--gene-source", default="brainseq-sczd",
                   help="an analysis (e.g. brainseq-sczd) or a bundle name (e.g. aging).")
    p.add_argument("--trait", default="scz", help="GWAS trait key (scz/ad/pd/lbd/als/...).")
    p.add_argument("--region", default=None)
    p.add_argument("--variant", default="standard")
    p.add_argument("--fdr", type=float, default=_QVAL)
    p.add_argument("--min-recurrence", type=int, default=3,
                   help="for a bundle gene-source: keep genes recurring in >= this "
                        "many analyses (core switch layer).")
    p.add_argument("--xqtl-dir", default=str(DEFAULT_XQTL_DIR))
    p.add_argument("--tissue", action="append", help="GTEx tissue(s); default all 13 brain.")
    args = p.parse_args()
    tissues = args.tissue or _GTEX_BRAIN
    run(args.gene_source, args.trait, args.region, args.variant, args.fdr,
        Path(args.xqtl_dir), tissues, args.min_recurrence)


if __name__ == "__main__":
    main()
