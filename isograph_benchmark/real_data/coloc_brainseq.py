"""Signal-level colocalization of the BrainSEQ switch and abundance QTLs (EA-only arm).

WHAT THIS ADDS
--------------
The GTEx layer (`coloc_signal_susie`) colocalizes the disease GWAS with GTEx sQTL and eQTL.
BrainSEQ is the discovery cohort, so the same question asked of BrainSEQ is same-tissue genetic
anchoring, never replication: does the disease signal share a causal variant with the QTL of
the switch coordinate S_g itself -- the phenotype IsoGraph modelled -- and with the gene's
abundance channel A_g, in the same donors, variants, cis window and covariates? Each gene has
one S_g and one A_g, so the axis contrast is paired, with no maximum over introns on either
side.

WHY ONLY THE EA-ONLY ARM
------------------------
The GWAS are European, and the GWAS side of every fit is the stage-A SuSiE cache fine-mapped on
the 1000G EUR panel. The all-samples arm is roughly half African-ancestry, so its QTL signals
sit on a different LD structure from the GWAS they would be colocalized with. That arm is
refused outright. `ea_only` is accepted only once `brainseq_qtl_checks` passed in every region:
the recomputed S_g reproduces the discovery axis, and the A_g arm recovers GTEx eGenes.

THE LD, WHICH IS WHAT BrainSEQ BUYS OVER GTEx
---------------------------------------------
GTEx ships no genotypes, so the GTEx QTL side is fine-mapped on reference LD and needs the
`cs_matches_gtex` agreement filter to keep reference-LD artefacts out. BrainSEQ's genotypes are
here, so the QTL side is fine-mapped on LD computed in-sample, from exactly the donors the region
was mapped in, as the correlation of ALT dosages residualized on that axis's covariate matrix --
the genotypes tensorQTL's z-scores actually came from. Raw-genotype LD is not a match at these
sizes (41-42 covariates on 169-229 donors): on DLPFC it doubled estimate_s_rss and made
susie_rss abort in 54 cells. There is no external credible set to agree with, and the reason the
filter existed does not apply, so every QTL signal is primary. The GWAS side is unchanged -- the
same cached fit the GTEx layer colocalizes against -- so the two layers differ in their QTL alone.

WHAT AN S_g COLOCALIZATION DOES NOT SAY
---------------------------------------
S_g is PC1 of a gene's within-gene composition. A colocalization on it says the disease signal
shares a variant with the gene's switch coordinate. It does not name an intron or an event --
that is the GTEx all-introns arm and the event audit -- and the sign of S_g plays no part,
because coloc is invariant to phenotype orientation.

The hierarchy is the GTEx layer's: coloc.susie where both traits fine-map, coloc.abf elsewhere,
flagged. The p12 sweep, prior-robustness and fallback-reason descriptors, the gene collapse and
the paired tests are that layer's functions, called unchanged, with the BrainSEQ region in the
`tissue` column and the axis in `modality`.

STAGES
------
  --stage prep   targets (the GTEx switch-gene grid, restricted to genes BrainSEQ mapped), the
                 per-region phenotype ids, donor lists and n, and the (analysis, region) work
                 list. Login-node safe.
  [R]            05_genetic_anchoring/_h/29.coloc_brainseq_susie.R <analysis> <region>
  --stage meta   hierarchy, per-region cells, gene collapse over regions, the paired S_g vs
                 A_g test, the GTEx nominations read in BrainSEQ, and BrainSEQ's own
                 nominations. Login-node safe.

Outputs under 05_genetic_anchoring/_m/coloc_brainseq/<arm>/.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data import coloc_modality_contrast as cmc
from isograph_benchmark.real_data import coloc_signal_susie as css
from isograph_benchmark.real_data.brainseq_qtl_checks import GTEX_TISSUE, arm_check_failures
from isograph_benchmark.real_data.brainseq_switch_qtl import REGIONS
from isograph_benchmark.real_data.brainseq_switch_qtl import out_dir as qtl_dir

ARMS: tuple[str, ...] = ("ea_only",)
AXES: tuple[str, ...] = ("S_g", "A_g")
FTYPE: dict[str, str] = {"S_g": "switch", "A_g": "abundance"}

# Identical to the GTEx layer, so the two are read on one scale.
PP4_CALL = css.PP4_CALL
P12_PRIMARY = css.P12_PRIMARY
MIN_SHARED_SNPS = css.MIN_SHARED_SNPS
# `tissue` carries the BrainSEQ region throughout, which is what lets the GTEx layer's cell
# functions apply without a parallel implementation.
CELL_KEYS = ["analysis", "trait", "LOCUS_ID", "gene", "tissue", "modality"]
_NOM_KEYS = ["analysis", "trait", "LOCUS_ID", "gene"]


def coloc_root(arm: str = "ea_only") -> Path:
    return stage_out("anchoring.coloc_brainseq") / arm


def check_arm(arm: str, failures: list[str] | None = None) -> None:
    """Refuse any arm colocalization cannot use defensibly. `failures` defaults to disk."""
    if arm not in ARMS:
        raise SystemExit(
            f"BrainSEQ colocalization runs on the `ea_only` arm only, not {arm!r}: the GWAS and "
            "the stage-A GWAS SuSiE LD are European, and the all-samples arm is roughly half "
            "African-ancestry.")
    failures = arm_check_failures(arm) if failures is None else failures
    if failures:
        raise SystemExit(f"BrainSEQ `{arm}` has not passed its QTL checks: "
                         + "; ".join(failures) + " (run 28.brainseq_qtl_checks.sh)")


# --------------------------------------------------------------------------- #
# Stage: prep
# --------------------------------------------------------------------------- #
def phenotype_table(arm: str, regions=REGIONS) -> pd.DataFrame:
    """(region, axis, gene, versioned phenotype id) for every phenotype the arm mapped."""
    rows = []
    for region in regions:
        for axis, ftype in FTYPE.items():
            f = qtl_dir(arm, region) / "qtl" / f"cis_qtl_{ftype}.parquet"
            if not f.exists():
                raise SystemExit(f"missing {f}; run 25.brainseq_switch_qtl.sh for {arm}")
            ids = pd.read_parquet(f, columns=["phenotype_id"])["phenotype_id"].astype(str)
            rows.append(pd.DataFrame({"region": region, "modality": axis,
                                      "gene": ids.str.split(".", n=1).str[0].to_numpy(),
                                      "phenotype_id": ids.to_numpy()}))
    return pd.concat(rows, ignore_index=True)


def region_table(arm: str, regions=REGIONS, dest: Path | None = None) -> pd.DataFrame:
    """Donors and n per region: the donors each region was actually mapped in.

    The in-sample LD must come from exactly these donors, and `n` for susie_rss and coloc.abf
    is their count. The covariate matrix tensorQTL was given is the record of them, so it is
    read rather than re-derived, and both axes must agree on it.
    """
    rows = []
    for region in regions:
        q = qtl_dir(arm, region) / "qtl"
        donors = {}
        for axis, ftype in FTYPE.items():
            f = q / f"covariates_used_{ftype}.txt"
            if not f.exists():
                raise SystemExit(f"missing {f}; run 25.brainseq_switch_qtl.sh for {arm}")
            donors[axis] = [str(x) for x in pd.read_csv(f, sep="\t", index_col=0).index]
        if set(donors["S_g"]) != set(donors["A_g"]):
            raise SystemExit(f"{region}: the S_g and A_g arms were mapped in different donors")
        n_map = set(pd.read_csv(q / "map_summary.tsv", sep="\t")["n_donors"].astype(int))
        if n_map != {len(donors["S_g"])}:
            raise SystemExit(f"{region}: map_summary n_donors {sorted(n_map)} does not match "
                             f"the {len(donors['S_g'])} donors in the covariates")
        name = f"donors_{region}.txt"
        if dest is not None:
            (dest / name).write_text("\n".join(sorted(donors["S_g"])) + "\n")
        rows.append({"region": region, "n_donors": len(donors["S_g"]), "donors_file": name,
                     "gtex_tissue": GTEX_TISSUE[region]})
    return pd.DataFrame(rows)


def work_list(targets: pd.DataFrame, regions=REGIONS) -> pd.DataFrame:
    """One array task per (analysis, region), in a fixed order so a task id names one cell."""
    return pd.DataFrame([(a, r) for a in sorted(targets["analysis"].unique()) for r in regions],
                        columns=["analysis", "region"])


def run_prep(arm: str = "ea_only") -> Path:
    check_arm(arm)
    src = css.signal_root("switch")
    need = [src / "targets.parquet", src / "gwas_meta.tsv"]
    missing = [str(p) for p in need if not p.exists()]
    if missing:
        raise SystemExit(f"missing {missing}; run `coloc_signal_susie --stage prep` first")
    dest = ensure_dir(coloc_root(arm))

    ph = phenotype_table(arm)
    targets = pd.read_parquet(need[0])
    t = targets[targets["gene"].isin(set(ph["gene"]))].reset_index(drop=True)
    gw = src / "gwas_susie"
    cached = {p.name for p in gw.iterdir() if p.is_dir()} if gw.exists() else set()
    uncached = sorted(set(t["analysis"]) - cached)
    if uncached:
        raise SystemExit(f"no stage-A GWAS SuSiE cache for {uncached}; run 22.coloc_gwas_susie.sh")

    t.to_parquet(dest / "targets.parquet", index=False)
    ph[ph["gene"].isin(set(t["gene"]))].to_parquet(dest / "phenotypes.parquet", index=False)
    reg = region_table(arm, dest=dest)
    reg.to_csv(dest / "regions.tsv", sep="\t", index=False)
    pd.read_csv(need[1], sep="\t").to_csv(dest / "gwas_meta.tsv", sep="\t", index=False)
    work = work_list(t)
    work.to_csv(dest / "work_list.tsv", sep="\t", index=False)

    print(f"  {len(t):,} of {len(targets):,} (analysis, locus, gene) targets have a BrainSEQ "
          f"{arm} phenotype ({t['gene'].nunique():,} genes)")
    for r in reg.itertuples(index=False):
        n_s = int(((ph["region"] == r.region) & (ph["modality"] == "S_g")
                   & ph["gene"].isin(set(t["gene"]))).sum())
        n_a = int(((ph["region"] == r.region) & (ph["modality"] == "A_g")
                   & ph["gene"].isin(set(t["gene"]))).sum())
        print(f"  {r.region}: {r.n_donors} donors; target genes with S_g {n_s:,}, A_g {n_a:,}")
    print(f"  work list: {len(work)} (analysis, region) tasks -> {dest / 'work_list.tsv'}")
    return dest


# --------------------------------------------------------------------------- #
# Stage: meta
# --------------------------------------------------------------------------- #
def _load(root: Path, kind: str, required: bool = True) -> pd.DataFrame:
    d = root / kind
    parts = [pd.read_parquet(p) for p in sorted(d.glob("*.parquet"))] if d.exists() else []
    parts = [p for p in parts if len(p)]
    if not parts:
        if required:
            raise SystemExit(f"no {kind}/ shards under {root}; run 29.coloc_brainseq_susie.sh")
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def abf_cells(abf: pd.DataFrame) -> pd.DataFrame:
    """Primary-prior coloc.abf cells at the shared-SNP floor, with their prior-survival range."""
    d = abf[(abf["p12"] == P12_PRIMARY) & (abf["nsnps"] >= MIN_SHARED_SNPS)].copy()
    sw = css.p12_survival(abf[abf["nsnps"] >= MIN_SHARED_SNPS], CELL_KEYS)
    keep = [*CELL_KEYS, "PP3", "PP4", "nsnps",
            *(c for c in ("symbol", "module_id", "go_invisible", "phenotype_id",
                          "n_qtl_donors") if c in d.columns)]
    return d[keep].merge(sw, on=CELL_KEYS, how="left")


def build_hierarchy(pairs: pd.DataFrame, abf: pd.DataFrame,
                    status: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(signal-level cells, full hierarchy). No GTEx agreement filter: the QTL LD is in-sample."""
    susie = (css.best_signal_pair(pairs, require_gtex_match=False) if len(pairs)
             else pd.DataFrame(columns=CELL_KEYS))
    merged = css.apply_hierarchy(susie, abf_cells(abf))
    merged["fallback_reason"] = css.fallback_reasons(merged, status, pairs,
                                                     require_gtex_match=False)
    merged["prior_robustness"] = [css.prior_robustness(x) for x in merged["p12_min_call"]]
    return susie, merged


def pivot_axes(merged: pd.DataFrame, require_both: bool = True) -> pd.DataFrame:
    """Wide (analysis, locus, gene, region) cells carrying both axes.

    `require_both` is what makes the S_g vs A_g contrast paired: a cell enters only when each
    axis cleared the shared-SNP floor. Cross-referencing the GTEx nominations uses the unpaired
    table, because a nomination tested on one axis only is still worth reporting.
    """
    keys = cmc._cell_keys(merged) + ["tissue"]
    d = cmc._fill_key_na(merged, keys)
    vals = [c for c in ("PP3", "PP4", "nsnps", "estimator") if c in d.columns]
    w = d.pivot_table(index=keys, columns="modality", values=vals, aggfunc="first")
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.reset_index()
    for ax in AXES:
        for c in ("PP3", "PP4", "nsnps", "estimator"):
            if f"{c}_{ax}" not in w.columns:
                w[f"{c}_{ax}"] = np.nan
    if require_both:
        w = w.dropna(subset=[f"{c}_{ax}" for c in ("PP3", "PP4") for ax in AXES])
    for ax in AXES:
        den = w[f"PP3_{ax}"].astype(float) + w[f"PP4_{ax}"].astype(float)
        w[f"cond_{ax}"] = np.where(den > 0, w[f"PP4_{ax}"].astype(float) / den, np.nan)
    return w.reset_index(drop=True)


def _to_cmc(d: pd.DataFrame) -> pd.DataFrame:
    return d.rename(columns=lambda c: c.replace("_S_g", "_sQTL").replace("_A_g", "_eQTL"))


def _from_cmc(name: str) -> str:
    """Relabel a GTEx-layer column: sQTL/eQTL become the axes, tissue becomes region."""
    for a, b in (("sQTL", "S_g"), ("eQTL", "A_g"), ("sqtl", "S_g"), ("eqtl", "A_g"),
                 ("splicing", "switch"), ("expression", "abundance"), ("tissue", "region")):
        name = name.replace(a, b)
    return name


def collapse_regions(wide: pd.DataFrame, call: float = PP4_CALL) -> pd.DataFrame:
    """Best region per (analysis, locus, gene), with the per-region consistency beside it."""
    return cmc.collapse_genes(_to_cmc(wide), call=call).rename(columns=_from_cmc)


def contrast_table(genes: pd.DataFrame, call: float = PP4_CALL) -> pd.DataFrame:
    """Exact McNemar and paired Wilcoxon, S_g against A_g, per trait and pooled."""
    g = genes.rename(columns={"PP4_S_g": "PP4_sQTL", "PP4_A_g": "PP4_eQTL",
                              "cond_S_g": "cond_sQTL", "cond_A_g": "cond_eQTL"})
    rows = [{"analysis": a, "trait": t, **cmc.paired_tests(sub, call=call)}
            for (a, t), sub in g.groupby(["analysis", "trait"])]
    rows.append({"analysis": "POOLED", "trait": "ALL", **cmc.paired_tests(g, call=call)})
    return pd.DataFrame(rows).rename(columns=_from_cmc)


def gtex_crossref(gtex_nom: pd.DataFrame, gtex_cells: pd.DataFrame, wide_all: pd.DataFrame,
                  regions=REGIONS, call: float = PP4_CALL) -> pd.DataFrame:
    """Every GTEx signal-level nomination, read in every BrainSEQ region on both axes.

    `gtex_call_in_matched_tissue` says whether the GTEx sQTL call holds in the region's
    tissue-matched GTEx brain tissue. A BrainSEQ region with no matched GTEx call is not a
    failed replication: the GTEx call was made elsewhere, and this layer is anchoring, not
    replication.
    """
    base = gtex_nom[[*_NOM_KEYS, "symbol"]].drop_duplicates(_NOM_KEYS)
    rows = base.merge(pd.DataFrame({"region": list(regions)}), how="cross")
    rows["gtex_tissue"] = rows["region"].map(GTEX_TISSUE)
    gc = (gtex_cells[[*_NOM_KEYS, "tissue", "PP4_sQTL"]]
          .rename(columns={"tissue": "gtex_tissue", "PP4_sQTL": "gtex_PP4_sQTL_matched"})
          .drop_duplicates([*_NOM_KEYS, "gtex_tissue"]))
    rows = rows.merge(gc, on=[*_NOM_KEYS, "gtex_tissue"], how="left")
    rows["gtex_call_in_matched_tissue"] = rows["gtex_PP4_sQTL_matched"].fillna(0) >= call
    cols = [f"{c}_{ax}" for c in ("PP4", "estimator") for ax in AXES]
    w = (wide_all.rename(columns={"tissue": "region"})[[*_NOM_KEYS, "region", *cols]]
         .drop_duplicates([*_NOM_KEYS, "region"]))
    rows = rows.merge(w, on=[*_NOM_KEYS, "region"], how="left")
    rows["brainseq_tested"] = rows["PP4_S_g"].notna() | rows["PP4_A_g"].notna()
    for ax in AXES:
        rows[f"{ax}_coloc"] = rows[f"PP4_{ax}"].astype(float).fillna(0) >= call
    return rows.sort_values([*_NOM_KEYS, "region"]).reset_index(drop=True)


def brainseq_nominations(genes: pd.DataFrame, gtex_nom: pd.DataFrame,
                         call: float = PP4_CALL) -> pd.DataFrame:
    """Loci where the switch coordinate colocalizes in at least one region."""
    nom = genes[genes["PP4_S_g"] >= call].copy()
    g = gtex_nom[_NOM_KEYS].drop_duplicates().assign(in_gtex_nominations=True)
    nom = nom.merge(g, on=_NOM_KEYS, how="left")
    nom["in_gtex_nominations"] = nom["in_gtex_nominations"].notna()
    return nom.sort_values("PP4_S_g", ascending=False).reset_index(drop=True)


def run_meta(arm: str = "ea_only") -> Path:
    root = coloc_root(arm)
    if not (root / "targets.parquet").exists():
        raise SystemExit(f"missing {root / 'targets.parquet'}; run --stage prep first")
    pairs = _load(root, "susie", required=False)
    abf = _load(root, "abf")
    fits = _load(root, "fits", required=False)
    print(f"  loaded {len(pairs):,} coloc.susie signal-pair rows, {len(abf):,} coloc.abf rows")
    pairs.to_parquet(root / "signal_pairs.parquet", index=False)

    status = css.load_gwas_status(css.signal_root("switch"))
    susie, merged = build_hierarchy(pairs, abf, status)
    susie.to_parquet(root / "cells_susie.parquet", index=False)
    merged.to_parquet(root / "cells_hierarchy.parquet", index=False)
    n_s = int((merged["estimator"] == "susie").sum())
    print(f"  hierarchy: {n_s:,} / {len(merged):,} (cell, axis) rows scored by coloc.susie")

    wide = pivot_axes(merged, require_both=True)
    wide.to_parquet(root / "cells.parquet", index=False)
    genes = collapse_regions(wide)
    genes.to_parquet(root / "genes.parquet", index=False)
    contrast = contrast_table(genes)
    contrast.to_parquet(root / "contrast.parquet", index=False)

    from isograph_benchmark.real_data.locus_event_audit import load_nominations

    gnom, gcells = load_nominations("susie", call=PP4_CALL, sqtl_arm="all")
    xref = gtex_crossref(gnom, gcells, pivot_axes(merged, require_both=False))
    xref.to_parquet(root / "gtex_nominations_in_brainseq.parquet", index=False)
    nom = brainseq_nominations(genes, gnom)
    nom.to_parquet(root / "nominations.parquet", index=False)

    _write_report(root, arm, pairs, merged, genes, contrast, xref, nom, fits)
    print(f"  {len(nom)} BrainSEQ S_g nominations "
          f"({int(nom['in_gtex_nominations'].sum())} also GTEx nominations)")
    print(f"  wrote {root}")
    return root


def _fmt(v, nd: int = 3) -> str:
    if v is None or (isinstance(v, (float, np.floating)) and not np.isfinite(v)):
        return "—"
    if isinstance(v, (float, np.floating)):
        return f"{v:.{nd}g}" if (abs(v) < 1e-3 and v != 0) else f"{v:.{nd}f}"
    return str(v)


def _est(v) -> str:
    return "s" if v == "susie" else ("a" if v == "abf" else "")


def _write_report(root: Path, arm: str, pairs: pd.DataFrame, merged: pd.DataFrame,
                  genes: pd.DataFrame, contrast: pd.DataFrame, xref: pd.DataFrame,
                  nom: pd.DataFrame, fits: pd.DataFrame) -> None:
    L: list[str] = []
    A = L.append
    A("# BrainSEQ signal-level colocalization: switch (S_g) and abundance (A_g) QTLs")
    A("")
    A(f"Arm: `{arm}` (European-ancestry donors; both BrainSEQ QTL checks passed in every "
      "region). BrainSEQ is the discovery cohort, so this is **same-tissue genetic "
      "anchoring, not replication**. GTEx is the external cohort.")
    A("")
    A("GWAS side: the stage-A SuSiE cache the GTEx layer uses (1000G EUR LD). QTL side: "
      "`susie_rss` on LD computed **in-sample** from the donors each region was mapped in, as "
      "the correlation of ALT dosages residualized on that axis's covariates (the genotypes the "
      "z-scores came from), so no GTEx-style reference-LD agreement filter is applied. Cells "
      "where either trait does not fine-map fall back to `coloc.abf` and are flagged.")
    A("")
    A("## What was fit")
    A("")
    A(f"- signal-pair posteriors: {len(pairs):,}")
    A(f"- (cell, axis) rows: {len(merged):,}; estimator by region and axis:")
    A("")
    A("| region | axis | coloc.susie | coloc.abf | PP4 >= 0.8 |")
    A("|---|---|---|---|---|")
    for (r, ax), sub in merged.groupby(["tissue", "modality"]):
        A(f"| {r} | {ax} | {int((sub['estimator'] == 'susie').sum()):,} | "
          f"{int((sub['estimator'] == 'abf').sum()):,} | {int((sub['PP4'] >= PP4_CALL).sum()):,} |")
    A("")
    fb = merged.loc[merged["estimator"] == "abf", "fallback_reason"].value_counts()
    if len(fb):
        A("Why coloc.abf, first place the cell left the signal-level pipeline: "
          + ", ".join(f"`{k}` {v:,}" for k, v in fb.items()) + ".")
        A("")
    if len(fits) and "outcome" in fits.columns:
        A("QTL-side fit outcomes: "
          + ", ".join(f"`{k}` {v:,}" for k, v in fits["outcome"].value_counts().items()) + ".")
        A("")
    called = merged[merged["PP4"] >= PP4_CALL]
    if len(called):
        pr = called.groupby(["modality", "estimator", "prior_robustness"]).size()
        A("Prior robustness of the calls at the primary prior: "
          + ", ".join(f"{m} {e} `{lb}` {n:,}" for (m, e, lb), n in pr.items()) + ".")
        A("")
    A("## Paired contrast: S_g against A_g")
    A("")
    A("Cells enter only when both axes cleared the shared-SNP floor; each gene contributes "
      "its best region on both axes symmetrically.")
    A("")
    A("| analysis | trait | genes | S_g coloc | A_g coloc | switch-only | abundance-only | "
      "McNemar P | Wilcoxon (conditional) P |")
    A("|---|---|---|---|---|---|---|---|---|")
    for r in contrast.itertuples(index=False):
        A(f"| {r.analysis} | {r.trait} | {r.n_genes} | {r.n_S_g_coloc} | {r.n_A_g_coloc} | "
          f"{r.switch_only} | {r.abundance_only} | {_fmt(r.mcnemar_p)} | "
          f"{_fmt(r.wilcoxon_cond_p)} |")
    A("")
    A("## The GTEx signal-level nominations, read in BrainSEQ")
    A("")
    A("Each cell is `S_g PP4 / A_g PP4`, suffixed `s` (coloc.susie) or `a` (coloc.abf); "
      "`*` marks a region whose tissue-matched GTEx tissue carries the GTEx sQTL call. "
      "`—` is untested (no BrainSEQ phenotype, or under the shared-SNP floor).")
    A("")
    regions = list(dict.fromkeys(xref["region"]))
    A("| gene | trait | " + " | ".join(regions) + " | S_g call in any region |")
    A("|---|---|" + "---|" * len(regions) + "---|")
    for (sym, trait), sub in xref.groupby(["symbol", "trait"], sort=True):
        cells = []
        for rg in regions:
            r = sub[sub["region"] == rg]
            if r.empty or not bool(r["brainseq_tested"].iloc[0]):
                cells.append("—" + ("*" if len(r) and r["gtex_call_in_matched_tissue"].iloc[0] else ""))
                continue
            r = r.iloc[0]
            cells.append(f"{_fmt(r.PP4_S_g)}{_est(r.estimator_S_g)} / "
                         f"{_fmt(r.PP4_A_g)}{_est(r.estimator_A_g)}"
                         + ("*" if r.gtex_call_in_matched_tissue else ""))
        A(f"| {sym} | {trait} | " + " | ".join(cells) + " | "
          + ("**yes**" if sub["S_g_coloc"].any() else "no") + " |")
    A("")
    A("## BrainSEQ's own switch-axis nominations")
    A("")
    A(f"{len(nom)} (analysis, locus, gene) cells reach PP4_S_g >= {PP4_CALL} in at least one "
      f"region; {int(nom['in_gtex_nominations'].sum())} are also GTEx signal-level nominations. "
      "The count is a set of locus nominations selected on S_g, not a switch-specificity "
      "estimate; the paired test above is the unbiased comparison.")
    A("")
    if len(nom):
        A("| gene | trait | PP4 S_g | PP4 A_g | regions coloc | best region | GTEx nomination |")
        A("|---|---|---|---|---|---|---|")
        for r in nom.head(50).itertuples(index=False):
            A(f"| {r.symbol} | {r.trait} | {_fmt(r.PP4_S_g)} | {_fmt(r.PP4_A_g)} | "
              f"{int(r.n_region_S_g_coloc)}/{int(r.n_region)} | {r.max_region_S_g} | "
              f"{'yes' if r.in_gtex_nominations else 'no'} |")
        A("")
    A("## Reading rules")
    A("")
    A("- An S_g colocalization shares a variant with a gene's switch coordinate (PC1 of its "
      "within-gene composition). It names no intron or event; the GTEx all-introns arm and the "
      "event audit do that.")
    A("- The PP4 per gene is a maximum over three regions; quote it with the region count.")
    A("- PP4_S_g >= 0.8 with low PP4_A_g is switch-preferential colocalization under "
      "prespecified thresholds, not a demonstrated switch-mediated mechanism.")
    A("- Region donors are 169-229, so a non-colocalizing A_g or S_g can be a power result.")
    A("")
    (root / "BRAINSEQ_COLOC.md").write_text("\n".join(L) + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("prep", "meta"), required=True)
    ap.add_argument("--arm", default="ea_only",
                    help="BrainSEQ genotype arm; only ea_only is accepted")
    args = ap.parse_args(argv)
    if args.stage == "prep":
        run_prep(args.arm)
    else:
        check_arm(args.arm)
        run_meta(args.arm)


if __name__ == "__main__":
    main()
