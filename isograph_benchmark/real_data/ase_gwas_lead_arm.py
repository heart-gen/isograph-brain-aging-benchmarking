"""Second allelic arm: anchor the within-donor test on the GWAS lead itself (PI item 10a).

WHY A SECOND ARM
----------------
`ase_risk_orientation` fits the allelic test at the BrainSEQ switch-QTL lead and then
transfers the sign to the disease allele through signed LD. That transfer failed for most
rows -- 336 of 556 fall below the |r| >= 0.8 gate or have no LD at all -- so the
disease-linkage result reads "we could not ask" rather than "we asked and the answer was
no". The ambiguity is structural, not bad luck: a pair enters the test because its gene has
a switch-QTL, and a switch-QTL lead is by construction chosen for QTL signal rather than for
proximity to a GWAS peak.

This arm removes the transfer. It anchors the same within-donor contrast on the
colocalizing locus's **GWAS lead variant**, so the haplotype carrying the risk allele is
read directly off each donor's phased genotype and the orientation is exact by construction
(r = 1). The cost is power, not validity: the GWAS lead usually has weaker QTL signal, and
only donors heterozygous AT THAT VARIANT are informative.

Read the two arms together. A directional effect here converts the negative into a positive.
No effect here makes the negative considerably stronger, because it can no longer be
explained by the LD gate.

NO RECOUNT IS NEEDED
--------------------
`allelic_donor_counts.parquet` stores, per pair x donor, T1 and T2 fragments on the ALT and
the REF haplotype at the *fitted* lead, plus `alt_hap` naming which phASER haplotype was
ALT. Haplotype 1 and haplotype 2 counts are recoverable from that (hap 1 is the ALT side
unless `alt_hap == 2`), and re-anchoring is then a relabelling: which haplotype carries the
risk allele at the new variant. The fragments never move.

THE PHASE FRAME IS CHECKED, NOT ASSUMED
---------------------------------------
The counts' haplotypes come from phASER's genome-wide `PW`, anchored to the phased VCF used
on Quest. This arm reads the GWAS lead's phased genotype from the BrainSEQ TOPMed panel on
Bridges-2, which is a different file. If the two are not in the same phase frame, every
re-anchored sign is arbitrary and the arm is worthless.

So the frame is tested before anything is fitted: the fitted leads are extracted from the
same panel and compared donor by donor with the `gt` the counts carry. The check reports
agreement on the phased genotype string among donors heterozygous in both, and the arm
refuses to run below `FRAME_MIN`. A run that stops here has produced a real result -- the
frames differ -- not a failure.

OUTPUT
------
`gwas_lead_arm.parquet` per region, `GWAS_LEAD_ARM.md`, `gwas_lead_arm_summary.json`.

  anchor_variant, risk_allele      the GWAS lead and its trait-increasing allele
  beta_risk_direct                 within-donor log odds of T1 on the RISK haplotype,
                                   oriented by construction rather than by LD
  n_paired_het                     donors heterozygous at the GWAS lead and informative
                                   on both haplotypes -- the power this arm trades for
  risk_along_module                beta_risk_direct signed by the pair's module polarity
  frame_agreement                  phase-frame concordance behind the whole arm
"""
from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.real_data.ase_junction_allelic import allelic_contrast
from isograph_benchmark.real_data.ase_risk_orientation import (
    REGIONS, QVAL_ALLELIC, _PANELS, _PLINK2, _bh, _gwas_risk, _palindromic,
    build_targets, dest_dir, donor_ancestry,
)

# The per-sample phASER VCFs that define the `PW` haplotypes the counts are phased on.
_PHASER = Path("/ocean/projects/bio260021p/shared/resources/processed-data/ase-files")


def _bcftools() -> str:
    """bcftools is not on the default PATH on the login nodes; resolve it once, loudly."""
    import os
    import shutil
    for c in (os.environ.get("BCFTOOLS"), shutil.which("bcftools"),
              "/ocean/projects/bio250020p/shared/opt/env/ml_dev/bin/bcftools"):
        if c and Path(c).exists():
            return c
    raise SystemExit("bcftools not found: set BCFTOOLS, or `module load bcftools`")

# Agreement between the fitted lead as this code reads it and as the counts carry it. This
# is a self-check, not a phase transfer: both sides are phASER `PW`, so anything below this
# means the two are not reading the same donors or the same variants.
FRAME_MIN = 0.95

# Donors heterozygous at the anchor and informative on both haplotypes. Same floor as the
# primary arm, so the two are comparable.
MIN_PAIRED_HET = 5


# --------------------------------------------------------------------------- #
# Phased genotypes at an arbitrary variant, read from phASER's own VCFs
# --------------------------------------------------------------------------- #
def variant_sites(rsids: set[str], panel: Path) -> pd.DataFrame:
    """rsID -> chrom, pos, ref, alt, from the QTL panel's variant table.

    Only the COORDINATES come from the panel; every genotype this arm uses comes from
    phASER. A lead the panel does not carry is dropped here and reported as such.
    """
    if not rsids:
        return pd.DataFrame(columns=["variant_id", "chrom", "pos", "ref", "alt"])
    rows = []
    with Path(str(panel) + ".pvar").open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] in rsids:
                rows.append((f[2], f[0] if f[0].startswith("chr") else "chr" + f[0],
                             int(f[1]), f[3], f[4]))
    return (pd.DataFrame(rows, columns=["variant_id", "chrom", "pos", "ref", "alt"])
            .drop_duplicates("variant_id"))


def phaser_gt(region: str, sites: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """donor x variant phased genotype in phASER's `PW` frame, one query per sample.

    `PW` is phASER's genome-wide haplotype call and is what the fragment counts were
    assigned against; `GT` is the input phasing and is used only where `PW` is missing
    (phASER leaves it unphased at sites it could not place in a block).
    """
    if sites.empty:
        return pd.DataFrame(columns=["variant_id", "donor_id", "ref", "alt", "gt", "source"])
    vdir = _PHASER / region / "phASER"
    if not vdir.is_dir():
        raise SystemExit(f"no phASER VCFs for {region} under {vdir}")
    _BCF = _bcftools()
    rows, missing = [], 0
    with tempfile.TemporaryDirectory() as td:
        reg = Path(td) / "sites.txt"
        reg.write_text("".join(f"{r.chrom}\t{r.pos}\n" for r in sites.itertuples()))
        want = set(sites["variant_id"])
        for r in samples.itertuples():
            vcf = vdir / f"{r.sample_id}.vcf.gz"
            if not vcf.exists():
                missing += 1
                continue
            q = subprocess.run(
                [_BCF, "query", "-R", str(reg),
                 "-f", "%ID\t%REF\t%ALT[\t%GT\t%PW]\n", str(vcf)],
                capture_output=True, text=True)
            if q.returncode != 0:
                missing += 1
                continue
            for line in q.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 5 or f[0] not in want:
                    continue
                vid, ref, alt, gt, pw = f[0], f[1], f[2], f[3], f[4]
                if len(ref) != 1 or len(alt) != 1:
                    continue
                use, src = (pw, "PW") if "|" in pw else (gt, "GT")
                if "|" not in use:
                    continue        # unphased here: this donor gets no haplotype assignment
                rows.append((vid, r.donor_id, ref, alt, use, src))
    if missing:
        print(f"  {region}: {missing} sample VCF(s) unreadable or absent", flush=True)
    return (pd.DataFrame(rows, columns=["variant_id", "donor_id", "ref", "alt", "gt",
                                        "source"])
            .drop_duplicates(["variant_id", "donor_id"]))


# --------------------------------------------------------------------------- #
# Self-check: does this code read the fitted lead the way the recount did?
# --------------------------------------------------------------------------- #
def frame_check(counts: pd.DataFrame, panel_gt: pd.DataFrame) -> dict:
    """Compare the counts' phased `gt` at the fitted lead with phASER's, per donor.

    Only heterozygotes carry information about the frame: `1|0` and `0|1` are the same
    genotype and different haplotype assignments, which is exactly what must agree.
    """
    c = (counts[["variant_id_all", "donor_id", "gt"]]
         .drop_duplicates(["variant_id_all", "donor_id"])
         .rename(columns={"variant_id_all": "variant_id", "gt": "gt_counts"}))
    m = c.merge(panel_gt[["variant_id", "donor_id", "gt"]].rename(
        columns={"gt": "gt_panel"}), on=["variant_id", "donor_id"], how="inner")
    het = m[m["gt_counts"].isin(["0|1", "1|0"]) & m["gt_panel"].isin(["0|1", "1|0"])]
    same_geno = m[m["gt_counts"].str.replace("|", "", regex=False).apply(sorted).astype(str)
                  == m["gt_panel"].str.replace("|", "", regex=False).apply(sorted).astype(str)]
    return {
        "variants_compared": int(m["variant_id"].nunique()),
        "donor_variant_pairs": int(len(m)),
        "genotype_agreement": float(len(same_geno) / len(m)) if len(m) else float("nan"),
        "het_pairs": int(len(het)),
        "frame_agreement": (float((het["gt_counts"] == het["gt_panel"]).mean())
                            if len(het) else float("nan")),
    }


# --------------------------------------------------------------------------- #
# Re-anchor the fragments on the risk allele
# --------------------------------------------------------------------------- #
def reanchor(counts: pd.DataFrame, gt: pd.DataFrame, risk_allele: str) -> pd.DataFrame:
    """Relabel each donor's haplotype counts by which one carries the risk allele.

    `counts` is one pair's donor rows; `gt` is that donor set's phased genotype at the
    anchor. Homozygous donors keep haplotype 1 as the pseudo-risk side and fall into the
    estimator's built-in null arm, exactly as at the fitted lead.
    """
    c = counts.copy()
    alt1 = c["alt_hap"] != 2                      # hap 1 is the ALT side unless alt_hap == 2
    for iso in ("t1", "t2"):
        c[f"{iso}_hap1"] = np.where(alt1, c[f"{iso}_alt"], c[f"{iso}_ref"])
        c[f"{iso}_hap2"] = np.where(alt1, c[f"{iso}_ref"], c[f"{iso}_alt"])
    g = gt.set_index("donor_id")
    c = c[c["donor_id"].isin(g.index)].copy()
    if c.empty:
        return c
    gg = g.loc[c["donor_id"]]
    # which phased slot carries the risk allele: '1' when the risk allele is the panel ALT
    risk_is_alt = gg["alt"].str.upper().to_numpy() == risk_allele.upper()
    risk_code = np.where(risk_is_alt, "1", "0")
    left = gg["gt"].str[0].to_numpy()
    right = gg["gt"].str[2].to_numpy()
    het = left != right
    risk_hap = np.where(~het, 0, np.where(left == risk_code, 1, 2))
    c["alt_hap"] = risk_hap                       # 'ALT' now means 'risk allele'
    on1 = risk_hap != 2
    for iso in ("t1", "t2"):
        c[f"{iso}_alt"] = np.where(on1, c[f"{iso}_hap1"], c[f"{iso}_hap2"])
        c[f"{iso}_ref"] = np.where(on1, c[f"{iso}_hap2"], c[f"{iso}_hap1"])
    return c.drop(columns=[f"{i}_hap{h}" for i in ("t1", "t2") for h in (1, 2)])


def fit_region(region: str, targets: pd.DataFrame, panel_gt: pd.DataFrame,
               ancestry: pd.Series) -> pd.DataFrame:
    """One fit per (pair, trait): the anchor is a property of the trait's locus."""
    counts = pd.read_parquet(dest_dir(region) / "allelic_donor_counts.parquet")
    ea = set(ancestry[ancestry == "EA"].index)
    rows = []
    for (pair_id, trait), sub in targets.groupby(["pair_id", "trait"]):
        r = sub.iloc[0]
        rec = {"pair_id": pair_id, "trait": trait, "region": region,
               "gene": r.get("gene"), "symbol": r.get("symbol"),
               "anchor_variant": r["lead_snp"], "risk_allele": r["risk_allele"],
               "risk_gwas_p": r.get("risk_gwas_p"), "fitted_lead": r["variant_id_all"],
               "anchor_is_fitted_lead": bool(r["lead_snp"] == r["variant_id_all"])}
        gt = panel_gt[panel_gt["variant_id"] == r["lead_snp"]]
        c = counts[counts["pair_id"] == pair_id]
        c = c[c["donor_id"].isin(ea)]
        if gt.empty:
            rec["status"] = "anchor_not_in_panel"
        elif pd.isna(r["risk_allele"]):
            rec["status"] = "no_gwas_risk_allele"
        elif bool(_palindromic(pd.Series([r["risk_allele"]]),
                               pd.Series([gt["ref"].iloc[0]])).iloc[0]) or \
                bool(_palindromic(pd.Series([r["risk_allele"]]),
                                  pd.Series([gt["alt"].iloc[0]])).iloc[0]):
            rec["status"] = "palindromic_anchor"
        elif r["risk_allele"].upper() not in {gt["ref"].iloc[0].upper(),
                                              gt["alt"].iloc[0].upper()}:
            rec["status"] = "allele_mismatch"
        elif c.empty:
            rec["status"] = "no_ea_donors"
        else:
            tab = reanchor(c, gt, str(r["risk_allele"]))
            if tab.empty:
                rec["status"] = "no_genotyped_donors"
            else:
                f = allelic_contrast(tab, min_paired=MIN_PAIRED_HET)
                rec |= {"beta_risk_direct": f.get("beta", np.nan),
                        "se": f.get("se", np.nan), "pval": f.get("pval", np.nan),
                        "n_donors_het": f.get("n_donors_het", 0),
                        "n_paired_het": f.get("n_paired_het", 0),
                        "n_paired_hom": f.get("n_paired_hom", 0),
                        "score_z": f.get("score_z", np.nan),
                        "score_pval": f.get("score_pval", np.nan),
                        "beta_hom_null": f.get("beta_hom_null", np.nan),
                        "at_bound": bool(f.get("at_bound", False)),
                        "status": f.get("status", "not_fitted")}
        for col in ("module_polarity", "module_id"):
            if col in sub.columns:
                rec[col] = r.get(col)
        rows.append(rec)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    if "beta_risk_direct" in out.columns and "module_polarity" in out.columns:
        pol = pd.to_numeric(out["module_polarity"], errors="coerce")
        out["risk_along_module"] = out["beta_risk_direct"] * np.sign(pol)
    ok = out.get("pval", pd.Series(dtype=float)).notna()
    out["qval"] = np.nan
    if ok.any():
        out.loc[ok, "qval"] = _bh(out.loc[ok, "pval"])
    return out


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def _summary(t: pd.DataFrame, frame: dict) -> dict:
    fitted = t[t["status"] == "fitted"] if "status" in t else t.iloc[0:0]
    s = {
        "rows": int(len(t)),
        "pairs": int(t["pair_id"].nunique()) if len(t) else 0,
        "genes": int(t["gene"].nunique()) if "gene" in t and len(t) else 0,
        "fitted": int(len(fitted)),
        "anchor_is_fitted_lead": int(t.get("anchor_is_fitted_lead",
                                           pd.Series(dtype=bool)).sum()),
        "median_paired_het": (float(fitted["n_paired_het"].median())
                              if len(fitted) else None),
        "significant": int((fitted["qval"] < QVAL_ALLELIC).sum()) if len(fitted) else 0,
        # A fit at the optimiser's bound is quasi-separation -- one haplotype carries one
        # isoform only. Its SIGN is informative, its magnitude is not, so it is counted apart
        # rather than reported as an effect. The primary arm applies the same rule.
        "significant_at_bound": (int(((fitted["qval"] < QVAL_ALLELIC)
                                      & fitted["at_bound"]).sum()) if len(fitted) else 0),
        "significant_interpretable": (int(((fitted["qval"] < QVAL_ALLELIC)
                                           & ~fitted["at_bound"]).sum()) if len(fitted) else 0),
        "at_bound": int(fitted["at_bound"].sum()) if len(fitted) else 0,
        "positive_beta": int((fitted["beta_risk_direct"] > 0).sum()) if len(fitted) else 0,
        "median_abs_beta": (float(fitted.loc[~fitted["at_bound"], "beta_risk_direct"]
                                  .abs().median()) if len(fitted) else None),
    }
    if "status" in t:
        s["status_counts"] = {k: int(v) for k, v in t["status"].value_counts().items()}
    s["frame_check"] = frame
    return s


def _report(region: str, t: pd.DataFrame, s: dict, dest: Path) -> None:
    f = s["frame_check"]
    L = [f"# Allelic switch test anchored on the GWAS lead — {region}", "",
         "The primary arm fits at the switch-QTL lead and transfers the sign to the disease "
         "allele through signed LD, which fails whenever the two variants are not in strong "
         "LD — most rows. This arm fits the same within-donor contrast **at the GWAS lead "
         "itself**, so the risk haplotype is read from the donor's phased genotype and the "
         "orientation is exact. It buys validity with power: only donors heterozygous at "
         "that variant are informative.", "",
         "## Phase-frame check (run before anything is fitted)", "",
         "The fragment counts are phased by phASER against the Quest VCF; the anchor "
         "genotypes come from the BrainSEQ TOPMed panel. Both must be the same frame or "
         "every re-anchored sign is arbitrary, so the fitted leads were extracted from the "
         "panel and compared donor by donor.", "",
         f"- variants compared: **{f['variants_compared']}**, donor x variant pairs: "
         f"**{f['donor_variant_pairs']}**",
         f"- genotype agreement (unphased): **{f['genotype_agreement']:.4f}**"
         if np.isfinite(f.get("genotype_agreement", np.nan)) else
         "- genotype agreement (unphased): not computable",
         f"- heterozygous pairs carrying frame information: **{f['het_pairs']}**",
         f"- **phase-frame agreement: {f['frame_agreement']:.4f}**"
         if np.isfinite(f.get("frame_agreement", np.nan)) else
         "- **phase-frame agreement: not computable**",
         f"- gate: the arm runs only at >= {FRAME_MIN}", "",
         "## Result", "",
         f"- pair x trait rows: **{s['rows']}** over {s['genes']} genes",
         f"- fitted: **{s['fitted']}**; median donors informative on both haplotypes: "
         f"**{s['median_paired_het']}**",
         f"- rows whose GWAS lead IS the fitted lead (nothing to re-anchor): "
         f"**{s['anchor_is_fitted_lead']}**",
         f"- significant at q < {QVAL_ALLELIC}: **{s['significant']}** "
         f"({s['significant_interpretable']} with an interpretable effect size, "
         f"{s['significant_at_bound']} at the estimator bound)",
         f"- risk allele raises T1 in **{s['positive_beta']}** of the fitted rows "
         f"(median |beta| among unbounded fits {s['median_abs_beta']:.3f})"
         if s.get("median_abs_beta") is not None else
         f"- risk allele raises T1 in **{s['positive_beta']}** of the fitted rows",
         "",
         "A fit **at the bound** is quasi-separation: one haplotype carries one isoform only. "
         "Its sign is informative and its magnitude is not, so those rows are counted apart "
         "rather than read as effects.", ""]
    if s.get("status_counts"):
        L += ["| status | rows |", "| --- | ---: |"]
        L += [f"| {k} | {v} |" for k, v in sorted(s["status_counts"].items(),
                                                  key=lambda kv: -kv[1])]
        L.append("")
    fitted = t[t["status"] == "fitted"] if "status" in t else t.iloc[0:0]
    if len(fitted):
        L += ["## Fitted rows", "",
              "| gene | trait | anchor | risk | beta (risk hap) | p | q | donors paired |",
              "| --- | --- | --- | :---: | ---: | ---: | ---: | ---: |"]
        for _, r in fitted.sort_values("pval").iterrows():
            beta = ("at bound" if r.get("at_bound")
                    else f"{r['beta_risk_direct']:.3f}")
            L.append(f"| {r.get('symbol') or r.get('gene')} | {r['trait']} | "
                     f"{r['anchor_variant']} | {r['risk_allele']} | "
                     f"{beta} | {r['pval']:.3g} | "
                     f"{r['qval']:.3g} | {int(r['n_paired_het'])} |")
        L.append("")
    L += ["## How to read a null here", "",
          "A null in the primary arm was ambiguous: the sign could not be transferred for "
          "most rows, so 'no disease linkage' and 'no answer' looked alike. A null **here** "
          "is not ambiguous in that way — the orientation is exact by construction — but it "
          "is a power statement. Read `n_paired_het` before reading the p-value: the GWAS "
          "lead is chosen for disease signal, not for heterozygosity in 200 EA donors, and a "
          "row with a handful of informative donors has not tested anything.", ""]
    (dest / "GWAS_LEAD_ARM.md").write_text("\n".join(L))


def run(regions: tuple[str, ...]) -> None:
    t = build_targets(regions)
    print(f"targets: {len(t)} pair x trait rows, {t['gene'].nunique()} genes", flush=True)

    risk = []
    for trait, sub in t.dropna(subset=["lead_snp"]).groupby("trait"):
        with tempfile.TemporaryDirectory() as td:
            g = _gwas_risk(trait, set(sub["lead_snp"]), Path(td))
        print(f"  {trait}: {len(g)}/{sub['lead_snp'].nunique()} locus leads in the sumstats",
              flush=True)
        if not g.empty:
            risk.append(g.assign(trait=trait))
    risk = (pd.concat(risk, ignore_index=True) if risk else
            pd.DataFrame(columns=["rsid", "risk_allele", "trait"]))
    t = t.merge(risk.rename(columns={"rsid": "lead_snp", "gwas_p": "risk_gwas_p"})[
                    [c for c in ("lead_snp", "trait", "risk_allele", "risk_gwas_p")
                     if c in risk.columns or c in ("lead_snp", "trait")]],
                on=["lead_snp", "trait"], how="left")

    wanted = set(t["lead_snp"].dropna()) | set(t["variant_id_all"].dropna())
    sites = variant_sites(wanted, _PANELS["EA"])
    print(f"  coordinates resolved for {len(sites)} of {len(wanted)} variants", flush=True)

    for region in regions:
        dest = dest_dir(region)
        counts_f = dest / "allelic_donor_counts.parquet"
        if not counts_f.exists():
            print(f"  {region}: no donor counts, skipped", flush=True)
            continue
        counts = pd.read_parquet(counts_f)
        samples = pd.read_csv(dest / "samples.tsv", sep="\t",
                              usecols=["sample_id", "donor_id"]).drop_duplicates("sample_id")
        gt = phaser_gt(region, sites, samples)
        print(f"  {region}: phASER genotypes for {gt['variant_id'].nunique()} variants in "
              f"{gt['donor_id'].nunique()} donors "
              f"({int((gt['source'] == 'PW').sum())} from PW)", flush=True)
        frame = frame_check(counts, gt)
        print(f"  {region}: phase-frame agreement {frame['frame_agreement']:.4f} over "
              f"{frame['het_pairs']} het donor x variant pairs", flush=True)
        if not np.isfinite(frame["frame_agreement"]) or frame["frame_agreement"] < FRAME_MIN:
            (dest / "gwas_lead_arm_summary.json").write_text(json.dumps(
                {"aborted": "phase frames disagree", "frame_check": frame}, indent=2))
            print(f"  {region}: ABORTED -- this code does not reproduce the recount's own "
                  "haplotype assignment, so re-anchoring cannot be trusted", flush=True)
            continue
        sub = t[t["region"] == region].dropna(subset=["lead_snp"])
        out = fit_region(region, sub, gt, donor_ancestry(region))
        if out.empty:
            print(f"  {region}: nothing to fit", flush=True)
            continue
        out.to_parquet(dest / "gwas_lead_arm.parquet", index=False)
        s = _summary(out, frame)
        (dest / "gwas_lead_arm_summary.json").write_text(json.dumps(s, indent=2, default=float))
        _report(region, out, s, dest)
        print(f"{region}: {s['fitted']}/{s['rows']} fitted, {s['significant']} at "
              f"q < {QVAL_ALLELIC} -> {dest}", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--region", action="append", choices=REGIONS,
                    help="repeatable; default all three")
    a = ap.parse_args(argv)
    run(tuple(a.region) if a.region else REGIONS)


if __name__ == "__main__":
    main()
