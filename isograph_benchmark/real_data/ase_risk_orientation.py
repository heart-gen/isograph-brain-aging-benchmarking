"""Orient the allele-specific switch effect to the GWAS risk allele (PI item 10a, step 5).

WHAT THIS ADDS
--------------
`ase_junction_allelic` fits, per switch pair, the within-donor log odds of isoform T1 on the
haplotype carrying the switch-QTL lead's ALT allele (`beta`). ALT is a sequencing convention,
so the sign says nothing about disease until it is tied to the allele that RAISES disease
risk. This stage does that tie, and only that: it is a sign flip plus its justification. The
GLMM is not refitted, which is why it runs on Bridges-2, where the GWAS sumstats live, and
not on Quest.

WHICH VARIANT CARRIES THE DISEASE EVIDENCE
------------------------------------------
Not the fitted lead. The lead is the BrainSEQ all_samples switch-QTL lead: chosen for its QTL
signal, and mostly unremarkable in the GWAS (median GWAS p ~ 0.2 over the nominated leads, 5
of 60 genome-wide significant). Reading a "risk allele" off a beta that is indistinguishable
from noise would assign signs at random.

The disease evidence sits at the locus the gene colocalizes in. `candidate_loci.tsv` defines
those loci and names each one's GWAS lead (`lead_snp`, `lead_p`); the loci are built around
genome-wide-significant signals, so this variant always carries real disease evidence. The
locus is the one where this gene x trait cell colocalizes at the sQTL layer (PP4 >= 0.8, the
nomination threshold), taken from the signal-level hierarchy.

ANCESTRY: WHY THE FIT AND THE LD ARE BOTH EA-ONLY HERE
------------------------------------------------------
BrainSEQ is roughly half African-American (44-49% of the tested donors per region); the GWAS
are European. Both halves of the orientation break on that.

  * LD does not transfer. Which lead allele carries the risk allele is a property of a
    haplotype background, and it differs between EA and AA donors. Worse, LD measured in the
    pooled panel is partly admixture LD -- correlation created by population structure rather
    than by proximity -- so a pooled estimate is not "the average of the two", it is a third
    thing that describes neither.
  * The pooled `beta` is a mixture. The Quest fit pools EA and AA donors, which is right for
    "is this switch under cis control" (phase is resolved within each donor) but not for
    "does the EUR risk allele push it", because the EA-derived pairing may be wrong in the AA
    donors.

So this stage refits the GLMM in EA donors only, from the donor-level counts the recount
already shipped (`allelic_donor_counts.parquet`), and takes the LD from the EA panel that the
ea_only QTL arm uses. `beta_ea` is what gets oriented. The AA arm is fitted and reported
beside it -- never pooled, never used to orient against a EUR GWAS -- because a sign that
agrees in both ancestries is worth more than one seen in either alone.

WHY SIGNED LD, AND WHY IN-SAMPLE
--------------------------------
The risk allele and the fitted lead are two different variants, so the flip is only defined
once we know which lead allele TRAVELS with the risk allele on a haplotype. That is the sign
of D (equivalently of r), not r^2 and not distance. plink2's `--ld` prints the two-locus
haplotype frequency table with its allele labels, which determines the pairing directly.

The panel is the BrainSEQ TOPMed genotypes the switch-QTL mapping used -- the same donors the
allelic test is fitted in, so the LD is the LD those haplotypes actually have, not a
population average. Where the two variants are in weak LD the pairing is not determined and
the pair is left unoriented (`r_gate` reason), because a flip decided at |r| ~ 0.3 is a coin
toss dressed as a result.

OUTPUT
------
`risk_orientation.parquet` per region, plus `ASE_RISK_ORIENTATION.md`. Key columns:

  risk_allele, risk_gwas_p, risk_beta_gwas  the locus GWAS lead and its trait-increasing allele
  r_lead_risk                               signed r between lead ALT and the risk allele (EA panel)
  r_lead_risk_aa                            the same in the AA panel, for disclosure
  beta_ea, pval_ea, n_paired_ea             the GLMM refitted in EA donors only
  beta_aa, pval_aa, n_paired_aa             the same in AA donors, never pooled
  beta_risk                                 beta_ea oriented to the risk allele (NA when ungated)
  risk_along_module                         beta_risk signed by the pair's module polarity:
                                            > 0 when the risk allele shifts the isoforms the
                                            way the module score rises

`risk_along_module` is the disease-linkage readout the allelic table could not carry.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel, stage_out
from isograph_benchmark.real_data.ase_junction_allelic import (
    QVAL_ALLELIC, module_annotation as _module_annotation, out_dir)
from isograph_benchmark.real_data.coloc_direction import _gwas_risk

REGIONS: tuple[str, ...] = ("caudate", "dlpfc", "hippocampus")

# The nomination threshold of the signal-level layer; the same cut the locus event audit uses.
PP4_MIN = 0.8
# |r| below this leaves the pair unoriented. The flip is a hard call, not a graded one, so the
# bar is the one used for "same signal" elsewhere in the project rather than a nominal r != 0.
R_MIN = 0.8

_PLINK2 = Path("/ocean/projects/bio260021p/shared/opt/envs/eqtl/bin/plink2")
_GENO = Path("/ocean/projects/bio260021p/shared/resources/processed-data/genotypes/qtl")
# EA is the GWAS-matched panel and the one the ea_only QTL arm uses; AA is reported beside it.
_PANELS = {"EA": _GENO / "ea_only" / "TOPMed_LIBD.EA-merge",
           "AA": _GENO / "aa_only" / "TOPMed_LIBD.AA-merge"}
_PANEL = _PANELS["EA"]
# Fallback only: a GWAS lead the QTL panel does not carry (imputation/MAF filtered) still has
# LD in the reference panel the coloc pipeline uses. Population LD, not these donors', so it
# is recorded as such in `ld_panel`.
_PANEL_1KG = Path("/ocean/projects/bio250020p/shared/resources/ldsc/1000G_EUR_Phase3_plink")


def dest_dir(region: str) -> Path:
    return ensure_dir(out_dir(region))


# --------------------------------------------------------------------------- #
# Targets: nominated pair -> colocalizing locus -> that locus's GWAS lead
# --------------------------------------------------------------------------- #
def _nominated_cells() -> pd.DataFrame:
    """sQTL cells at or above the nomination threshold, best cell per gene x trait."""
    cells = pd.read_parquet(
        stage_out("anchoring", "coloc_signal_susie", "cells_hierarchy.parquet"),
        columns=["trait", "gene", "tissue", "modality", "PP4", "LOCUS_ID", "symbol",
                 "estimator"])
    cells = cells[(cells["modality"] == "sQTL") & (cells["PP4"] >= PP4_MIN)]
    best = (cells.sort_values("PP4", ascending=False)
                 .drop_duplicates(["gene", "trait"])
                 .rename(columns={"tissue": "coloc_tissue", "PP4": "coloc_PP4",
                                  "estimator": "coloc_estimator"}))
    # A pair can be nominated by a layer that is not in this table (the abf arm of the audit,
    # or the BrainSEQ coloc). Those carry a locus too, and the locus is what we need.
    extra = []
    for f in (stage_out("anchoring", "locus_event_audit", "abf", "locus_event_audit.parquet"),
              stage_out("anchoring", "locus_event_audit", "susie_all_introns",
                        "locus_event_audit.parquet"),
              stage_out("anchoring", "coloc_brainseq", "ea_only", "nominations.parquet")):
        if f.exists():
            extra.append(pd.read_parquet(f, columns=["trait", "gene", "LOCUS_ID", "symbol"]))
    if extra:
        e = pd.concat(extra, ignore_index=True).drop_duplicates(["gene", "trait"])
        e = e[~e.set_index(["gene", "trait"]).index.isin(
            best.set_index(["gene", "trait"]).index)]
        best = pd.concat([best, e], ignore_index=True)
    return best


def _locus_leads() -> pd.DataFrame:
    """LOCUS_ID -> GWAS lead variant, per trait, from the coloc locus definitions."""
    rows = []
    for d in sorted((stage_out("anchoring", "coloc")).glob("*__*")):
        f = d / "candidate_loci.tsv"
        if not f.exists() or "__" not in d.name:
            continue
        t = pd.read_csv(f, sep="\t", usecols=["LOCUS_ID", "lead_snp", "lead_p"])
        # both the aging and the BrainSEQ-SCZD analyses define loci; a locus id is unique
        # within a trait, and the disease lead is a property of the trait, not the analysis
        rows.append(t.assign(trait=d.name.split("__")[1]))
    if not rows:
        raise SystemExit("no candidate_loci.tsv found under the coloc outputs")
    return pd.concat(rows, ignore_index=True).drop_duplicates(["LOCUS_ID", "trait"])


def build_targets(regions: tuple[str, ...]) -> pd.DataFrame:
    """One row per (region, pair, trait) for the coloc-nominated pairs."""
    cells, loci = _nominated_cells(), _locus_leads()
    rows = []
    for region in regions:
        f = dest_dir(region) / "allelic_test.parquet"
        if not f.exists():
            print(f"  {region}: no allelic_test.parquet, skipped", flush=True)
            continue
        a = pd.read_parquet(f)
        # the Quest table carries empty placeholders for these; this stage fills them
        a = a.drop(columns=["risk_allele", "beta_risk", "risk_along_module"],
                   errors="ignore")
        a = a[a["coloc_nominated"] & a["coloc_traits"].notna()].copy()
        a["gene"] = a["gene_id"].str.split(".").str[0]
        a["region"] = region
        # coloc_traits is a comma-separated set: one orientation per trait.
        a["trait"] = a["coloc_traits"].str.split(",")
        rows.append(a.explode("trait").assign(trait=lambda d: d["trait"].str.strip()))
    if not rows:
        raise SystemExit("no allelic tables to orient")
    t = pd.concat(rows, ignore_index=True)
    t = t.merge(cells, on=["gene", "trait"], how="left").merge(
        loci, on=["LOCUS_ID", "trait"], how="left")
    return t


# --------------------------------------------------------------------------- #
# Ancestry: label the donors, refit each arm on its own
# --------------------------------------------------------------------------- #
def _panel_donors(panel: Path) -> set[str]:
    # the panel stem contains dots (TOPMed_LIBD.EA-merge), so with_suffix would eat them
    psam = Path(str(panel) + ".psam")
    return {ln.split()[0] for ln in psam.read_text().splitlines()
            if ln and not ln.startswith("#")}


def donor_ancestry(region: str) -> pd.Series:
    """donor_id -> EA / AA, from membership of the two QTL genotype panels."""
    s = pd.read_csv(dest_dir(region) / "samples.tsv", sep="\t")
    gid = s.set_index("donor_id")["genotype_sample_id"].astype(str)
    gid = gid[~gid.index.duplicated()]
    ea, aa = _panel_donors(_PANELS["EA"]), _panel_donors(_PANELS["AA"])
    return gid.map(lambda g: "EA" if g in ea else ("AA" if g in aa else None))


def refit_by_ancestry(region: str, pair_ids: set[str]) -> pd.DataFrame:
    """Refit the pair GLMM within each ancestry, from the shipped donor-level counts.

    Same estimator as the Quest step (`allelic_contrast`), only the donor set changes, so
    `beta_ea` is comparable with the pooled `beta` and is the one a European GWAS can orient.
    """
    from isograph_benchmark.real_data.ase_junction_allelic import allelic_contrast

    counts = pd.read_parquet(dest_dir(region) / "allelic_donor_counts.parquet")
    counts = counts[counts["pair_id"].isin(pair_ids)]
    anc = donor_ancestry(region)
    counts = counts.assign(ancestry=counts["donor_id"].map(anc))
    rows = []
    for pair_id, sub in counts.groupby("pair_id"):
        rec = {"pair_id": pair_id}
        for arm in ("EA", "AA"):
            tab, k = sub[sub["ancestry"] == arm], arm.lower()
            if tab.empty:
                rec |= {f"beta_{k}": np.nan, f"pval_{k}": np.nan, f"n_paired_{k}": 0,
                        f"status_{k}": "no_donors"}
                continue
            f = allelic_contrast(tab)
            rec |= {f"beta_{k}": f.get("beta", np.nan), f"pval_{k}": f.get("pval", np.nan),
                    f"n_paired_{k}": f.get("n_paired_het", 0),
                    f"status_{k}": f.get("status", "not_fitted"),
                    # a bound fit is quasi-separation: its sign is informative, its
                    # magnitude is not, so it is flagged rather than reported as an effect
                    f"at_bound_{k}": bool(f.get("at_bound", False))}
        rows.append(rec)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    both = out["beta_ea"].notna() & out["beta_aa"].notna()
    out["sign_agrees_ancestry"] = pd.Series(pd.NA, index=out.index, dtype="boolean")
    out.loc[both, "sign_agrees_ancestry"] = (np.sign(out.loc[both, "beta_ea"])
                                             == np.sign(out.loc[both, "beta_aa"]))
    return out


# --------------------------------------------------------------------------- #
# Signed LD between the fitted lead and the locus GWAS lead
# --------------------------------------------------------------------------- #
def _subset_panel(rsids: set[str], tmp: Path, panel: Path, tag: str) -> Path:
    """One pass over the panel; every --ld call then reads a small pfile."""
    keep = tmp / "keep.txt"
    keep.write_text("\n".join(sorted(rsids)) + "\n")
    out = tmp / f"panel_{tag}"
    subprocess.run([str(_PLINK2), "--pfile", str(panel), "--extract", str(keep),
                    "--make-pgen", "--out", str(out)],
                   check=True, capture_output=True, text=True)
    return out


def _run_ld(pfile: list[str], a: str, b: str, tmp: Path) -> dict | None:
    p = subprocess.run([str(_PLINK2), *pfile, "--ld", a, b, "--out", str(tmp / "ld")],
                       capture_output=True, text=True)
    log = tmp / "ld.log"
    return _parse_ld(log.read_text() if log.exists() else p.stdout + p.stderr)


def _ld_1kg(a: str, b: str, chrom: str | None, tmp: Path) -> dict | None:
    """Same two variants in 1000G EUR, for a variant the QTL panel does not carry."""
    if not chrom:
        return None
    bfile = _PANEL_1KG / f"1000G.EUR.QC.{chrom}"
    if not Path(str(bfile) + ".bed").exists():
        return None
    return _run_ld(["--bfile", str(bfile)], a, b, tmp)


def _parse_ld(log_text: str) -> dict | None:
    """Pull the allele labels and the 2x2 haplotype frequency table out of `--ld` output.

    plink2 prints MAJOR/MINOR per variant and the four haplotype frequencies. That table
    fixes which allele of one variant sits with which allele of the other, which r^2 alone
    cannot.
    """
    lines = [ln.rstrip() for ln in log_text.splitlines()]
    alleles: list[dict[str, str]] = []
    freqs: list[float] = []
    block: dict[str, str] | None = None
    for ln in lines:
        s = ln.strip()
        # "<rsid> alleles:" opens a block of "MAJOR = REF = C" / "MINOR = T" / "(REF = C)"
        if s.endswith("alleles:"):
            if block and {"major", "minor"} <= set(block):
                alleles.append(block)
            block = {}
            continue
        if block is not None and s and not s.startswith("("):
            if s.startswith("MAJOR") and "0." not in s:
                block.setdefault("major", s.split("=")[-1].strip())
                continue
            if s.startswith("MINOR") and "0." not in s:
                block.setdefault("minor", s.split("=")[-1].strip())
                continue
            if not s.startswith(("MAJOR", "MINOR")):
                if {"major", "minor"} <= set(block):
                    alleles.append(block)
                block = None
        # the haplotype frequency rows carry two numbers each, MAJOR then MINOR
        if s.startswith("MAJOR") and "0." in s:
            freqs += [float(x) for x in s.split()[1:3]]
        if s.startswith("MINOR") and "0." in s:
            freqs += [float(x) for x in s.split()[1:3]]
    if block and {"major", "minor"} <= set(block):
        alleles.append(block)
    if len(alleles) != 2 or len(freqs) != 4:
        return None
    # freqs order: (A major, B major), (A major, B minor), (A minor, B major), (A minor, B minor)
    return {"a": alleles[0], "b": alleles[1],
            "f": {("major", "major"): freqs[0], ("major", "minor"): freqs[1],
                  ("minor", "major"): freqs[2], ("minor", "minor"): freqs[3]}}


def signed_r(ld: dict, allele_a: str, allele_b: str) -> float | None:
    """Signed r for the specific allele pair, from the haplotype frequency table."""
    lab_a = "major" if allele_a == ld["a"]["major"] else (
        "minor" if allele_a == ld["a"]["minor"] else None)
    lab_b = "major" if allele_b == ld["b"]["major"] else (
        "minor" if allele_b == ld["b"]["minor"] else None)
    if lab_a is None or lab_b is None:
        return None
    f = ld["f"]
    p_a = f[("major", "major")] + f[("major", "minor")] if lab_a == "major" else \
        f[("minor", "major")] + f[("minor", "minor")]
    p_b = f[("major", "major")] + f[("minor", "major")] if lab_b == "major" else \
        f[("major", "minor")] + f[("minor", "minor")]
    d = f[(lab_a, lab_b)] - p_a * p_b
    denom = p_a * (1 - p_a) * p_b * (1 - p_b)
    return float(d / np.sqrt(denom)) if denom > 0 else None


def ld_table(pairs: pd.DataFrame, tmp: Path) -> pd.DataFrame:
    """Signed r between (lead ALT) and (risk allele) for each distinct variant pair."""
    rsids = set(pairs["variant_id_all"]) | set(pairs["lead_snp"])
    panels = {tag: _subset_panel({r for r in rsids if isinstance(r, str)}, tmp, path, tag)
              for tag, path in _PANELS.items()}
    rows = []
    for (lead, gwas), sub in pairs.groupby(["variant_id_all", "lead_snp"]):
        rec = {"variant_id_all": lead, "lead_snp": gwas, "r_lead_risk": np.nan,
               "r_lead_risk_aa": np.nan, "ld_status": "ok", "ld_panel": "brainseq_ea"}
        if lead == gwas:
            # the fitted lead IS the locus lead: orientation is direct, r = 1 by definition
            rec["r_lead_risk"] = rec["r_lead_risk_aa"] = 1.0
            rec["ld_status"] = "same_variant"
            rows.append(rec)
            continue
        alt, risk = sub["lead_alt"].iloc[0], sub["risk_allele"].iloc[0]
        for tag, pf in panels.items():
            ld = _run_ld(["--pfile", str(pf)], lead, gwas, tmp)
            if ld is None and tag == "EA":
                # the QTL panels are MAF/imputation filtered; the EUR reference can still
                # carry the pair, and it is the ancestry the GWAS was run in
                # LOCUS_ID is "<locus>_chr<N>", the one chromosome label every source carries
                loc = str(sub["LOCUS_ID"].iloc[0])
                ld = _ld_1kg(lead, gwas, loc.split("_chr")[-1] if "_chr" in loc else None,
                             tmp)
                rec["ld_panel"] = "1000G_EUR" if ld is not None else "none"
            if ld is None:
                continue
            r = (signed_r(ld, str(alt).upper(), str(risk).upper())
                 if pd.notna(risk) else None)
            if r is not None:
                rec["r_lead_risk" if tag == "EA" else "r_lead_risk_aa"] = r
            elif tag == "EA":
                rec["ld_status"] = "allele_mismatch"
        if np.isnan(rec["r_lead_risk"]) and rec["ld_status"] == "ok":
            rec["ld_status"] = "not_in_panel"
        rows.append(rec)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Orientation
# --------------------------------------------------------------------------- #
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _palindromic(a: pd.Series, b: pd.Series) -> pd.Series:
    """A/T or C/G GWAS leads: the allele labels match on either strand.

    For any other variant a strand flip makes the risk allele fail to match the panel's
    alleles and the pair drops out as `allele_mismatch`. A palindromic one matches either
    way, so a strand disagreement between the sumstats and the panel would flip the sign
    silently. These are left unoriented instead.
    """
    a, b = a.astype("string").str.upper(), b.astype("string").str.upper()
    return pd.Series([x in _COMPLEMENT and _COMPLEMENT.get(x) == y
                      for x, y in zip(a.fillna(""), b.fillna(""))], index=a.index)


def projection_guard(t: pd.DataFrame) -> pd.DataFrame:
    """Decide, per pair, whether the module switch axis means anything for it.

    The pair's polarity is a difference of two transcript-eigengene correlations. It is only
    an axis worth signing an allelic effect along when (a) the gene is in its module through
    SWITCHING rather than abundance, and (b) at least one of the two transcripts actually
    correlates with the module score. KLC1 failed both and produced four significant rows whose
    `risk_along_module` signs disagreed while the underlying cis effect was a single, coherent
    one (stage 06a, section 5a). The upstream Quest step now writes
    `module_projection_status`; this recomputes it where an older table lacks it, so the guard
    holds without waiting for a re-run on Quest.
    """
    from isograph_benchmark.real_data.ase_junction_allelic import (
        POLARITY_Q, _PROJECTABLE_ROLES)

    if "module_projection_status" in t.columns:
        return t
    needed = {"module_id", "module_role", "module_polarity",
              "transcript_id_1", "transcript_id_2"}
    if not needed.issubset(t.columns):
        # A caller that does not carry the module annotation cannot be judged; say so rather
        # than silently guarding or silently projecting.
        t["module_projection_status"] = "not_evaluated"
        return t
    q_lookup: dict[str, object] = {}
    for region in sorted(set(t.get("region", pd.Series(dtype=str)).dropna())):
        try:
            _, polarity, _ = _module_annotation(region)
        except Exception as exc:                     # noqa: BLE001 - reported, not raised
            print(f"  {region}: polarity q unavailable ({exc}); projection guard falls back "
                  "to the role test only", flush=True)
            polarity = None
        q_lookup[region] = polarity

    def _q(region, module, tx):
        pol = q_lookup.get(region)
        if pol is None or not isinstance(module, str):
            return np.nan
        try:
            return float(pol.loc[(module, tx), "qvalue"])
        except (KeyError, TypeError):
            return np.nan

    q1 = np.array([_q(r, m, a) for r, m, a in
                   zip(t.get("region", pd.Series(index=t.index)), t["module_id"],
                       t["transcript_id_1"])], dtype=float)
    q2 = np.array([_q(r, m, b) for r, m, b in
                   zip(t.get("region", pd.Series(index=t.index)), t["module_id"],
                       t["transcript_id_2"])], dtype=float)
    anchored = (np.nan_to_num(q1, nan=1.0) <= POLARITY_Q) | \
               (np.nan_to_num(q2, nan=1.0) <= POLARITY_Q)
    t["pair_polarity_anchored"] = anchored
    role_ok = t["module_role"].isin(_PROJECTABLE_ROLES)
    t["module_projection_status"] = np.select(
        [t["module_id"].isna(), ~role_ok, ~anchored,
         t["module_polarity"].isna() | (t["module_polarity"] == 0)],
        ["no_module", "role_not_switching", "pair_polarity_not_anchored", "no_polarity"],
        default="projected")
    return t


def orient(t: pd.DataFrame, beta_col: str = "beta") -> pd.DataFrame:
    """beta -> beta_risk -> risk_along_module, with the reason when it stays NA.

    `beta_col` is the EA-only refit in the pipeline; the pooled `beta` mixes ancestries and
    cannot be oriented against a European GWAS (see the module docstring).
    """
    b = t[beta_col]
    on_alt = t["r_lead_risk"] > 0
    gated = t["r_lead_risk"].abs() >= R_MIN
    other = t["gwas_other_allele"] if "gwas_other_allele" in t else pd.Series(
        pd.NA, index=t.index, dtype="string")
    pal = _palindromic(t["risk_allele"], other)
    t["gwas_lead_palindromic"] = pal
    ok = gated & ~pal
    t["beta_risk"] = np.where(ok, b * np.where(on_alt, 1.0, -1.0), np.nan)
    t = projection_guard(t)
    s = np.sign(t["module_polarity"]).replace(0, np.nan).where(
        t["module_projection_status"].isin(("projected", "not_evaluated")))
    t["risk_along_module"] = t["beta_risk"] * s
    t["orientation_status"] = np.select(
        [t["risk_allele"].isna(), t["r_lead_risk"].isna(), pal, ~gated, b.isna()],
        ["no_gwas_risk_allele", "no_ld", "palindromic_gwas_lead", "r_gate",
         "not_fitted_in_ea"], default="oriented")
    return t


def _bh(p: pd.Series) -> pd.Series:
    p = p.astype(float)
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.clip(q, 0, 1)
    return pd.Series(out, index=p.index)


def _summary(t: pd.DataFrame, region: str) -> dict:
    o = t[t["orientation_status"] == "oriented"]
    sig = o[(o["qval_ea"] < QVAL_ALLELIC) & ~o["at_bound_ea"].eq(True)]
    vp = t.dropna(subset=["r_lead_risk"]).drop_duplicates(["variant_id_all", "lead_snp"])
    vp_aa = t.dropna(subset=["r_lead_risk_aa"]).drop_duplicates(["variant_id_all",
                                                                 "lead_snp"])
    anc = donor_ancestry(region)
    # sign agreement is a property of the pair, not of the pair x trait rows
    pair = t.drop_duplicates("pair_id")
    cmp_anc = pair["sign_agrees_ancestry"].notna() if "sign_agrees_ancestry" in pair else \
        pd.Series(False, index=pair.index)
    return {
        "donors_ea": int((anc == "EA").sum()),
        "donors_aa": int((anc == "AA").sum()),
        "donors_other": int(anc.isna().sum()),
        "median_abs_r_aa": (float(vp_aa["r_lead_risk_aa"].abs().median())
                            if len(vp_aa) else None),
        "sign_compared_ancestry": int(cmp_anc.sum()),
        "sign_agrees_ancestry": int(pair.loc[cmp_anc, "sign_agrees_ancestry"].eq(True).sum()),
        "variant_pairs_with_ld": int(len(vp)),
        "median_abs_r": float(vp["r_lead_risk"].abs().median()) if len(vp) else None,
        "variant_pairs_at_r_gate": int((vp["r_lead_risk"].abs() >= R_MIN).sum()),
        "pairs": int(len(t)), "genes": int(t["gene"].nunique()),
        "fitted": int((t["status"] == "fitted").sum()),
        "oriented": int(len(o)), "oriented_genes": int(o["gene"].nunique()),
        "status_counts": t["orientation_status"].value_counts().to_dict(),
        "significant_oriented": int(len(sig)),
        "significant_genes": int(sig["gene"].nunique()),
        "risk_raises_module_score": int((sig["risk_along_module"] > 0).sum()),
        "risk_lowers_module_score": int((sig["risk_along_module"] < 0).sum()),
    }


def _report(region: str, t: pd.DataFrame, s: dict, dest: Path) -> None:
    L = [f"# Risk-allele orientation of the allelic switch test — {region}", "",
         "`beta` from the allelic test is the within-donor log odds of T1 on the haplotype "
         "carrying the switch-QTL lead's ALT allele. ALT is arbitrary, so this stage re-signs "
         "it to the allele that raises disease risk. The risk allele is the trait-increasing "
         "allele at the GWAS lead of the locus where the gene colocalizes (sQTL PP4 >= "
         f"{PP4_MIN}); the flip follows the signed LD between that variant and the fitted "
         f"lead, in the BrainSEQ genotypes the QTL mapping used, gated at |r| >= {R_MIN}.", "",
         "`risk_along_module` > 0 means the risk allele shifts the isoforms the way the "
         "module score rises.", "",
         "## Coverage", "",
         f"- nominated pair x trait rows: **{s['pairs']}** over {s['genes']} genes "
         f"({s['fitted']} fitted)",
         f"- oriented: **{s['oriented']}** ({s['oriented_genes']} genes)", ""]
    L += ["| orientation status | rows |", "|---|---:|"]
    L += [f"| {k} | {v} |" for k, v in sorted(s["status_counts"].items(),
                                              key=lambda kv: -kv[1])]
    L += ["", "## LD between the fitted lead and the disease lead", "",
          "This is what limits the arm. The switch-QTL lead the test is fitted on and the "
          "GWAS lead of the locus are usually not on the same haplotype, so the disease "
          "allele cannot be placed on the fitted lead at all.", "",
          f"- distinct lead x disease-lead variant pairs with LD: **{s['variant_pairs_with_ld']}**",
          f"- median |r|: **{s['median_abs_r']:.3f}**" if s["median_abs_r"] is not None
          else "- median |r|: n/a",
          f"- at |r| >= {R_MIN}: **{s['variant_pairs_at_r_gate']}**", "",
          "## Ancestry", "",
          "BrainSEQ is roughly half African-American and the GWAS are European, so a pooled "
          "fit and a pooled LD estimate both mix ancestries. `beta_ea` (EA donors only) is "
          "what gets oriented, against EA-panel LD; the AA arm is fitted beside it and "
          "never pooled.", "",
          f"- donors: EA {s['donors_ea']}, AA {s['donors_aa']}, neither {s['donors_other']}",
          "- median |r| lead-to-risk-allele: EA "
          + (f"{s['median_abs_r']:.3f}" if s["median_abs_r"] is not None else "n/a")
          + (f", AA {s['median_abs_r_aa']:.3f}" if s.get("median_abs_r_aa") is not None
             else ""),
          f"- EA and AA betas agree in sign in {s['sign_agrees_ancestry']} of "
          f"{s['sign_compared_ancestry']} pairs fitted in both", "",
          "## Oriented pairs at q < 0.05 (EA refit, BH within this family)", "",
          f"- {s['significant_oriented']} pairs over {s['significant_genes']} genes",
          f"- risk allele raises the module score in {s['risk_raises_module_score']}, "
          f"lowers it in {s['risk_lowers_module_score']}", ""]
    o = t[(t["orientation_status"] == "oriented") & (t["qval_ea"] < QVAL_ALLELIC)
          & ~t["at_bound_ea"].eq(True)]
    if not o.empty:
        L += ["| gene | trait | lead | risk variant | risk allele | GWAS p | r | beta_risk | "
              "risk_along_module | q |", "|---|---|---|---|---|---:|---:|---:|---:|---:|"]
        for _, r in o.sort_values("qval_ea").head(25).iterrows():
            L.append(f"| {r.get('symbol') or r['gene']} | {r['trait']} | "
                     f"{r['variant_id_all']} | {r['lead_snp']} | {r['risk_allele']} | "
                     f"{r['risk_gwas_p']:.2g} | {r['r_lead_risk']:+.2f} | "
                     f"{r['beta_risk']:+.2f} | {r['risk_along_module']:+.2f} | "
                     f"{r['qval_ea']:.2g} |")
        L.append("")
    L += ["## Caveats", "",
          "- The locus GWAS lead is not necessarily the causal variant; orientation inherits "
          "that assumption, and imperfect LD with the fitted lead attenuates `beta` itself.",
          "- Every caveat of the allelic test carries over (cell-type-specific allelic "
          "effects, residual mapping bias, phase switches).",
          "- Pairs left unoriented are not negative results: the disease allele could not be "
          "placed on a haplotype with the fitted lead.", ""]
    (dest / "ASE_RISK_ORIENTATION.md").write_text("\n".join(L))


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
            pd.DataFrame(columns=["rsid", "risk_allele", "risk_beta", "gwas_p", "trait"]))
    t = t.merge(risk.rename(columns={"rsid": "lead_snp", "risk_beta": "risk_beta_gwas",
                                     "gwas_p": "risk_gwas_p"})[
                    ["lead_snp", "trait", "risk_allele", "gwas_other_allele",
                     "risk_beta_gwas", "risk_gwas_p"]],
                on=["lead_snp", "trait"], how="left")

    have = t.dropna(subset=["lead_snp", "variant_id_all"])
    with tempfile.TemporaryDirectory() as td:
        ld = ld_table(have, Path(td))
    t = t.merge(ld, on=["variant_id_all", "lead_snp"], how="left")

    # the pooled beta mixes EA and AA donors; refit each ancestry on its own
    fits = []
    for region in sorted(t["region"].unique()):
        ids = set(t.loc[t["region"] == region, "pair_id"])
        f = refit_by_ancestry(str(region), ids)
        if not f.empty:
            fits.append(f.assign(region=region))
        print(f"  {region}: refit {len(f)} pairs in EA and AA donors separately", flush=True)
    if fits:
        t = t.merge(pd.concat(fits, ignore_index=True), on=["pair_id", "region"], how="left")
    # the Quest q is BH over the pooled gate-family fits; the EA refits are a different
    # family and need their own correction before anything here is called significant
    t["qval_ea"] = np.nan
    ok = t["pval_ea"].notna()
    uniq = t.loc[ok].drop_duplicates("pair_id")[["pair_id", "pval_ea"]]
    uniq["q"] = _bh(uniq["pval_ea"])
    t.loc[ok, "qval_ea"] = t.loc[ok, "pair_id"].map(uniq.set_index("pair_id")["q"])
    t = orient(t, beta_col="beta_ea")

    for region, sub in t.groupby("region"):
        dest = dest_dir(str(region))
        sub.to_parquet(dest / "risk_orientation.parquet", index=False)
        s = _summary(sub, str(region))
        (dest / "risk_orientation_summary.json").write_text(json.dumps(s, indent=2))
        _report(str(region), sub, s, dest)
        print(f"{region}: {s['oriented']}/{s['pairs']} oriented, "
              f"{s['significant_oriented']} at q < {QVAL_ALLELIC} -> {dest}", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--region", action="append", choices=REGIONS,
                    help="repeatable; default all three")
    a = ap.parse_args(argv)
    run(tuple(a.region) if a.region else REGIONS)


if __name__ == "__main__":
    main()
