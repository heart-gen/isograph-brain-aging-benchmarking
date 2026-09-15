"""SMR + HEIDI on the signal-level colocalization nominations (Analysis 5).

WHAT THIS ADDS, AND HOW FAR ITS NUMBERS MAY BE READ
---------------------------------------------------
`coloc.susie` asks whether a disease association and a QTL share a causal variant. It gives
no effect direction, and it separates linkage from sharing only through the balance of PP3
and PP4. SMR adds a signed, ratio-type estimate relating the genetically predicted molecular
phenotype to disease (`b_SMR`), and HEIDI asks whether the association pattern across the cis
region is compatible with one shared causal variant rather than distinct variants in LD. That
makes it orthogonal corroboration for loci coloc already nominated -- positioned beneath coloc,
never beside it. Three rules govern every sentence written from this module's output:

  * `b_SMR` is a ratio estimate under SMR's assumptions (one causal variant per probe, no
    horizontal pleiotropy). It does not establish causal direction in the biological sense
    and cannot distinguish causality from pleiotropy.
  * Failing to reject HEIDI is not proof of a shared variant. HEIDI is underpowered at GTEx
    brain sample sizes, and more so at BrainSEQ's, so p_HEIDI >= HEIDI_REJECT reads "not
    rejected", never "shared".
  * A HEIDI rejection does not overrule a strong signal-level colocalization, and a significant
    SMR estimate does not promote a locus coloc did not support. Disagreements are reported as
    disagreements, in `agreement`, and neither method adjudicates the other.

SCOPE
-----
GTEx (`--qtl-source gtex`): only the loci the signal-level layer nominated (PP4_sQTL >= 0.8,
all-introns arm), only in the tissues where that sQTL call holds, for every intron phenotype of
the gene and for its eQTL. The colocalizing intron is the pre-specified primary sQTL probe; the
gene's other introns are reported beside it so the SMR evidence is never a maximum chosen over
introns after the fact.

BrainSEQ (`--qtl-source brainseq --arm ea_only`): every GTEx signal-level nomination, plus every
locus BrainSEQ's own colocalization nominates on the switch axis (PP4_S_g >= 0.8 in any region,
`coloc_brainseq`), in all three BrainSEQ regions and on both axes (A_g, S_g). BrainSEQ regions
are not GTEx tissues, so a nomination is tested in every region rather than only where a GTEx
call holds; `gtex_call_in_matched_tissue` records whether the region's tissue-matched GTEx
tissue carries the GTEx sQTL call. Each gene has exactly one S_g and one A_g phenotype, so every
BrainSEQ probe is primary. The coloc posterior an SMR result is compared with is BrainSEQ's own,
for the same region and axis.

Either way this is a scoped follow-up on nominated loci, not a transcriptome-wide scan, and its
multiple-testing family is the probes it actually instrumented.

PRE-SPECIFIED SETTINGS
----------------------
SMR 1.4.2's own defaults, read from its source (src/SMR.cpp) and passed explicitly so a later
release with different defaults cannot change the analysis silently: instrument p < 5e-8, HEIDI
SNPs at p < 1.5654e-3, 3-20 HEIDI SNPs, a 2 Mb cis window, and the allele-frequency QC (0.2 per
SNP, at most 5% failing). Interpretation, fixed before any result existed: SMR significance is
Bonferroni over instrumented probes within (analysis, QTL modality), and p_HEIDI < 0.01 rejects
a single shared causal variant.

SENSITIVITY ARMS
----------------
Two, each writing to its own directory under `sensitivity/` so neither can overwrite the
primary result. The ESD/BESD files depend on neither and are shared, so both run on a
primary array's BESD without rebuilding it.

  --peqtl-smr <p>  A relaxed instrument threshold. One pre-specified arm at 1e-6; this is
                   not a sweep, and the relaxed arm is where weak instruments (see WEAK_F)
                   should be looked for.
  --smr-multi      SMR's multi-SNP test, which combines the cis SNPs surviving LD pruning
                   at r2 LD_MULTI_SNP instead of testing the top SNP alone. It asks whether
                   a signal rests on one lead SNP or on the cis signal more broadly, and is
                   a robustness arm for the SMR estimate -- never a second discovery pass.
                   SMR writes `.msmr` here, not `.smr`, with `p_SMR_multi` added. `p_SMR`
                   and `smr_status` stay the single-SNP quantities so the arms compare row
                   for row; the multi verdict is `smr_multi_status`, corrected over the
                   probes the multi test actually ran on. Where too few cis SNPs survive
                   pruning SMR skips the test, and that probe is `multi_unavailable` -- an
                   absence of evidence, not a null.

INPUTS, AND THE BUILD
---------------------
LD reference: the 1000G EUR Phase 3 PLINK panel (hg19) the coloc layer fine-maps on, for both
sources. For BrainSEQ that is a choice, not a constraint: it keeps the two SMR sources different
in their QTL alone, and HEIDI's LD-based test is better served by the reference panel than by a
region's 169-229 donors. BrainSEQ's in-sample LD is used where it is decisive instead, on the QTL
side of `coloc_brainseq`.
GWAS: the per-locus summary statistics the coloc layer used, with the same long-range-LD
exclusions, written as SMR `.ma`. Its frequency column is NA: the per-locus files carry none,
and with a missing GWAS frequency SMR checks the QTL frequency against the reference panel
alone (`freq_check`, src/SMR_data.cpp) rather than against a number invented here.
QTL, GTEx: v11 all-pairs nominal statistics. GRCh38 variant ids are bridged to rsIDs and each
ESD position is the panel's hg19 position for that rsID, so GWAS, QTL and LD reference agree on
position. The effect allele is the GTEx ALT allele, which is the allele `af` and `slope` refer
to. Strand-ambiguous and panel-mismatched variants are dropped, as in the coloc layer. Probe
positions are hg19 gene TSSs from MAGMA's NCBI37.3 gene.loc.
QTL, BrainSEQ: the ea_only nominal cis statistics from `brainseq_switch_qtl`. Variant ids are
already rsIDs; REF/ALT come from the arm's .pvar and are checked against the panel exactly as for
GTEx, and the effect allele is ALT, which tensorQTL counts (the positive control confirmed it).
S_g slopes are sign-pinned to the discovery fit before they enter the ESD, so b_SMR on S_g points
along the discovery switch axis; a gene the pin cannot orient keeps its mapped sign and is
flagged (`sign_pinned`, `axis_unstable`).

The BrainSEQ source is refused unless the arm is `ea_only` AND both `brainseq_qtl_checks` gates
passed in every region. The mixed-ancestry arm is never accepted: the GWAS and the LD reference
are European.

STAGES
------
  --stage prep   targets, the (tissue-or-region, chr) work list, one `.ma` per analysis.
                 Login-node safe. BrainSEQ needs `coloc_brainseq --stage meta` first.
  --stage besd   per work-list task: nominal statistics -> ESD -> `smr --make-besd`.
  --stage smr    per work-list task: `smr` for every analysis with targets in that cell.
                 `--peqtl-smr` / `--smr-multi` select a sensitivity arm; both reuse the BESD.
  --stage meta   assemble, apply the multiple-testing family, join coloc per probe, classify,
                 report.

Outputs under 05_genetic_anchoring/_m/smr_heidi/gtex/ and .../smr_heidi/brainseq/<arm>/.
"""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data import coloc_modality_contrast as cmc
from isograph_benchmark.real_data.coloc_signal_susie import results_dir, signal_root

SMR_BIN = Path("/ocean/projects/bio260021p/shared/opt/SMR/build/Release/smr")
QTL_SOURCES: tuple[str, ...] = ("gtex", "brainseq")
BRAINSEQ_ARMS: tuple[str, ...] = ("ea_only",)
# The QTL axes each source carries: GTEx expression and intron excision; BrainSEQ's per-gene
# abundance channel and switch coordinate.
SOURCE_MODALITIES: dict[str, tuple[str, ...]] = {"gtex": ("eQTL", "sQTL"),
                                                 "brainseq": ("A_g", "S_g")}
MODALITIES: tuple[str, ...] = SOURCE_MODALITIES["gtex"]
BRAINSEQ_FTYPE: dict[str, str] = {"A_g": "abundance", "S_g": "switch"}

# SMR 1.4.2 defaults (src/SMR.cpp), passed explicitly on every call.
PEQTL_SMR = 5e-8
PEQTL_HEIDI = 1.5654e-3
HEIDI_MIN_M = 3
HEIDI_MAX_M = 20
CIS_WIND_KB = 2000
DIFF_FREQ = 0.2
DIFF_FREQ_PROP = 0.05
# Multi-SNP SMR (`--smr-multi`). SMR prunes the cis SNPs it combines at this LD r2 before
# testing; 0.1 is SMR 1.4.2's own default (src/SMR.cpp:122), passed explicitly like the rest.
LD_MULTI_SNP = 0.1

# Interpretation, pre-specified.
SMR_ALPHA = 0.05
HEIDI_REJECT = 0.01
PP4_CALL = cmc.PP4_CALL
P12_PRIMARY = cmc.P12_PRIMARY
MIN_SHARED_SNPS = cmc.MIN_SHARED_SNPS

# File formats, as SMR 1.4.2 parses them: `.ma` in read_gwas_data (8 columns, header starting
# SNP), the ESD in get_esi_info (9 columns, header starting Chr; beta/se/p are columns 7-9),
# the flist in read_probeinfolst (7 columns, header starting Chr).
MA_COLS = ["SNP", "A1", "A2", "freq", "b", "se", "p", "n"]
ESD_COLS = ["Chr", "SNP", "Bp", "A1", "A2", "Freq", "Beta", "se", "p"]
FLIST_COLS = ["Chr", "ProbeID", "GeneticDistance", "ProbeBp", "Gene", "Orientation", "PathOfEsd"]
SMR_OUT_COLS = ["probeID", "ProbeChr", "Gene", "Probe_bp", "topSNP", "topSNP_chr", "topSNP_bp",
                "A1", "A2", "Freq", "b_GWAS", "se_GWAS", "p_GWAS", "b_eQTL", "se_eQTL",
                "p_eQTL", "b_SMR", "se_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]
# `--smr-multi` writes `.msmr`, not `.smr`, with one extra column inserted before p_HEIDI.
# Verified against the three writers in src/SMR_data_p1.cpp (lines 2039, 2632, 4811), which
# emit an identical header, so the format does not depend on which dispatch path runs.
SMR_MULTI_COL = "p_SMR_multi"
SMR_MULTI_OUT_COLS = [*SMR_OUT_COLS[:SMR_OUT_COLS.index("p_HEIDI")], SMR_MULTI_COL,
                      *SMR_OUT_COLS[SMR_OUT_COLS.index("p_HEIDI"):]]

AGREEMENT: tuple[str, ...] = (
    "coloc_and_smr_heidi_not_rejected",   # coloc call, SMR significant, HEIDI not rejected
    "coloc_and_smr_heidi_rejected",       # coloc call, SMR significant, HEIDI rejects
    "coloc_and_smr_heidi_untestable",     # coloc call, SMR significant, < HEIDI_MIN_M SNPs
    "coloc_smr_not_significant",          # coloc call, SMR instrumented but not significant
    "coloc_no_instrument",                # coloc call, no cis-QTL reaches the SMR threshold
    "smr_without_coloc",                  # SMR significant where coloc scored below the call
    "smr_coloc_not_scored",               # SMR significant at a probe coloc never scored
    "neither_tested_null",                # no coloc call; instrumented, SMR not significant
    "neither_no_instrument",              # no coloc call and never instrumented -- UNTESTED
)

# The SMR-side verdict, carried beside `agreement` because the two axes answer different
# questions. The old vocabulary collapsed "never instrumented" into `neither`, which made the
# largest cell of both reports unreadable: of the 2026-09-11 `neither` rows, 1,761/1,823 (GTEx)
# and 219/231 (BrainSEQ) had no instrument at all and so were never tested for anything.
SMR_STATUS: tuple[str, ...] = (
    "no_instrument",                  # no cis-QTL clears --peqtl-smr; the probe is UNTESTED
    "instrumented_tested_null",       # tested, SMR p above the family threshold
    "smr_signal_heidi_unavailable",   # SMR significant, HEIDI skipped (< HEIDI_MIN_M SNPs)
    "smr_signal_heidi_rejects",       # SMR significant, HEIDI rejects one shared variant
    "smr_heidi_supported",            # SMR significant, HEIDI not rejected
)

# The multi-SNP arm's own verdict. Kept apart from `smr_status` rather than replacing it:
# `p_SMR` is still emitted in `.msmr`, so both arms carry a comparable single-SNP verdict, and
# the interesting quantity is whether a probe significant on its top SNP survives the
# multi-SNP test. `multi_unavailable` is NOT a null -- SMR skips the test when too few cis
# SNPs survive LD pruning, and folding that into "tested and null" would invent evidence.
SMR_MULTI_STATUS: tuple[str, ...] = (
    "no_instrument",                    # no cis-QTL clears --peqtl-smr; UNTESTED either way
    "multi_unavailable",                # instrumented, but the multi-SNP test produced no p
    "multi_tested_null",                # multi-SNP p above its own family threshold
    "multi_signal_heidi_unavailable",   # multi-SNP significant, HEIDI skipped
    "multi_signal_heidi_rejects",       # multi-SNP significant, HEIDI rejects
    "multi_heidi_supported",            # multi-SNP significant, HEIDI not rejected
)

# Instrument strength F = (b/se)^2 for the top cis-QTL SNP. Reported on every row and never
# used to exclude one: a post hoc F filter on an already-thresholded instrument set is its own
# selection. 10 is the conventional weak-instrument marker; at --peqtl-smr 5e-8, z ~ 5.45 and
# F ~ 30, so weak instruments should be rare in the primary arm and common only if a relaxed
# threshold is ever run.
WEAK_F = 10.0

# The two multiple-testing families. A gene's pre-designated probe is confirmatory; its other
# introns localize the event and are corrected separately, so neither the primary threshold is
# inflated by introns nor the introns escape correction when they are discussed.
PROBE_FAMILIES: tuple[str, ...] = ("primary", "secondary")

_AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
_LOCUS_CHR = re.compile(r"_chr(\d+)$")


# --------------------------------------------------------------------------- #
# Layout and the source guard
# --------------------------------------------------------------------------- #
def check_qtl_source(qtl_source: str, arm: str | None = None,
                     failures: list[str] | None = None) -> None:
    """Refuse any QTL source this module cannot use defensibly.

    `failures` is the BrainSEQ QTL-check verdict (`brainseq_qtl_checks.arm_check_failures`);
    left None it is read from disk.
    """
    if qtl_source not in QTL_SOURCES:
        raise SystemExit(f"unknown --qtl-source {qtl_source!r}; choose from {QTL_SOURCES}")
    if qtl_source != "brainseq":
        return
    if arm not in BRAINSEQ_ARMS:
        raise SystemExit(
            f"BrainSEQ SMR must use the `ea_only` arm, not {arm!r}: the cohort is roughly half "
            "African-ancestry while the GWAS and the 1000G EUR LD reference are European, so "
            "an SMR estimate or a HEIDI test on the mixed-ancestry QTL would rest on mismatched LD.")
    if failures is None:
        from isograph_benchmark.real_data.brainseq_qtl_checks import arm_check_failures

        failures = arm_check_failures(arm)
    if failures:
        raise SystemExit(
            f"BrainSEQ `{arm}` has not passed its QTL checks, so no SMR estimate is built on it: "
            + "; ".join(failures) + " (run 02g.brainseq_qtl_checks.sh)")


def modalities(qtl_source: str) -> tuple[str, ...]:
    return SOURCE_MODALITIES[qtl_source]


def source_root(qtl_source: str = "gtex", arm: str | None = None) -> Path:
    """Targets, GWAS `.ma` and BESD for one QTL source (independent of the SMR threshold)."""
    base = stage_out("anchoring.smr") / qtl_source
    if qtl_source == "brainseq":
        if not arm:
            raise SystemExit("the BrainSEQ QTL source needs --arm")
        base = base / arm
    return base


def run_root(qtl_source: str = "gtex", peqtl_smr: float = PEQTL_SMR,
             arm: str | None = None, multi: bool = False) -> Path:
    """SMR runs and meta outputs.

    The primary arm is the bare source root. A non-default instrument threshold and the
    multi-SNP test are each sensitivity arms with their own directory, so neither can
    overwrite the primary result; requesting both composes one directory rather than
    silently collapsing onto either.
    """
    base = source_root(qtl_source, arm)
    tags: list[str] = []
    if not np.isclose(peqtl_smr, PEQTL_SMR):
        tags.append(f"peqtl_smr_{peqtl_smr:g}")
    if multi:
        tags.append("smr_multi")
    return base if not tags else base / "sensitivity" / "__".join(tags)


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(s))


def locus_chr(locus_id: str) -> int:
    m = _LOCUS_CHR.search(str(locus_id))
    if not m:
        raise SystemExit(f"cannot read a chromosome from LOCUS_ID {locus_id!r}")
    return int(m.group(1))


def probe_gene(probe_id: str) -> str:
    """Bare ENSG of a probe: `ENSG...v` (eQTL, A_g, S_g) or `chr:start:end:clu_N_s:ENSG...v` (sQTL)."""
    return str(probe_id).split(":")[-1].split(".")[0]


# --------------------------------------------------------------------------- #
# Stage: prep
# --------------------------------------------------------------------------- #
def build_targets(nom: pd.DataFrame, cells: pd.DataFrame, call: float = PP4_CALL) -> pd.DataFrame:
    """(nominated locus, gene, tissue) cells where the sQTL call holds."""
    keys = ["analysis", "trait", "LOCUS_ID", "gene"]
    c = cells.merge(nom[keys].drop_duplicates(), on=keys, how="inner")
    c = c[c["PP4_sQTL"] >= call].copy()
    c["chr"] = c["LOCUS_ID"].map(locus_chr)
    c = c.rename(columns={"estimator_sQTL": "coloc_estimator_sQTL",
                          "PP4_sQTL": "coloc_PP4_sQTL", "PP4_eQTL": "coloc_PP4_eQTL",
                          "phenotype_id": "coloc_phenotype_id"})
    keep = [*keys, "chr", "symbol", "tissue", "coloc_estimator_sQTL", "coloc_PP4_sQTL",
            "coloc_PP4_eQTL", "coloc_phenotype_id"]
    for k in keep:
        if k not in c.columns:
            c[k] = None
    return c[keep].sort_values(["analysis", "LOCUS_ID", "gene", "tissue"]).reset_index(drop=True)


def build_brainseq_targets(gtex_nom: pd.DataFrame, gtex_cells: pd.DataFrame,
                           bs_nom: pd.DataFrame, bs_cells: pd.DataFrame,
                           regions, gtex_tissue: dict[str, str],
                           call: float = PP4_CALL) -> pd.DataFrame:
    """(nominated locus, gene, BrainSEQ region) cells for the BrainSEQ source.

    The union of the GTEx signal-level nominations and BrainSEQ's own switch-axis nominations,
    crossed with every region. `target_source` says which layer nominated the locus, and
    `coloc_PP4_<axis>` carries BrainSEQ's coloc posterior for that region and axis.
    """
    keys = ["analysis", "trait", "LOCUS_ID", "gene"]
    g = gtex_nom[[*keys, "symbol"]].drop_duplicates(keys).assign(_g=True)
    b = bs_nom[[*keys, "symbol"]].drop_duplicates(keys).assign(_b=True)
    u = g.merge(b, on=keys, how="outer", suffixes=("", "_bs"))
    u["symbol"] = u["symbol"].fillna(u.pop("symbol_bs"))
    in_g = u.pop("_g").notna()
    in_b = u.pop("_b").notna()
    u["target_source"] = np.select([in_g & in_b, in_g], ["both", "gtex"], "brainseq")

    t = u.merge(pd.DataFrame({"tissue": list(regions)}), how="cross")
    t["chr"] = t["LOCUS_ID"].map(locus_chr)
    t["_gt"] = t["tissue"].map(gtex_tissue)
    gc = (gtex_cells[[*keys, "tissue", "PP4_sQTL"]]
          .rename(columns={"tissue": "_gt", "PP4_sQTL": "gtex_PP4_sQTL_matched"})
          .drop_duplicates([*keys, "_gt"]))
    t = t.merge(gc, on=[*keys, "_gt"], how="left").drop(columns="_gt")
    t["gtex_call_in_matched_tissue"] = t["gtex_PP4_sQTL_matched"].fillna(0) >= call

    for ax in SOURCE_MODALITIES["brainseq"]:
        sub = (bs_cells[bs_cells["modality"] == ax][[*keys, "tissue", "PP4", "estimator"]]
               .rename(columns={"PP4": f"coloc_PP4_{ax}", "estimator": f"coloc_estimator_{ax}"})
               .drop_duplicates([*keys, "tissue"]))
        t = t.merge(sub, on=[*keys, "tissue"], how="left")
    t["coloc_phenotype_id"] = None
    return t.sort_values(["analysis", "LOCUS_ID", "gene", "tissue"]).reset_index(drop=True)


def ma_from_loci(loci: list[tuple[int, pd.DataFrame]], excl: pd.DataFrame, n: int) -> pd.DataFrame:
    """Per-locus GWAS (rsid, pos, a1, a2, beta, se, p) as SMR `.ma` rows.

    `a1` is the effect allele in these files (stage 22 re-signs z by it). The long-range-LD
    exclusions are the coloc layer's, applied the same way, so both methods see one GWAS.
    """
    parts = []
    for ch, g in loci:
        g = g.copy()
        for r in excl[excl["chr"] == ch].itertuples(index=False):
            g = g[~((g["pos"] >= r.start) & (g["pos"] <= r.stop))]
        parts.append(g)
    if not parts:
        return pd.DataFrame(columns=MA_COLS)
    g = pd.concat(parts, ignore_index=True)
    g = g[np.isfinite(g["beta"]) & np.isfinite(g["se"]) & (g["se"] > 0)]
    g = g.drop_duplicates("rsid")
    return pd.DataFrame({"SNP": g["rsid"].values, "A1": g["a1"].str.upper().values,
                         "A2": g["a2"].str.upper().values, "freq": np.nan,
                         "b": g["beta"].values, "se": g["se"].values, "p": g["p"].values,
                         "n": int(n)})[MA_COLS]


def write_ma(ma: pd.DataFrame, path: Path) -> None:
    ma.to_csv(path, sep="\t", index=False, na_rep="NA")


def _write_prep(t: pd.DataFrame, dest: Path) -> Path:
    t.to_parquet(dest / "targets.parquet", index=False)
    work = t[["tissue", "chr"]].drop_duplicates().sort_values(["tissue", "chr"])
    work.to_csv(dest / "work_list.tsv", sep="\t", index=False)

    gm = pd.read_csv(cmc.out_dir("switch") / "gwas_meta.tsv", sep="\t").set_index("trait")
    gdir = ensure_dir(dest / "gwas")
    for analysis, sub in t.groupby("analysis"):
        cdir = cmc.coloc_dir(analysis)
        ef = cdir / "exclude_regions.tsv"
        excl = pd.read_csv(ef, sep="\t") if ef.exists() else pd.DataFrame(
            columns=["chr", "start", "stop"])
        loci = [(locus_chr(lid), pd.read_csv(cdir / "susie" / f"{lid}.gwas.tsv", sep="\t"))
                for lid in sorted(set(sub["LOCUS_ID"]))]
        ma = ma_from_loci(loci, excl, n=int(gm.at[sub["trait"].iloc[0], "n_gwas"]))
        write_ma(ma, gdir / f"{analysis}.ma")
        print(f"  {analysis}: {len(ma):,} GWAS SNPs over {len(loci)} loci")
    print(f"  {len(t)} (locus, gene, tissue) targets, {t['gene'].nunique()} genes, "
          f"{len(work)} (tissue, chr) tasks -> {dest}")
    return dest


def run_prep(qtl_source: str = "gtex", arm: str | None = None, call: float = PP4_CALL) -> Path:
    check_qtl_source(qtl_source, arm)
    from isograph_benchmark.real_data.locus_event_audit import load_nominations

    nom, cells = load_nominations("susie", call=call, sqtl_arm="all")
    if qtl_source == "gtex":
        t = build_targets(nom, cells, call=call)
    else:
        from isograph_benchmark.real_data.brainseq_qtl_checks import GTEX_TISSUE
        from isograph_benchmark.real_data.coloc_brainseq import REGIONS, coloc_root

        croot = coloc_root(arm)
        need = [croot / "nominations.parquet", croot / "cells_hierarchy.parquet"]
        missing = [str(p) for p in need if not p.exists()]
        if missing:
            raise SystemExit(f"missing {missing}; run `coloc_brainseq --stage meta` first -- "
                             "BrainSEQ SMR is read against BrainSEQ's own colocalization")
        t = build_brainseq_targets(nom, cells, pd.read_parquet(need[0]),
                                   pd.read_parquet(need[1]), REGIONS, GTEX_TISSUE, call=call)
    if t.empty:
        raise SystemExit("no nominated (locus, gene, tissue) cells; nothing to test")
    return _write_prep(t, ensure_dir(source_root(qtl_source, arm)))


# --------------------------------------------------------------------------- #
# Stage: besd
# --------------------------------------------------------------------------- #
def esd_rows(q: pd.DataFrame, bridge: pd.DataFrame, bim: pd.DataFrame) -> pd.DataFrame:
    """GTEx nominal statistics as SMR ESD rows, keyed by `phenotype_id`.

    `q`: phenotype_id, variant_id, af, pval_nominal, slope, slope_se. `bridge`: variant_id ->
    rsid. `bim`: the reference panel's bchr, rsid, bp, A1, A2 (hg19).
    """
    cols = ["phenotype_id", *ESD_COLS]
    if q.empty:
        return pd.DataFrame(columns=cols)
    d = (q.merge(bridge[["variant_id", "rsid"]], on="variant_id")
          .merge(bim[["rsid", "bchr", "bp", "A1", "A2"]], on="rsid"))
    if d.empty:
        return pd.DataFrame(columns=cols)
    parts = d["variant_id"].str.split("_", expand=True)
    ref, alt = parts[2].str.upper(), parts[3].str.upper()
    p1, p2 = d["A1"].astype(str).str.upper(), d["A2"].astype(str).str.upper()
    same = ((ref == p1) & (alt == p2)) | ((ref == p2) & (alt == p1))
    ambiguous = pd.Series([(a, b) in _AMBIGUOUS for a, b in zip(ref, alt)], index=d.index)
    # SMR refuses a frequency of exactly 0 or 1, and a zero SE has no ratio estimate.
    ok = same & ~ambiguous & (d["af"] > 0) & (d["af"] < 1) & (d["slope_se"] > 0)
    d, ref, alt = d[ok], ref[ok], alt[ok]
    out = pd.DataFrame({"phenotype_id": d["phenotype_id"].values,
                        "Chr": d["bchr"].astype(int).values, "SNP": d["rsid"].values,
                        "Bp": d["bp"].astype(int).values, "A1": alt.values, "A2": ref.values,
                        "Freq": d["af"].values, "Beta": d["slope"].values,
                        "se": d["slope_se"].values, "p": d["pval_nominal"].values})
    return out.drop_duplicates(["phenotype_id", "SNP"]).reset_index(drop=True)


def esd_rows_brainseq(q: pd.DataFrame, alleles: pd.DataFrame, bim: pd.DataFrame,
                      pin: pd.DataFrame | None = None) -> pd.DataFrame:
    """BrainSEQ nominal statistics as SMR ESD rows, keyed by `phenotype_id`.

    `q`: phenotype_id, variant_id (an rsID), af, pval_nominal, slope, slope_se. `alleles`: the
    arm's .pvar REF/ALT per rsID (`bs_ref`, `bs_alt`). `bim`: the reference panel (hg19).
    `pin`: for S_g, the sign pin per gene (`gene`, `sign`, `pinnable`). A pinnable gene's slope
    is multiplied by its pin so b_SMR points along the discovery switch axis; a gene the pin
    cannot orient keeps its mapped sign, and the meta stage flags it.
    """
    cols = ["phenotype_id", *ESD_COLS]
    if q.empty:
        return pd.DataFrame(columns=cols)
    d = (q.merge(alleles[["variant_id", "bs_ref", "bs_alt"]], on="variant_id")
          .merge(bim[["rsid", "bchr", "bp", "A1", "A2"]], left_on="variant_id", right_on="rsid"))
    if d.empty:
        return pd.DataFrame(columns=cols)
    ref, alt = d["bs_ref"].astype(str).str.upper(), d["bs_alt"].astype(str).str.upper()
    p1, p2 = d["A1"].astype(str).str.upper(), d["A2"].astype(str).str.upper()
    same = ((ref == p1) & (alt == p2)) | ((ref == p2) & (alt == p1))
    ambiguous = pd.Series([(a, b) in _AMBIGUOUS for a, b in zip(ref, alt)], index=d.index)
    ok = same & ~ambiguous & (d["af"] > 0) & (d["af"] < 1) & (d["slope_se"] > 0)

    beta = d["slope"].astype(float).to_numpy()
    if pin is not None and len(pin):
        p = pin.drop_duplicates("gene").set_index("gene")
        sign = p["sign"].astype(float).where(p["pinnable"].astype(bool), 1.0)
        beta = beta * d["phenotype_id"].map(probe_gene).map(sign).fillna(1.0).to_numpy()

    keep = ok.to_numpy()
    d = d[keep]
    out = pd.DataFrame({"phenotype_id": d["phenotype_id"].values,
                        "Chr": d["bchr"].astype(int).values, "SNP": d["rsid"].values,
                        "Bp": d["bp"].astype(int).values, "A1": alt[keep].values,
                        "A2": ref[keep].values, "Freq": d["af"].values, "Beta": beta[keep],
                        "se": d["slope_se"].values, "p": d["pval_nominal"].values})
    return out.drop_duplicates(["phenotype_id", "SNP"]).reset_index(drop=True)


def gene_tss_hg19(genes: pd.DataFrame) -> pd.DataFrame:
    """(gene, symbol) -> hg19 TSS and strand from MAGMA's NCBI37.3 gene.loc."""
    from isograph_benchmark.real_data.coloc_prep import GENE_LOC_HG19

    loc = pd.read_csv(GENE_LOC_HG19, sep=r"\s+", header=None,
                      names=["entrez", "chr", "start", "stop", "strand", "symbol"],
                      dtype={"chr": str})
    loc = loc.drop_duplicates("symbol")
    loc["tss"] = np.where(loc["strand"] == "-", loc["stop"], loc["start"])
    return (genes[["gene", "symbol"]].drop_duplicates("gene")
            .merge(loc[["symbol", "strand", "tss"]], on="symbol", how="left"))


def flist_rows(esd: pd.DataFrame, genes: pd.DataFrame, esd_dir: Path) -> pd.DataFrame:
    """One `--eqtl-flist` row per probe.

    A probe whose gene has no hg19 TSS falls back to the median ESD position with orientation
    NA. SMR uses the probe position only to centre its 2 Mb cis window, and the QTL cis variants
    (within 1 Mb of the TSS in both sources) already sit inside it.
    """
    g = genes.drop_duplicates("gene").set_index("gene")
    rows = []
    for pid, blk in esd.groupby("phenotype_id", sort=True):
        gene = probe_gene(pid)
        has = gene in g.index
        sym = str(g.at[gene, "symbol"]) if has and pd.notna(g.at[gene, "symbol"]) else gene
        if has and "tss" in g.columns and pd.notna(g.at[gene, "tss"]):
            bp, orient = int(g.at[gene, "tss"]), str(g.at[gene, "strand"])
        else:
            bp, orient = int(blk["Bp"].median()), "NA"
        rows.append({"Chr": int(blk["Chr"].iloc[0]), "ProbeID": pid, "GeneticDistance": 0,
                     "ProbeBp": bp, "Gene": sym, "Orientation": orient,
                     "PathOfEsd": str(esd_dir / f"{_safe(pid)}.esd")})
    return pd.DataFrame(rows, columns=FLIST_COLS)


def read_gtex_qtl(tissue: str, chrom: int, genes: set[str], modality: str) -> pd.DataFrame:
    """Nominal cis statistics for the target genes, by parquet predicate pushdown."""
    import pyarrow.compute as pc
    import pyarrow.dataset as pds

    if modality == "eQTL":
        f = (cmc.GTEX_V11 / "GTEx_Analysis_v11_eQTL_all_associations"
             / f"{tissue}.v11.allpairs.chr{chrom}.parquet")
        col = "gene_id"
    else:
        f = (cmc.GTEX_V11 / "GTEx_Analysis_v11_sQTL_all_associations"
             / f"{tissue}.v11.cis_sqtl.allpairs.chr{chrom}.parquet")
        col = "phenotype_id"
    cols = ["phenotype_id", "variant_id", "af", "pval_nominal", "slope", "slope_se"]
    if not f.exists() or not genes:
        return pd.DataFrame(columns=cols)
    gs = sorted(genes)
    expr = pc.match_substring(pds.field(col), gs[0])
    for g in gs[1:]:
        expr = expr | pc.match_substring(pds.field(col), g)
    t = (pds.dataset(str(f))
         .to_table(columns=[col, "variant_id", "af", "pval_nominal", "slope", "slope_se"],
                   filter=expr)
         .to_pandas().rename(columns={col: "phenotype_id"}))
    # A substring hit is not an identity: keep only probes whose gene IS a target.
    return t[t["phenotype_id"].map(probe_gene).isin(genes)][cols]


def read_brainseq_qtl(arm: str, region: str, chrom: int, genes: set[str],
                      modality: str) -> pd.DataFrame:
    """BrainSEQ nominal cis statistics for the target genes on one axis and chromosome."""
    import pyarrow.dataset as pds

    from isograph_benchmark.real_data.brainseq_switch_qtl import out_dir as qtl_dir

    cols = ["phenotype_id", "variant_id", "af", "pval_nominal", "slope", "slope_se"]
    qdir = qtl_dir(arm, region) / "qtl"
    ftype = BRAINSEQ_FTYPE[modality]
    perm = qdir / f"cis_qtl_{ftype}.parquet"
    files = sorted(qdir.glob(f"{ftype}.chr{chrom}.cis_qtl_pairs.*.parquet"))
    if not genes or not perm.exists() or not files:
        return pd.DataFrame(columns=cols)
    ids = pd.read_parquet(perm, columns=["phenotype_id"])["phenotype_id"].astype(str)
    want = sorted(ids[ids.map(probe_gene).isin(genes)])
    if not want:
        return pd.DataFrame(columns=cols)
    return (pds.dataset([str(f) for f in files], format="parquet")
            .to_table(columns=cols, filter=pds.field("phenotype_id").isin(want)).to_pandas())


def load_sign_pin(arm: str, region: str) -> pd.DataFrame:
    f = stage_out("anchoring.brainseq_qtl", arm) / "checks" / f"sign_pin_{region}.parquet"
    if not f.exists():
        raise SystemExit(f"missing {f}; run 02g.brainseq_qtl_checks.sh for {arm}")
    return pd.read_parquet(f)[["gene", "r", "sign", "flipped", "axis_unstable", "pinnable"]]


def _read_bim(chrom: int) -> pd.DataFrame:
    from isograph_benchmark.real_data.coloc_prep import PANEL_DIR

    return pd.read_csv(PANEL_DIR / f"1000G.EUR.QC.{chrom}.bim", sep=r"\s+", header=None,
                       names=["bchr", "rsid", "cm", "bp", "A1", "A2"],
                       dtype={"A1": str, "A2": str})


def _run(cmd: list[str], log: Path) -> int:
    with open(log, "w") as fh:
        fh.write(" ".join(cmd) + "\n\n")
        fh.flush()
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode


def _tasks(root: Path, task: int | None) -> pd.DataFrame:
    wf = root / "work_list.tsv"
    if not wf.exists():
        raise SystemExit(f"missing {wf}; run --stage prep first")
    work = pd.read_csv(wf, sep="\t")
    if task is None:
        return work
    if not 0 <= task < len(work):
        raise SystemExit(f"--task {task} out of range (0-{len(work) - 1})")
    return work.iloc[[task]]


def run_besd(qtl_source: str = "gtex", task: int | None = None, arm: str | None = None) -> None:
    check_qtl_source(qtl_source, arm)
    root = source_root(qtl_source, arm)
    targets = pd.read_parquet(root / "targets.parquet")
    for w in _tasks(root, task).itertuples(index=False):
        tissue, ch = str(w.tissue), int(w.chr)
        sub = targets[(targets["tissue"] == tissue) & (targets["chr"] == ch)]
        genes = set(sub["gene"])
        bim = _read_bim(ch)
        tss = gene_tss_hg19(sub[["gene", "symbol"]])
        if qtl_source == "gtex":
            bridge = pd.read_parquet(cmc.BRIDGE_DIR / f"chr{ch}.parquet")
            esds = {mod: esd_rows(read_gtex_qtl(tissue, ch, genes, mod), bridge, bim)
                    for mod in modalities(qtl_source)}
        else:
            from isograph_benchmark.real_data.brainseq_qtl_checks import read_pvar_alleles

            qs = {mod: read_brainseq_qtl(arm, tissue, ch, genes, mod)
                  for mod in modalities(qtl_source)}
            rsids = set().union(*(set(q["variant_id"]) for q in qs.values()))
            alleles = (read_pvar_alleles(arm, ch, rsids) if rsids else
                       pd.DataFrame(columns=["variant_id", "bs_ref", "bs_alt"]))
            esds = {mod: esd_rows_brainseq(q, alleles, bim,
                                           load_sign_pin(arm, tissue) if mod == "S_g" else None)
                    for mod, q in qs.items()}
        for mod, esd in esds.items():
            if esd.empty:
                print(f"  {tissue} chr{ch} {mod}: no usable QTL rows")
                continue
            prefix = ensure_dir(root / "besd" / tissue) / f"{mod}.chr{ch}"
            esd_dir = ensure_dir(Path(f"{prefix}.esd"))
            for pid, blk in esd.groupby("phenotype_id"):
                blk[ESD_COLS].to_csv(esd_dir / f"{_safe(pid)}.esd", sep="\t", index=False)
            flist = Path(f"{prefix}.flist")
            flist_rows(esd, tss, esd_dir).to_csv(flist, sep="\t", index=False)
            rc = _run([str(SMR_BIN), "--eqtl-flist", str(flist), "--make-besd",
                       "--out", str(prefix)], Path(f"{prefix}.make_besd.log"))
            print(f"  {tissue} chr{ch} {mod}: {esd['phenotype_id'].nunique()} probes, "
                  f"{len(esd):,} rows -> {prefix}.besd (rc={rc})")
            if rc != 0:
                raise SystemExit(f"smr --make-besd failed (rc={rc}); see {prefix}.make_besd.log")


# --------------------------------------------------------------------------- #
# Stage: smr
# --------------------------------------------------------------------------- #
def smr_out_suffix(multi: bool = False) -> str:
    """SMR names its table `.msmr` under `--smr-multi` and `.smr` otherwise."""
    return ".msmr" if multi else ".smr"


def smr_command(bfile: Path, ma: Path, besd: Path, probes: Path, out: Path,
                peqtl_smr: float = PEQTL_SMR, threads: int = 4,
                multi: bool = False) -> list[str]:
    cmd = [str(SMR_BIN), "--bfile", str(bfile), "--gwas-summary", str(ma),
           "--beqtl-summary", str(besd), "--extract-probe", str(probes), "--out", str(out),
           "--peqtl-smr", f"{peqtl_smr:g}", "--peqtl-heidi", f"{PEQTL_HEIDI:g}",
           "--heidi-min-m", str(HEIDI_MIN_M), "--heidi-max-m", str(HEIDI_MAX_M),
           "--cis-wind", str(CIS_WIND_KB), "--diff-freq", f"{DIFF_FREQ:g}",
           "--diff-freq-prop", f"{DIFF_FREQ_PROP:g}", "--thread-num", str(threads)]
    if multi:
        cmd += ["--smr-multi", "--ld-multi-snp", f"{LD_MULTI_SNP:g}"]
    return cmd


def run_smr(qtl_source: str = "gtex", task: int | None = None, peqtl_smr: float = PEQTL_SMR,
            threads: int = 4, arm: str | None = None, multi: bool = False) -> None:
    check_qtl_source(qtl_source, arm)
    from isograph_benchmark.real_data.coloc_prep import PANEL_DIR

    root, dest = source_root(qtl_source, arm), run_root(qtl_source, peqtl_smr, arm, multi)
    targets = pd.read_parquet(root / "targets.parquet")
    for w in _tasks(root, task).itertuples(index=False):
        tissue, ch = str(w.tissue), int(w.chr)
        cell = targets[(targets["tissue"] == tissue) & (targets["chr"] == ch)]
        for analysis, sub in cell.groupby("analysis"):
            ma = root / "gwas" / f"{analysis}.ma"
            for mod in modalities(qtl_source):
                besd = root / "besd" / tissue / f"{mod}.chr{ch}"
                if not Path(f"{besd}.besd").exists():
                    continue
                epi = pd.read_csv(f"{besd}.epi", sep="\t", header=None,
                                  names=["chr", "probeID", "cm", "bp", "gene", "orientation"])
                probes = epi.loc[epi["probeID"].map(probe_gene).isin(set(sub["gene"])),
                                 "probeID"]
                if probes.empty:
                    continue
                odir = ensure_dir(dest / "smr" / analysis / tissue)
                pf = odir / f"{mod}.chr{ch}.probes"
                probes.to_csv(pf, index=False, header=False)
                out = odir / f"{mod}.chr{ch}"
                cmd = smr_command(PANEL_DIR / f"1000G.EUR.QC.{ch}", ma, besd, pf, out,
                                  peqtl_smr=peqtl_smr, threads=threads, multi=multi)
                rc = _run(cmd, Path(f"{out}.log"))
                # SMR exits non-zero when no probe clears the instrument threshold; that is
                # an outcome (`coloc_no_instrument`), so it is recorded, not raised.
                print(f"  {analysis} {tissue} chr{ch} {mod}: {len(probes)} probes (rc={rc})")


# --------------------------------------------------------------------------- #
# Stage: meta
# --------------------------------------------------------------------------- #
def read_smr(path: Path) -> pd.DataFrame:
    """One SMR table. `.msmr` (from `--smr-multi`) carries one extra column; both are checked
    against the exact header SMR 1.4.2 writes, so a format change fails loudly here."""
    d = pd.read_csv(path, sep="\t")
    want = SMR_MULTI_OUT_COLS if Path(path).suffix == ".msmr" else SMR_OUT_COLS
    missing = [c for c in want if c not in d.columns]
    if missing:
        raise SystemExit(f"{path} lacks SMR output columns {missing}")
    return d


def smr_threshold(n_instrumented: int, alpha: float = SMR_ALPHA) -> float:
    """Bonferroni threshold for ONE family; the caller decides which family a probe is in."""
    return alpha / n_instrumented if n_instrumented > 0 else np.nan


def instrumented_family_sizes(d: pd.DataFrame) -> pd.DataFrame:
    """Instrumented probe count per (analysis, modality, probe_family).

    A distinct function because this count IS the multiple-testing denominator, and it was
    the thing that went wrong: pooling the primary and secondary families gave AD sQTL a
    denominator of 77 while the report printed the primary count, 31. A probe is counted once
    per (tissue, probeID) -- the same probe in two tissues is two tests -- and only if it was
    instrumented, because an untested probe carries no p-value to correct.
    """
    ins = d[d["p_SMR"].notna()]
    if ins.empty:
        return pd.DataFrame(columns=["analysis", "modality", "probe_family",
                                     "n_instrumented_family"])
    return (ins.drop_duplicates(["analysis", "modality", "probe_family", "tissue", "probeID"])
               .groupby(["analysis", "modality", "probe_family"]).size()
               .rename("n_instrumented_family").reset_index())


def multi_family_sizes(d: pd.DataFrame) -> pd.DataFrame:
    """Probes the multi-SNP test actually produced a p-value for, per family.

    Its own denominator, not the single-SNP one: SMR skips the multi-SNP test where too few
    cis SNPs survive LD pruning, so correcting it over every instrumented probe would divide
    by tests that were never run and make the arm conservative for the wrong reason.
    """
    ins = d[d[SMR_MULTI_COL].notna()] if SMR_MULTI_COL in d.columns else d.iloc[:0]
    if ins.empty:
        return pd.DataFrame(columns=["analysis", "modality", "probe_family", "n_multi_family"])
    return (ins.drop_duplicates(["analysis", "modality", "probe_family", "tissue", "probeID"])
               .groupby(["analysis", "modality", "probe_family"]).size()
               .rename("n_multi_family").reset_index())


def smr_status(p_smr, p_heidi, nsnp_heidi, threshold) -> str:
    """The SMR-side verdict alone, with no coloc input (see `SMR_STATUS`).

    HEIDI is reported as *unavailable* rather than folded into a null when SMR ran but HEIDI
    could not: SMR skips the test below `HEIDI_MIN_M` SNPs, and a handful of SNPs is a weak
    basis for either verdict.
    """
    if pd.isna(p_smr):
        return "no_instrument"
    if not (pd.notna(threshold) and p_smr <= threshold):
        return "instrumented_tested_null"
    if pd.isna(p_heidi) or pd.isna(nsnp_heidi) or nsnp_heidi < HEIDI_MIN_M:
        return "smr_signal_heidi_unavailable"
    return "smr_signal_heidi_rejects" if p_heidi < HEIDI_REJECT else "smr_heidi_supported"


def smr_multi_status(p_smr, p_smr_multi, p_heidi, nsnp_heidi, threshold) -> str:
    """The multi-SNP arm's verdict (see `SMR_MULTI_STATUS`).

    `p_smr` decides only whether the probe was instrumented at all; the verdict itself is the
    multi-SNP p-value against the multi-SNP family's own threshold. A probe SMR instrumented
    but could not run the multi-SNP test on is `multi_unavailable`, never a null.
    """
    if pd.isna(p_smr):
        return "no_instrument"
    if pd.isna(p_smr_multi):
        return "multi_unavailable"
    if not (pd.notna(threshold) and p_smr_multi <= threshold):
        return "multi_tested_null"
    if pd.isna(p_heidi) or pd.isna(nsnp_heidi) or nsnp_heidi < HEIDI_MIN_M:
        return "multi_signal_heidi_unavailable"
    return "multi_signal_heidi_rejects" if p_heidi < HEIDI_REJECT else "multi_heidi_supported"


def classify(pp4, p_smr, p_heidi, nsnp_heidi, threshold, call: float = PP4_CALL) -> str:
    """The coloc x SMR cross-classification, built on `smr_status` so the two agree."""
    coloc = pd.notna(pp4) and pp4 >= call
    st = smr_status(p_smr, p_heidi, nsnp_heidi, threshold)
    if st == "no_instrument":
        return "coloc_no_instrument" if coloc else "neither_no_instrument"
    if st == "instrumented_tested_null":
        return "coloc_smr_not_significant" if coloc else "neither_tested_null"
    if not coloc:
        # A probe coloc never scored is not a coloc negative; keep it apart from one that was
        # scored and fell short.
        return "smr_without_coloc" if pd.notna(pp4) else "smr_coloc_not_scored"
    return {"smr_signal_heidi_unavailable": "coloc_and_smr_heidi_untestable",
            "smr_signal_heidi_rejects": "coloc_and_smr_heidi_rejected",
            "smr_heidi_supported": "coloc_and_smr_heidi_not_rejected"}[st]


# --------------------------------------------------------------------------- #
# SNP attrition into HEIDI
# --------------------------------------------------------------------------- #
# Every SMR run prints what it kept at each harmonization step, and nothing read those logs
# until 2026-09-11, when a BrainSEQ run turned out to have retained 2,677 of 16,212 BESD SNPs
# (16.5%) against a ~92% norm elsewhere. HEIDI is computed on whatever survives, so the
# retention IS a validity statistic, not a curiosity: a thin or unrepresentative surviving set
# can look like "HEIDI underpowered" when the real cause is allele, id or panel-coverage loss.
_ATTRITION = {
    "n_besd_snps": re.compile(r"(\d+) SNPs to be included from \[[^\]]*\.esi\]"),
    "n_panel_snps": re.compile(r"(\d+) SNPs to be included from \[[^\]]*\.bim\]"),
    # SMR prints this as "included after allele checking", but it is the count AFTER the
    # three-way BESD n reference n GWAS intersection. Measured 2026-09-11: the GWAS side does
    # essentially all the cutting, and allele checking proper removes 0-4 SNPs per run. Naming
    # it after alleles invites reading a GWAS-coverage number as an allele-QC failure.
    "n_smr_shared": re.compile(r"(\d+) SNPs are included after allele checking"),
    "n_gwas_snps": re.compile(r"GWAS summary data of (\d+) SNPs"),
    "n_probes_epi": re.compile(r"(\d+) Probes to be included from \[[^\]]*\.epi\]"),
    "n_probes_besd": re.compile(r"summary data of (\d+) Probes to be included from "
                                r"\[[^\]]*\.besd\]"),
    # SMR drops SNPs whose reference and QTL allele frequencies disagree by more than
    # --diff-freq in more than --diff-freq-prop of probes; it says so only when it happens.
    "n_freq_mismatch": re.compile(r"(\d+) SNPs? (?:are|were|is) (?:excluded|removed)[^.\n]*"
                                  r"(?:frequenc|freq)", re.IGNORECASE),
}


def parse_smr_log(path: Path) -> dict:
    """The per-run harmonization counts SMR prints. Missing lines stay NaN, never 0."""
    text = Path(path).read_text(errors="replace")
    out: dict[str, float] = {}
    for key, rx in _ATTRITION.items():
        m = rx.search(text)
        out[key] = int(m.group(1)) if m else np.nan
    return out


def _attrition_steps(source_dir: Path, tissue: str, mod: str, ch: int, analysis: str,
                     bims: dict, mas: dict) -> dict:
    """The harmonization steps SMR does not separate, computed from its own inputs."""
    esi = source_dir / "besd" / tissue / f"{mod}.chr{ch}.esi"
    ma = source_dir / "gwas" / f"{analysis}.ma"
    if not (esi.exists() and ma.exists()):
        return {}
    e = pd.read_csv(esi, sep="\t", header=None,
                    names=["chr", "rsid", "cm", "bp", "A1", "A2", "freq"],
                    dtype={"rsid": str, "A1": str, "A2": str})
    if ch not in bims:
        bims[ch] = _read_bim(ch)[["rsid", "A1", "A2"]]
    if analysis not in mas:
        mas[analysis] = set(pd.read_csv(ma, sep=r"\s+", usecols=["SNP"])["SNP"].astype(str))
    m = e.merge(bims[ch], on="rsid", suffixes=("_e", "_b"))
    pair_e = list(zip(m["A1_e"].str.upper(), m["A2_e"].str.upper()))
    pair_b = list(zip(m["A1_b"].str.upper(), m["A2_b"].str.upper()))
    ok = [set(a) == set(b) for a, b in zip(pair_e, pair_b)]
    matched = m[ok]
    return {
        "n_besd_snps_counted": int(len(e)),
        "n_panel_shared": int(len(m)),
        "n_allele_match": int(len(matched)),
        "n_strand_ambiguous": int(sum(p in _AMBIGUOUS for p in
                                      (pair_e[i] for i, k in enumerate(ok) if k))),
        "n_gwas_shared": int(matched["rsid"].isin(mas[analysis]).sum()),
    }


def collect_attrition(run_dir: Path, source_dir: Path | None = None) -> pd.DataFrame:
    """Per (analysis, tissue, modality, chr) attrition trace into HEIDI.

    Per-RUN resolution, which is what the inputs carry; the per-PROBE end of the trace is
    `nsnp_HEIDI` in `smr_results.parquet` (itself capped at `HEIDI_MAX_M`).

    The steps are computed from the BESD, the LD panel and the GWAS `.ma` directly rather than
    read off SMR's log, because the log collapses them into one line whose name
    ("after allele checking") describes the smallest of the three filters. Measured over all 200
    BrainSEQ and 60 GTEx runs on 2026-09-11: rsID overlap with the panel and allele agreement are
    both 100%, and the residual unexplained by the GWAS intersection is 0-4 SNPs per run. What
    looks like "attrition" is overwhelmingly GWAS locus coverage, which is by construction --
    the `.ma` files hold only the SNPs of the selected GWAS loci.
    """
    rows: list[dict] = []
    bims: dict = {}
    mas: dict = {}
    for log in sorted((run_dir / "smr").rglob("*.log")):
        m = _PROBE_FILE.match(log.name.replace(".log", ".probes"))
        if not m:
            continue
        analysis, tissue = log.parent.parent.name, log.parent.name
        mod, ch = m.group(1), int(m.group(2))
        rec = {"analysis": analysis, "tissue": tissue, "modality": mod, "chr": ch}
        rec.update(parse_smr_log(log))
        if source_dir is not None:
            rec.update(_attrition_steps(source_dir, tissue, mod, ch, analysis, bims, mas))
        rows.append(rec)
    if not rows:
        return pd.DataFrame()
    d = pd.DataFrame(rows)
    # The QC fractions are the panel and allele steps. BESD-SNP "retention" is NOT one: its
    # denominator is imputation density, so BrainSEQ (TOPMed, dense) scores far below GTEx
    # (pre-filtered all-pairs) with nothing wrong in either -- 43.8% against 94.7% on the
    # 2026-09-11 runs. Reporting that ratio as quality would manufacture a defect.
    if "n_panel_shared" in d.columns:
        d["frac_panel_shared"] = d["n_panel_shared"] / d["n_besd_snps_counted"]
        d["frac_allele_match"] = d["n_allele_match"] / d["n_panel_shared"]
        # The only genuine harmonization loss: SNPs in BESD n panel n GWAS that SMR still
        # dropped (frequency-discrepancy filtering, duplicates).
        d["n_harmonization_residual"] = d["n_gwas_shared"] - d["n_smr_shared"]
    return d.sort_values(["analysis", "tissue", "modality", "chr"]).reset_index(drop=True)


def tss_fallback_genes(targets: pd.DataFrame) -> set[str]:
    """Genes that took `flist_rows`' median-ESD probe position for want of an hg19 TSS.

    Recomputed here from the same inputs rather than recorded at BESD time, so the flag works
    on runs already on disk. A fallback moves the 2 Mb cis window off the promoter, which
    changes the eligible instrument set -- QC-failed until someone looks at the gene.
    """
    if targets.empty or "symbol" not in targets.columns:
        return set()
    t = gene_tss_hg19(targets[["gene", "symbol"]].drop_duplicates("gene"))
    return set(t.loc[t["tss"].isna(), "gene"].astype(str))


def probe_coloc(pairs: pd.DataFrame, hierarchy: pd.DataFrame, srep: pd.DataFrame) -> pd.DataFrame:
    """The GTEx coloc posterior for each SMR probe, at the resolution the probe has.

    SMR tests one phenotype per probe, so the comparison is per intron, not per cell: a
    signal-level PP4 for that intron (GTEx-matched pairs, primary prior, shared-SNP floor)
    where SuSiE scored it, else the abf posterior -- which exists only for GTEx's
    representative intron, so no other intron of an abf cell inherits it.
    """
    k = ["analysis", "trait", "LOCUS_ID", "tissue", "modality"]
    sp = pairs[(pairs["p12"] == P12_PRIMARY)
               & pairs["cs_matches_gtex"].fillna(False).astype(bool)
               & (pairs["n_shared"] >= MIN_SHARED_SNPS)].copy()
    sp["probe_key"] = np.where(sp["modality"] == "eQTL", sp["gene"].astype(str),
                               sp["phenotype_id"].astype(str))
    s = (sp.groupby([*k, "probe_key"], as_index=False)["PP4"].max()
           .assign(coloc_estimator_probe="susie"))
    ha = hierarchy[hierarchy["estimator"] == "abf"]
    e = ha[ha["modality"] == "eQTL"].assign(probe_key=lambda x: x["gene"].astype(str))
    q = ha[ha["modality"] == "sQTL"].merge(
        srep[["tissue", "gene", "phenotype_id"]].rename(columns={"phenotype_id": "probe_key"}),
        on=["tissue", "gene"], how="inner")
    a = pd.concat([e, q], ignore_index=True)
    a = a[[*k, "probe_key", "PP4"]].assign(coloc_estimator_probe="abf")
    out = pd.concat([s, a], ignore_index=True).rename(columns={"PP4": "coloc_PP4_probe"})
    return out.drop_duplicates([*k, "probe_key"], keep="first").reset_index(drop=True)


def probe_coloc_brainseq(hierarchy: pd.DataFrame) -> pd.DataFrame:
    """BrainSEQ's coloc posterior per (locus, region, axis, gene): one phenotype per gene, so
    the hierarchy cell IS the probe."""
    k = ["analysis", "trait", "LOCUS_ID", "tissue", "modality"]
    return (hierarchy[[*k, "gene", "PP4", "estimator"]]
            .rename(columns={"gene": "probe_key", "PP4": "coloc_PP4_probe",
                             "estimator": "coloc_estimator_probe"})
            .drop_duplicates([*k, "probe_key"]).reset_index(drop=True))


def sign_pin_flags(arm: str, regions) -> pd.DataFrame:
    """Per (gene, region): whether an S_g probe's b_SMR is oriented to the discovery axis."""
    parts = []
    for r in regions:
        p = load_sign_pin(arm, r)
        parts.append(pd.DataFrame({
            "gene": p["gene"], "tissue": r, "modality": "S_g", "pin_r": p["r"],
            "sign_pinned": p["pinnable"].astype(bool),
            "axis_unstable": p["axis_unstable"].astype(bool)}))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        columns=["gene", "tissue", "modality", "pin_r", "sign_pinned", "axis_unstable"])


_PROBE_FILE = re.compile(r"^(eQTL|sQTL|A_g|S_g)\.chr(\d+)\.probes$")


def run_meta(qtl_source: str = "gtex", peqtl_smr: float = PEQTL_SMR,
             arm: str | None = None, multi: bool = False) -> Path:
    check_qtl_source(qtl_source, arm)
    root = source_root(qtl_source, arm)
    dest = run_root(qtl_source, peqtl_smr, arm, multi)
    targets = pd.read_parquet(root / "targets.parquet")
    suffix = smr_out_suffix(multi)

    tested, results, runs = [], [], []
    for pf in sorted((dest / "smr").rglob("*.probes")):
        m = _PROBE_FILE.match(pf.name)
        if not m:
            continue
        analysis, tissue, mod = pf.parent.parent.name, pf.parent.name, m.group(1)
        pr = pd.read_csv(pf, header=None, names=["probeID"]).assign(
            analysis=analysis, tissue=tissue, modality=mod)
        tested.append(pr)
        sf = Path(str(pf)[:-len(".probes")] + suffix)
        runs.append({"analysis": analysis, "tissue": tissue, "modality": mod,
                     "chr": int(m.group(2)), "n_probes": len(pr), "smr_file": sf.exists()})
        if sf.exists():
            s = read_smr(sf)
            if len(s):  # a header-only table is a run where no probe was instrumented
                results.append(s.assign(analysis=analysis, tissue=tissue, modality=mod))
    if not tested:
        raise SystemExit(f"no SMR runs under {dest / 'smr'}; run --stage smr first"
                         + (" --smr-multi" if multi else ""))

    d = pd.concat(tested, ignore_index=True)
    stat = ["topSNP", "A1", "A2", "b_GWAS", "p_GWAS", "b_eQTL", "se_eQTL", "p_eQTL",
            "b_SMR", "se_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]
    if multi:
        stat.append(SMR_MULTI_COL)
    if results:
        r = pd.concat(results, ignore_index=True)[["analysis", "tissue", "modality", "probeID",
                                                   *stat]]
        d = d.merge(r, on=["analysis", "tissue", "modality", "probeID"], how="left")
    else:
        for c in stat:
            d[c] = np.nan
    d["gene"] = d["probeID"].map(probe_gene)
    tcols = ["analysis", "trait", "LOCUS_ID", "gene", "symbol", "tissue", "coloc_phenotype_id",
             *(c for c in ("target_source", "gtex_call_in_matched_tissue") if c in targets.columns)]
    d = d.merge(targets[tcols], on=["analysis", "gene", "tissue"], how="inner")

    keys = ["analysis", "trait", "LOCUS_ID", "tissue", "modality", "probe_key"]
    if qtl_source == "gtex":
        d["probe_key"] = np.where(d["modality"] == "eQTL", d["gene"], d["probeID"])
        sig = results_dir(signal_root("switch"), "all")
        pairs = pd.read_parquet(sig / "signal_pairs.parquet")
        hier = pd.read_parquet(sig / "cells_hierarchy.parquet")
        srep = pd.read_parquet(cmc.out_dir("switch") / "sqtl_representative.parquet")
        d = d.merge(probe_coloc(pairs, hier, srep), on=keys, how="left")
        d["primary_probe"] = (d["modality"] == "eQTL") | (
            (d["probeID"] == d["coloc_phenotype_id"])
            | (d["coloc_phenotype_id"].isna() & (d["coloc_estimator_probe"] == "abf")))
    else:
        from isograph_benchmark.real_data.coloc_brainseq import coloc_root

        d["probe_key"] = d["gene"]
        hier = pd.read_parquet(coloc_root(arm) / "cells_hierarchy.parquet")
        d = d.merge(probe_coloc_brainseq(hier), on=keys, how="left")
        d["primary_probe"] = True
        d = d.merge(sign_pin_flags(arm, sorted(d["tissue"].unique())),
                    on=["gene", "tissue", "modality"], how="left")

    # TWO multiple-testing families, corrected apart (2026-09-11 review). Pooling them and then
    # printing only the primary count made the reported threshold irreproducible: AD sQTL showed
    # "31 instrumented" beside a threshold of 0.05/77, because the denominator silently included
    # the 46 instrumented non-primary introns. Now each probe is corrected inside its own family.
    d["probe_family"] = np.where(d["primary_probe"], "primary", "secondary")
    fam = instrumented_family_sizes(d)
    d = d.merge(fam, on=["analysis", "modality", "probe_family"], how="left")
    d["n_instrumented_family"] = d["n_instrumented_family"].fillna(0).astype(int)
    d["smr_threshold"] = d["n_instrumented_family"].map(smr_threshold)
    if multi:
        d = d.merge(multi_family_sizes(d), on=["analysis", "modality", "probe_family"],
                    how="left")
        d["n_multi_family"] = d["n_multi_family"].fillna(0).astype(int)
        d["smr_multi_threshold"] = d["n_multi_family"].map(smr_threshold)

    # Instrument strength on every row; flagged, never filtered (see WEAK_F).
    with np.errstate(invalid="ignore", divide="ignore"):
        d["F_instrument"] = (d["b_eQTL"].astype(float) / d["se_eQTL"].astype(float)) ** 2
    d["weak_instrument"] = d["F_instrument"] < WEAK_F

    # Probe-position fallbacks: QC-failed until inspected, so carried as a column.
    d["tss_fallback"] = d["gene"].astype(str).isin(tss_fallback_genes(targets))

    d["smr_status"] = [smr_status(r.p_SMR, r.p_HEIDI, r.nsnp_HEIDI, r.smr_threshold)
                       for r in d.itertuples(index=False)]
    if multi:
        # `smr_status` stays the single-SNP verdict in this arm too, so the two arms are
        # comparable row for row; the multi-SNP verdict is a separate column.
        d["smr_multi_status"] = [
            smr_multi_status(r.p_SMR, getattr(r, SMR_MULTI_COL), r.p_HEIDI, r.nsnp_HEIDI,
                             r.smr_multi_threshold)
            for r in d.itertuples(index=False)]
    d["agreement"] = [classify(r.coloc_PP4_probe, r.p_SMR, r.p_HEIDI, r.nsnp_HEIDI,
                               r.smr_threshold)
                      for r in d.itertuples(index=False)]

    attrition = collect_attrition(dest, root)
    if not attrition.empty:
        attrition.to_parquet(dest / "snp_attrition.parquet", index=False)

    d.to_parquet(dest / "smr_results.parquet", index=False)
    pd.DataFrame(runs).to_csv(dest / "smr_runs.tsv", sep="\t", index=False)
    _write_report(dest, d, qtl_source, peqtl_smr, arm, attrition, multi)
    print(f"  {len(d):,} probe rows; agreement: {d['agreement'].value_counts().to_dict()}")
    if multi:
        print(f"  multi-SNP: {d['smr_multi_status'].value_counts().to_dict()}")
    print(f"  wrote {dest}")
    return dest


def _fmt(v, nd=3) -> str:
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "—"
    if isinstance(v, (float, np.floating)):
        return f"{v:.{nd}g}" if (abs(v) < 1e-3 and v != 0) else f"{v:.{nd}f}"
    return str(v)


def _write_report(dest: Path, d: pd.DataFrame, qtl_source: str, peqtl_smr: float,
                  arm: str | None = None, attrition: pd.DataFrame | None = None,
                  multi: bool = False) -> None:
    brainseq = qtl_source == "brainseq"
    L: list[str] = []
    A = L.append
    A("# SMR + HEIDI on the signal-level colocalization nominations")
    A("")
    A(f"QTL source: `{qtl_source}`" + (f" (`{arm}` arm)" if brainseq else "")
      + f". Instrument threshold `--peqtl-smr {peqtl_smr:g}`"
      + ("" if np.isclose(peqtl_smr, PEQTL_SMR) else " -- **SENSITIVITY ARM**")
      + (f". Multi-SNP SMR (`--smr-multi`, LD pruned at r2 {LD_MULTI_SNP:g}) "
         "-- **SENSITIVITY ARM**" if multi else "") + ". "
      f"HEIDI SNPs at p < {PEQTL_HEIDI:g}, {HEIDI_MIN_M}-{HEIDI_MAX_M} SNPs, "
      f"{CIS_WIND_KB} kb cis window. Significance: Bonferroni at {SMR_ALPHA} over the "
      f"instrumented probes of the row's OWN family (primary confirmatory / secondary "
      f"event-localization), never pooled across the two. HEIDI rejects a single shared "
      f"variant at p < {HEIDI_REJECT}. LD reference: 1000G EUR.")
    A("")
    A("## How these numbers may be read")
    A("")
    A("- `b_SMR` is a signed ratio estimate relating genetically predicted phenotype to disease "
      "under a single-causal-variant, no-pleiotropy model. It does **not** establish causal "
      "direction and cannot distinguish causality from horizontal pleiotropy.")
    if brainseq:
        A("- BrainSEQ is the discovery cohort: agreement here is same-tissue genetic anchoring, "
          "not replication. Each gene has one `A_g` (abundance) and one `S_g` (switch) probe, "
          "compared with BrainSEQ's own coloc posterior for the same region and axis.")
        A("- `S_g` is PC1 of the gene's within-gene composition. Its `b_SMR` is sign-pinned to "
          "the discovery switch axis, so its sign is a direction along that axis, not \"more "
          "splicing\"; a probe with `sign_pinned = False` keeps its mapped orientation, and "
          "`axis_unstable` marks a gene whose recomputed axis differs from discovery. An `S_g` "
          "result names no intron or event.")
        A("- Failing to reject HEIDI is **not** evidence of a shared variant; at 169-229 donors "
          "per region HEIDI is weaker still than at GTEx. It reads \"not rejected\".")
    else:
        A("- An sQTL probe is a LeafCutter intron-excision ratio, which is compositional within "
          "its cluster: introns sharing a splice site trade usage, so sibling probes carry "
          "opposite `b_SMR` signs by construction. A sign is read relative to its cluster, "
          "never alone.")
        A("- Failing to reject HEIDI is **not** evidence of a shared variant; HEIDI is "
          "underpowered at GTEx brain sample sizes. It reads \"not rejected\".")
    A("- A HEIDI rejection does **not** overrule a strong signal-level colocalization, and SMR "
      "significance does not promote a locus coloc did not support. Disagreements stay "
      "disagreements.")
    A("")
    prim = d[d["primary_probe"]]
    A("## What was testable, by family")
    A("")
    A("Two families, corrected apart. The **primary confirmatory family** holds one "
      "pre-designated probe per gene; a gene's other introns form a **secondary "
      "event-localization family** with its own Bonferroni correction, so the primary "
      "threshold is not inflated by introns and the introns do not escape correction when they "
      "are discussed. `F` is the instrument strength `(b_eQTL/se_eQTL)^2` of the top cis-QTL "
      "SNP, reported on every row and never used to exclude one. It is summarized by its 5th "
      "percentile rather than a count below the conventional F < 10: an instrument that clears "
      "p < 5e-8 has |z| ≳ 5.4 and so F ≳ 30 (≈ 24 at the relaxed 1e-6 arm), so a weak-instrument "
      "count is zero by construction and carries no information; the per-row `weak_instrument` "
      "flag stays in `smr_results.parquet` for any run at a looser threshold.")
    A("")
    A("| analysis | modality | family | probes | instrumented | threshold | "
      "F median [min-max] | F p5 | "
      + " | ".join(f"`{s}`" for s in SMR_STATUS) + " |")
    A("|---|---|---|---|---|---|---|---|" + "---|" * len(SMR_STATUS))
    for (an, mod, fam), sub in d.groupby(["analysis", "modality", "probe_family"]):
        cnt = sub["smr_status"].value_counts()
        ins = sub[sub["p_SMR"].notna()]
        f = ins["F_instrument"].dropna()
        frange = ("—" if f.empty else
                  f"{f.median():.0f} [{f.min():.0f}-{f.max():.0f}]")
        fp5 = "—" if f.empty else f"{f.quantile(0.05):.0f}"
        A(f"| {an} | {mod} | {fam} | {len(sub)} | {len(ins)} | "
          f"{_fmt(sub['smr_threshold'].iloc[0])} | {frange} | {fp5} | "
          + " | ".join(str(int(cnt.get(s, 0))) for s in SMR_STATUS) + " |")
    A("")
    A("**`no_instrument` is not a negative result** — the probe was never tested, because no "
      "cis-QTL reached the instrument threshold. It is the largest cell in every arm here and "
      "must never be read as evidence against a locus.")
    A("")
    A("## Agreement with coloc, primary confirmatory family")
    A("")
    A("| analysis | modality | probes | instrumented | threshold | "
      + " | ".join(f"`{a}`" for a in AGREEMENT) + " |")
    A("|---|---|---|---|---|" + "---|" * len(AGREEMENT))
    for (an, mod), sub in prim.groupby(["analysis", "modality"]):
        cnt = sub["agreement"].value_counts()
        A(f"| {an} | {mod} | {len(sub)} | {int(sub['p_SMR'].notna().sum())} | "
          f"{_fmt(sub['smr_threshold'].iloc[0])} | "
          + " | ".join(str(int(cnt.get(a, 0))) for a in AGREEMENT) + " |")
    A("")
    sec = d[~d["primary_probe"]]
    if len(sec):
        A("### Secondary family (event localization)")
        A("")
        A("A gene's non-primary introns, corrected within their own family. These localize an "
          "event; they are not additional confirmatory evidence for a locus.")
        A("")
        A("| analysis | modality | probes | instrumented | threshold | `smr_heidi_supported` | "
          "`smr_signal_heidi_rejects` |")
        A("|---|---|---|---|---|---|---|")
        for (an, mod), sub in sec.groupby(["analysis", "modality"]):
            cnt = sub["smr_status"].value_counts()
            A(f"| {an} | {mod} | {len(sub)} | {int(sub['p_SMR'].notna().sum())} | "
              f"{_fmt(sub['smr_threshold'].iloc[0])} | "
              f"{int(cnt.get('smr_heidi_supported', 0))} | "
              f"{int(cnt.get('smr_signal_heidi_rejects', 0))} |")
        A("")
    if attrition is not None and not attrition.empty and "frac_panel_shared" in attrition:
        A("## SNP attrition into HEIDI")
        A("")
        A("HEIDI is computed on the SNPs shared by the BESD, the LD reference and the GWAS, so "
          "what survives is worth stating. Per SMR run; the per-probe end of the trace is "
          f"`nsnp_HEIDI` in the results table, itself capped at {HEIDI_MAX_M}.")
        A("")
        A("**The dominant filter is GWAS locus coverage, by construction** — the `.ma` files "
          "carry only the SNPs of the selected GWAS loci, so a gene whose cis window extends "
          "past its locus loses the remainder. That is not a QC failure and says nothing about "
          "the QTL data.")
        A("")
        A("| step | median across runs |")
        A("|---|---|")
        A(f"| BESD SNPs in the probe's cis window | {attrition['n_besd_snps_counted'].median():,.0f} |")
        A(f"| ... found in the LD panel by rsID | {attrition['frac_panel_shared'].median():.1%} |")
        A(f"| ... with matching alleles | {attrition['frac_allele_match'].median():.1%} |")
        A(f"| ... also carried by the GWAS | {attrition['n_gwas_shared'].median():,.0f} |")
        A(f"| dropped by SMR beyond that (frequency check, duplicates) | "
          f"{attrition['n_harmonization_residual'].median():,.0f} |")
        A("")
        worst = int(attrition["n_harmonization_residual"].max())
        A(f"Genuine harmonization loss — SNPs present in all three and still dropped — peaks at "
          f"**{worst}** SNPs in any run. The identifier and allele steps are the QC ones, and "
          "both sit at 100%.")
        A("")
        A("**Do not compare BESD-SNP retention across QTL sources.** Its denominator is "
          "imputation density: BrainSEQ's TOPMed cis windows carry several times the SNPs of "
          "GTEx's pre-filtered all-pairs, so BrainSEQ retains a far smaller *fraction* of a far "
          "larger set against the same GWAS loci. The ratio measures panel density, not quality.")
        A("")
    nfb = int(d["tss_fallback"].sum())
    if nfb:
        genes = sorted(set(d.loc[d["tss_fallback"], "symbol"].astype(str)))
        A("## Probe-position fallbacks (QC)")
        A("")
        A(f"{nfb:,} probe rows over {len(genes)} genes had no hg19 TSS and took the median-ESD "
          "position instead, which moves the 2 Mb cis window off the promoter and can change "
          "which SNPs are eligible as instruments. **Treat these as QC-failed until inspected**: "
          + ", ".join(f"`{g}`" for g in genes[:20])
          + (" ..." if len(genes) > 20 else "") + ".")
        A("")
    if brainseq and "target_source" in prim.columns:
        A("## GTEx nominations on the BrainSEQ switch axis")
        A("")
        A("Per GTEx signal-level nomination, the `S_g` agreement class in each region "
          "(`*` = the region's tissue-matched GTEx tissue carries the GTEx sQTL call).")
        A("")
        s = prim[(prim["modality"] == "S_g") & prim["target_source"].isin(["gtex", "both"])]
        regions = sorted(prim["tissue"].unique())
        A("| gene | trait | " + " | ".join(regions) + " |")
        A("|---|---|" + "---|" * len(regions))
        for (sym, trait), sub in s.groupby(["symbol", "trait"], sort=True):
            cells = []
            for rg in regions:
                r = sub[sub["tissue"] == rg]
                if r.empty:
                    cells.append("—")
                    continue
                r = r.iloc[0]
                cells.append(f"`{r.agreement}`" + ("*" if bool(r.gtex_call_in_matched_tissue) else ""))
            A(f"| {sym} | {trait} | " + " | ".join(cells) + " |")
        A("")
    A("## Primary probes")
    A("")
    if brainseq:
        A("| gene | trait | region | axis | coloc PP4 | b_SMR (se) | p_SMR | p_HEIDI (nsnp) | "
          "pin | agreement |")
        A("|---|---|---|---|---|---|---|---|---|---|")
    else:
        A("| gene | trait | tissue | modality | coloc PP4 | b_SMR (se) | p_SMR | p_HEIDI (nsnp) | "
          "agreement |")
        A("|---|---|---|---|---|---|---|---|---|")
    show = prim if not brainseq else prim[prim["p_SMR"].notna()
                                          | (prim["coloc_PP4_probe"] >= PP4_CALL)]
    for r in show.sort_values(["trait", "symbol", "tissue", "modality"]).itertuples(index=False):
        est = r.coloc_estimator_probe if pd.notna(r.coloc_estimator_probe) else "—"
        row = (f"| {r.symbol} | {r.trait} | {r.tissue} | {r.modality} | "
               f"{_fmt(r.coloc_PP4_probe)} ({est}) | {_fmt(r.b_SMR)} ({_fmt(r.se_SMR)}) | "
               f"{_fmt(r.p_SMR)} | {_fmt(r.p_HEIDI)} ({_fmt(r.nsnp_HEIDI, 0)}) | ")
        if brainseq:
            if r.modality != "S_g":
                pin = ""
            elif pd.isna(r.sign_pinned):
                pin = "no pin"
            else:
                pin = ("pinned" if bool(r.sign_pinned) else "unpinned") + (
                    ", axis unstable" if bool(r.axis_unstable) else "")
            row += f"{pin} | "
        A(row + f"`{r.agreement}` |")
    A("")
    if multi:
        A("## The multi-SNP arm")
        A("")
        A("`--smr-multi` combines the cis SNPs surviving LD pruning at r2 "
          f"{LD_MULTI_SNP:g} instead of testing the top SNP alone. It is a robustness arm "
          "for the SMR estimate, not a second discovery pass: it answers whether a signal "
          "rests on one lead SNP or on the cis signal more broadly.")
        A("")
        A(f"`p_SMR` and `smr_status` in this table are still the single-SNP quantities, so "
          f"rows match the primary arm one for one. The multi-SNP verdict is "
          f"`{SMR_MULTI_COL}` / `smr_multi_status`, Bonferroni-corrected over the probes the "
          "multi-SNP test actually ran on (`n_multi_family`) -- a smaller denominator than "
          "the instrumented count, because SMR skips the test where too few cis SNPs survive "
          "pruning. Those probes are `multi_unavailable`, which is **not** a null result.")
        A("")
        if "smr_multi_status" in d.columns:
            vc = d["smr_multi_status"].value_counts()
            A("| smr_multi_status | probes |")
            A("|---|---|")
            for k in SMR_MULTI_STATUS:
                if k in vc.index:
                    A(f"| `{k}` | {int(vc[k]):,} |")
            A("")
            both = d[d["p_SMR"].notna() & d["smr_threshold"].notna()]
            sig1 = both[both["p_SMR"] <= both["smr_threshold"]]
            if len(sig1):
                surv = sig1[sig1["smr_multi_status"].isin(
                    ("multi_signal_heidi_unavailable", "multi_signal_heidi_rejects",
                     "multi_heidi_supported"))]
                unav = sig1[sig1["smr_multi_status"] == "multi_unavailable"]
                A(f"Of the {len(sig1):,} probes significant on the single-SNP test, "
                  f"{len(surv):,} are also significant under the multi-SNP test and "
                  f"{len(unav):,} could not be tested. A probe that does not survive is a "
                  "signal carried by its lead SNP alone; that is a caveat on the SMR "
                  "estimate, not a refutation of the colocalization.")
            A("")
    if brainseq:
        A("Probes with neither an SMR instrument nor a coloc call are in `smr_results.parquet` "
          "and not listed here.")
    else:
        A("Non-primary sQTL probes (the gene's other introns) are in `smr_results.parquet` with "
          "`primary_probe = False`; they are reported so the SMR evidence is not a maximum over "
          "introns, and they are not summarized here.")
    A("")
    (dest / "SMR_HEIDI.md").write_text("\n".join(L) + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("prep", "besd", "smr", "meta"), required=True)
    ap.add_argument("--qtl-source", choices=QTL_SOURCES, default="gtex")
    ap.add_argument("--arm", default=None, help="BrainSEQ genotype arm (must be ea_only)")
    ap.add_argument("--task", type=int, default=None,
                    help="0-based row of work_list.tsv (besd/smr); default: every row")
    ap.add_argument("--peqtl-smr", type=float, default=PEQTL_SMR,
                    help="instrument threshold; anything but 5e-8 is a sensitivity arm")
    ap.add_argument("--smr-multi", action="store_true",
                    help="multi-SNP SMR sensitivity arm; writes under sensitivity/smr_multi/")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args(argv)
    if args.stage == "prep":
        run_prep(args.qtl_source, arm=args.arm)
    elif args.stage == "besd":
        run_besd(args.qtl_source, task=args.task, arm=args.arm)
    elif args.stage == "smr":
        run_smr(args.qtl_source, task=args.task, peqtl_smr=args.peqtl_smr,
                threads=args.threads, arm=args.arm, multi=args.smr_multi)
    else:
        run_meta(args.qtl_source, peqtl_smr=args.peqtl_smr, arm=args.arm,
                 multi=args.smr_multi)


if __name__ == "__main__":
    main()
