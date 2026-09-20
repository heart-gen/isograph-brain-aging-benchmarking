"""One table for every gene whose colocalizing junction resolves to an IsoGraph switch.

The concordant set -- 30 genes over 76 events on the 2026-09-19 switching-filter re-run --
had no single place to read its evidence from. Genetics was spread over four tables with
different keys (coloc.abf per gene x trait, eCAVIAR CLPP per event, SMR per probe, the
concordance layer per junction) and orthogonal validation over three more with different
units (whole transcripts for long-read, PSI events for the LIBD arm, switch pairs for the
allele-aware junction recount). Comparing two genes meant joining seven files by hand, and
the summaries written that way went stale silently.

**coloc is primary, CLPP is secondary.** eCAVIAR CLPP is a per-variant colocalization
posterior product; it is noisy at this scale and, for the two largest values in the set,
contradicted by everything else. TMEM175's CLPP of 0.483 sits beside a coloc.abf PP4_sQTL
of ~0 in every trait and 6 of 36 SMR probes supported. So the table sorts on coloc and
carries CLPP as a trailing column, not the other way round, and `genetics_call` never reads
CLPP at all.

What the columns mean
---------------------
``coloc`` (primary)
    ``pp4_sqtl`` / ``pp4_eqtl`` are the best coloc.abf posteriors over GTEx brain tissues
    for that gene x trait, and ``n_tissue_sqtl_coloc`` / ``n_tissue_eqtl_coloc`` how many
    tissues actually reach a coloc call. Splicing-led means a high sQTL posterior AND more
    tissues colocalizing on splicing than on expression -- the tissue counts matter because
    GTEx eQTL power exceeds sQTL power everywhere, so a single-tissue sQTL hit is weak
    evidence on its own.

``smr`` (primary, causal-direction)
    Counts of instrumented probes that SMR supports with HEIDI not rejecting. A gene with
    sQTL support and NO eQTL instrument is the cleanest splicing-led case; HEIDI rejections
    are reported separately because they mean linkage rather than a shared causal variant.

``clpp`` (secondary)
    Max eCAVIAR CLPP over the gene's concordant events. Reported for continuity with the
    earlier ledger; it does not enter `genetics_call`.

``orthogonal``
    Three assays that do NOT share a failure mode. Long-read asks whether the switching
    transcripts exist and trade off on an independent platform, lab and quantifier (ONT
    DLPFC BA9/46, n = 12, Bambu). The allele-aware junction recount asks whether the exact
    junction is used, in BrainSEQ BAMs at ~n = 500. The LIBD PSI arm asks the same question
    through a curated event catalogue, which is the narrowest of the three -- a
    ``junction_not_measured`` there is a gap in that catalogue, not a negative.

Neither call is a test. ``genetics_call`` and ``orthogonal_call`` are transparent
summaries of columns that are all present beside them, so a reader can disagree with the
thresholds without re-deriving the evidence.

Writes ``anchored_gene_summary.{parquet,tsv}`` and ``ANCHORED_GENE_SUMMARY.md`` to
``08_integration/_m/anchored_gene_summary/``.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import subprocess

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, stage_out

# A gene is called splicing-led on coloc when the sQTL posterior is strong AND splicing
# colocalizes in more tissues than expression does. Both halves are needed: GTEx eQTL power
# exceeds sQTL power in every brain tissue, so a lone high PP4_sQTL is not yet a contrast.
PP4_STRONG = 0.8

# Minor-form usage is a DESCRIPTOR here, never a pass/fail gate.
#
# BrainSEQ and GTEx are neurotypical / control-enriched tissue. An isoform that matters in
# disease is often a minority form in normal brain precisely because the disease is what
# raises it, so a low minor-form fraction in controls is a statement about the tissue
# sampled, not about whether the switch is real or relevant. Gating on it would
# preferentially discard exactly the disease-relevant switches this layer exists to find.
# (PI, 2026-09-20.)
#
# What does matter is DETECTION: are both forms actually observed, in how many donors, at
# what depth. That is a power statement and it is reported as counts below. The
# pre-registered 0.05 threshold in junction_coloc_confirm was set against a different
# failure -- SNCA's anchored isoform at 0.29% of ONT gene output, where the form was
# essentially absent, not merely rare -- and is left in that CLI as the pre-registered
# rule. It is deliberately NOT applied here.
MIN_USAGE_REFERENCE = 0.05

# Donors that must carry BOTH forms before a pair counts as observed. A switch seen in a
# handful of donors is not evidence of a switch; one seen in hundreds at 4% is.
MIN_DONORS_BOTH = 30

# Fragments a donor must carry across the pair before its usage ratio is computed at all.
MIN_PAIR_FRAGS = 10

BRAINSEQ_REGIONS = ("caudate", "dlpfc", "hippocampus")


def out_dir():
    return ensure_dir(stage_out("integration", "anchored_gene_summary"))


def _bare(s: pd.Series) -> pd.Series:
    return s.astype(str).str.split(".").str[0]


# --------------------------------------------------------------------------- #
# Sources
# --------------------------------------------------------------------------- #
def concordant_events() -> pd.DataFrame:
    """The gene x trait spine: every event whose junction is in the switch pair."""
    ev = pd.read_parquet(
        stage_out("anchoring.coloc", "coloc_isoform_events_combined.parquet"))
    c = ev[ev["concordant"].eq(True)].copy()
    c["ens"] = _bare(c["gene"])
    return c


def _spine(c: pd.DataFrame) -> pd.DataFrame:
    """One row per (gene, trait), carrying the concordance evidence itself."""
    rows = []
    for (sym, ens, trait), g in c.groupby(["gene_name", "ens", "trait"]):
        rows.append({
            "gene": sym, "ens": ens, "trait": trait,
            "tissue": g["tissue"].mode().iat[0].replace("Brain_", ""),
            "n_tissues": int(g["tissue"].nunique()),
            "n_events": int(len(g)),
            "n_junctions": int(g["junction"].nunique()),
            "go_invisible": bool(g["go_invisible"].eq(True).any()),
            "polarity_r": float(g["junction_transcript_polarity_r"].abs().max()),
            "brainseq_switch_replicates": bool(
                g["brainseq_replicates_switch"].eq(True).any()),
            "clpp": float(g["clpp"].max()),
        })
    return pd.DataFrame(rows)


def _coloc(spine: pd.DataFrame) -> pd.DataFrame:
    """coloc.abf posteriors and tissue counts, per gene x trait."""
    g = pd.read_parquet(stage_out("anchoring.coloc_signal", "genes.parquet"))
    g["ens"] = _bare(g["gene"])
    agg = g.groupby(["ens", "trait"]).agg(
        pp4_sqtl=("PP4_sQTL", "max"), pp4_eqtl=("PP4_eQTL", "max"),
        median_pp4_sqtl=("median_PP4_sQTL", "max"),
        n_tissue_sqtl_coloc=("n_tissue_sQTL_coloc", "max"),
        n_tissue_eqtl_coloc=("n_tissue_eQTL_coloc", "max"),
    ).reset_index()
    return spine.merge(agg, on=["ens", "trait"], how="left")


def _smr(spine: pd.DataFrame) -> pd.DataFrame:
    """Instrumented SMR probes supported / rejected, split by modality."""
    s = pd.read_parquet(stage_out("anchoring.smr", "gtex", "smr_results.parquet"))
    s = s[s["smr_status"].ne("no_instrument")].copy()
    s["ens"] = _bare(s["gene"])
    out = spine.copy()
    for mod, tag in (("sQTL", "sqtl"), ("eQTL", "eqtl")):
        m = s[s["modality"].eq(mod)]
        agg = m.groupby(["ens", "trait"]).agg(
            **{f"smr_{tag}_tested": ("smr_status", "size"),
               f"smr_{tag}_supported": (
                   "smr_status", lambda x: int((x == "smr_heidi_supported").sum())),
               f"smr_{tag}_heidi_rejected": (
                   "smr_status", lambda x: int((x == "smr_signal_heidi_rejects").sum()))},
        ).reset_index()
        out = out.merge(agg, on=["ens", "trait"], how="left")
    cols = [c for c in out.columns if c.startswith("smr_")]
    out[cols] = out[cols].fillna(0).astype(int)
    return out


def _longread(spine: pd.DataFrame) -> pd.DataFrame:
    """ONT long-read: are the switching transcripts there, and do they trade off?"""
    base = stage_out("mechanism", "longread_switch_confirm")
    gc = pd.read_parquet(base / "gene_confirmation.parquet")
    pc = pd.read_parquet(base / "pair_confirmation.parquet")
    g = gc.groupby("gene").agg(
        lr_expressed=("expressed_longread", lambda x: bool(x.eq(True).any())),
        lr_confirmed=("confirmed", lambda x: bool(x.eq(True).any())),
        lr_top2_rho=("top2_usage_spearman", "min"),
    ).reset_index().rename(columns={"gene": "ens"})
    p = pc.groupby("gene").agg(
        lr_pairs_detected=("pair_detected", lambda x: int(x.eq(True).sum())),
        lr_pairs_switchlike=("switch_like", lambda x: int(x.eq(True).sum())),
    ).reset_index().rename(columns={"gene": "ens"})
    out = spine.merge(g, on="ens", how="left").merge(p, on="ens", how="left")
    for c in ("lr_expressed", "lr_confirmed"):
        out[c] = out[c].eq(True)
    for c in ("lr_pairs_detected", "lr_pairs_switchlike"):
        out[c] = out[c].fillna(0).astype(int)
    return out


def _ase(spine: pd.DataFrame) -> pd.DataFrame:
    """BrainSEQ allele-aware junction recount: is the junction actually used?

    Usage is computed from the pair's own junction fragments summed over haplotypes, so it
    is the switch ratio directly. Donors below ``MIN_PAIR_FRAGS`` across the pair are
    dropped before the median rather than contributing a ratio built on a handful of reads.
    """
    want = set(spine["ens"])
    feas, usage = [], []
    for region in BRAINSEQ_REGIONS:
        base = stage_out("mechanism", "ase_junction_switch", region)
        fp = base / "pair_feasibility.parquet"
        if not fp.exists():
            continue
        f = pd.read_parquet(fp)
        f["ens"] = _bare(f["gene_id"])
        f = f[f["ens"].isin(want)]
        if f.empty:
            continue
        f["region"] = region
        feas.append(f)

        pairs = set(f.loc[f["testable"].eq(True), "pair_id"])
        cp = base / "junction_allelic_counts.parquet"
        if not pairs or not cp.exists():
            continue
        cnt = pd.read_parquet(cp, columns=["donor_id", "pair_id", "isoform", "n_frag"])
        cnt = cnt[cnt["pair_id"].isin(pairs)]
        if cnt.empty:
            continue
        tot = (cnt.groupby(["pair_id", "donor_id", "isoform"])["n_frag"].sum()
               .unstack(fill_value=0).reindex(columns=[1, 2], fill_value=0))
        tot = tot[tot.sum(axis=1) >= MIN_PAIR_FRAGS]
        if tot.empty:
            continue
        frac = tot[1] / tot.sum(axis=1)
        both = (tot[1] > 0) & (tot[2] > 0)
        per_pair = pd.DataFrame({
            "med": frac.groupby(level="pair_id").median(),
            "n_donors": frac.groupby(level="pair_id").size(),
            "n_both": both.groupby(level="pair_id").sum(),
        })
        per_pair["frac_both"] = per_pair["n_both"] / per_pair["n_donors"]
        per_pair["minor"] = np.minimum(per_pair["med"], 1 - per_pair["med"])
        per_pair = per_pair.reset_index()
        per_pair["region"] = region
        per_pair["ens"] = per_pair["pair_id"].map(
            dict(zip(f["pair_id"], f["ens"], strict=True)))
        usage.append(per_pair)

    out = spine.copy()
    if feas:
        fa = pd.concat(feas, ignore_index=True)
        agg = fa.groupby("ens").agg(
            ase_pairs=("pair_id", "nunique"),
            ase_testable=("testable", lambda x: int(x.eq(True).sum())),
            ase_donors=("n_donors", "max"),
        ).reset_index()
        out = out.merge(agg, on="ens", how="left")
    if usage:
        ua = pd.concat(usage, ignore_index=True)
        # Select on DETECTION, not on usage: a gene with several switch pairs is observed
        # if any pair shows both forms in many donors. Selecting the highest minor-form
        # usage instead would rank a rare-but-real disease isoform below a balanced
        # housekeeping one, which is the bias this table is built to avoid.
        best = ua.sort_values(["n_both", "n_donors"], ascending=False).drop_duplicates("ens")
        out = out.merge(
            best[["ens", "minor", "n_donors", "n_both", "frac_both", "region"]].rename(
                columns={"minor": "ase_minor_usage", "n_donors": "ase_usage_donors",
                         "n_both": "ase_donors_both_forms",
                         "frac_both": "ase_frac_donors_both",
                         "region": "ase_region"}),
            on="ens", how="left")
    for c in ("ase_pairs", "ase_testable", "ase_donors", "ase_usage_donors",
              "ase_donors_both_forms"):
        if c in out:
            out[c] = out[c].fillna(0).astype(int)
        else:
            out[c] = 0
    for c in ("ase_minor_usage", "ase_frac_donors_both"):
        if c not in out:
            out[c] = np.nan
    if "ase_region" not in out:
        out["ase_region"] = ""
    return out


def _psi(spine: pd.DataFrame) -> pd.DataFrame:
    """The LIBD PSI arm's verdict, as the narrowest of the three orthogonal assays."""
    f = stage_out("mechanism", "junction_coloc_confirm", "junction_confirm.parquet")
    if not f.exists():
        return spine.assign(psi_verdict="not run")
    r = pd.read_parquet(f)

    def call(g: pd.DataFrame) -> str:
        v = set(g["verdict"])
        if "validated" in v:
            return "validated"
        if "not_validated" in v:
            return "not validated"
        if "junction_not_measured" in v:
            return "junction not in catalogue"
        return "no BrainSEQ region"

    agg = r.groupby(["gene_name", "trait"], group_keys=False).apply(
        call, include_groups=False).rename("psi_verdict").reset_index()
    agg = agg.rename(columns={"gene_name": "gene"})
    return spine.merge(agg, on=["gene", "trait"], how="left")


# --------------------------------------------------------------------------- #
# Calls
# --------------------------------------------------------------------------- #
def _genetics_call(r: pd.Series) -> str:
    """Transparent summary of the coloc + SMR columns. CLPP is deliberately not read."""
    strong = (r["pp4_sqtl"] >= PP4_STRONG) if pd.notna(r["pp4_sqtl"]) else False
    more_splicing = r["n_tissue_sqtl_coloc"] > r["n_tissue_eqtl_coloc"]
    smr_s = r["smr_sqtl_supported"] > 0
    no_eqtl_instrument = r["smr_eqtl_tested"] == 0
    if strong and more_splicing and smr_s and no_eqtl_instrument:
        return "splicing-specific"
    if strong and more_splicing and smr_s:
        return "splicing-led"
    if strong and more_splicing:
        return "splicing-led (coloc only)"
    if strong:
        return "colocalizes, not splicing-specific"
    return "weak coloc"


def _orthogonal_call(r: pd.Series) -> str:
    """How many independent assays OBSERVE the switch.

    Observation, not abundance: the junction-recount arm counts as confirming when both
    forms are seen in enough donors, whatever the ratio between them. A 4% minority form
    detected in 495 of 498 neurotypical donors is a well-measured switch, and may be
    exactly the form a disease raises.
    """
    lr = bool(r["lr_confirmed"]) and r["lr_pairs_switchlike"] > 0
    ase = r["ase_donors_both_forms"] >= MIN_DONORS_BOTH
    psi = r.get("psi_verdict") == "validated"
    n = int(lr) + int(ase) + int(psi)
    if n >= 2:
        return f"confirmed ({n} assays)"
    if n == 1:
        return "confirmed (1 assay)"
    if bool(r["lr_expressed"]) or r["ase_testable"] > 0:
        return "tested, not confirmed"
    return "not measurable"


# --------------------------------------------------------------------------- #
# Provenance
# --------------------------------------------------------------------------- #
def _inputs() -> list:
    """Every file this table is derived from, in read order."""
    return [
        stage_out("anchoring.coloc", "coloc_isoform_events_combined.parquet"),
        stage_out("anchoring.coloc_signal", "genes.parquet"),
        stage_out("anchoring.smr", "gtex", "smr_results.parquet"),
        stage_out("mechanism", "longread_switch_confirm", "gene_confirmation.parquet"),
        stage_out("mechanism", "longread_switch_confirm", "pair_confirmation.parquet"),
        stage_out("mechanism", "junction_coloc_confirm", "junction_confirm.parquet"),
        *[stage_out("mechanism", "ase_junction_switch", r, n)
          for r in BRAINSEQ_REGIONS
          for n in ("pair_feasibility.parquet", "junction_allelic_counts.parquet")],
    ]


def _digest(path, cap: int = 64 << 20) -> str:
    """SHA-256 of the file, or of its first ``cap`` bytes for the large count tables.

    A truncated digest is marked as such rather than silently passed off as a whole-file
    hash: junction_allelic_counts is ~16.5M rows per region and hashing it in full on every
    run would cost more than the table it verifies.
    """
    h, n = hashlib.sha256(), 0
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            if n + len(chunk) > cap:
                h.update(chunk[: cap - n])
                return f"sha256-first{cap}:{h.hexdigest()}"
            h.update(chunk)
            n += len(chunk)
    return f"sha256:{h.hexdigest()}"


def provenance() -> dict:
    """What the numbers were built from, so a reader can tell if they are still current."""
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                                text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain"], capture_output=True,
                                    text=True, check=True).stdout.strip())
    except (subprocess.CalledProcessError, FileNotFoundError):
        commit, dirty = "unknown", False
    files = []
    for f in _inputs():
        if f.exists():
            st = f.stat()
            files.append({
                "path": str(f), "bytes": st.st_size,
                "mtime": _dt.datetime.fromtimestamp(
                    st.st_mtime, _dt.timezone.utc).isoformat(timespec="seconds"),
                "digest": _digest(f)})
        else:
            files.append({"path": str(f), "missing": True})
    return {
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "git_commit": commit, "git_dirty": dirty,
        "generator": "isograph_benchmark/real_data/anchored_gene_summary.py",
        "pandas": pd.__version__, "numpy": np.__version__,
        "thresholds": {"PP4_STRONG": PP4_STRONG, "MIN_DONORS_BOTH": MIN_DONORS_BOTH,
                       "MIN_PAIR_FRAGS": MIN_PAIR_FRAGS,
                       "MIN_USAGE_REFERENCE": MIN_USAGE_REFERENCE},
        "inputs": files,
    }


# --------------------------------------------------------------------------- #
# Build
# --------------------------------------------------------------------------- #
def build() -> pd.DataFrame:
    c = concordant_events()
    t = _spine(c)
    t = _coloc(t)
    t = _smr(t)
    t = _longread(t)
    t = _ase(t)
    t = _psi(t)
    for col in ("n_tissue_sqtl_coloc", "n_tissue_eqtl_coloc"):
        t[col] = t[col].fillna(0).astype(int)
    t["genetics_call"] = t.apply(_genetics_call, axis=1)
    t["orthogonal_call"] = t.apply(_orthogonal_call, axis=1)
    # Sort on coloc, never on CLPP: strongest splicing evidence first.
    order = {"splicing-specific": 0, "splicing-led": 1, "splicing-led (coloc only)": 2,
             "colocalizes, not splicing-specific": 3, "weak coloc": 4}
    t["_o"] = t["genetics_call"].map(order)
    t = t.sort_values(["_o", "pp4_sqtl"], ascending=[True, False]).drop(columns="_o")
    return t.reset_index(drop=True)


def _f(v, d=3):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return "--"
    return "--" if not np.isfinite(v) else f"{v:.{d}f}"


def report(t: pd.DataFrame, prov: dict | None = None) -> str:
    L = [
        "# Genetics and orthogonal validation for the anchored switch genes", "",
        f"**{t['gene'].nunique()}** genes ({len(t)} gene x trait rows) whose colocalizing "
        f"sQTL junction resolves to the tissue-matched IsoGraph switch pair.", "",
        "Sorted on **coloc**, which is the primary genetic evidence. eCAVIAR **CLPP is "
        "secondary** and is carried as a trailing column only: it is noisy per locus and, "
        "for the two largest values in this set, contradicted by every other layer. It "
        "does not enter `genetics_call`.", "",]
    if prov:
        miss = [f["path"] for f in prov["inputs"] if f.get("missing")]
        L += [
            f"Generated {prov['generated_utc']} from commit `{prov['git_commit'][:8]}`"
            + (" **with a dirty working tree -- not reproducible from that commit alone**"
               if prov["git_dirty"] else "")
            + f", pandas {prov['pandas']}, numpy {prov['numpy']}, over "
              f"{len(prov['inputs']) - len(miss)} input files. Digests, sizes and mtimes "
              f"are in `provenance.json`; re-run "
              f"`python -m isograph_benchmark.real_data.anchored_gene_summary` to "
              f"regenerate. Do not edit any number here by hand.",
        ]
        if miss:
            L += ["", f"**{len(miss)} input(s) missing at build time**: "
                      + ", ".join(miss) + ". Columns derived from them are blank."]
        L += [""]
    L += [
        "## Genetics", "",
        "`PP4` are coloc.abf posteriors over GTEx brain; `tis` counts tissues reaching a "
        "coloc call, which matters because GTEx eQTL power exceeds sQTL power everywhere. "
        "`SMR` is supported/tested instrumented probes, `!n` marks HEIDI rejections. A "
        "gene with sQTL support and no eQTL instrument at all is the cleanest "
        "splicing-specific case.", "",
        "| gene | trait | tissue | PP4 sQTL | PP4 eQTL | tis s/e | SMR sQTL | SMR eQTL | "
        "GO-inv | call | CLPP |",
        "| --- | --- | --- | ---: | ---: | :---: | :---: | :---: | :---: | --- | ---: |",
    ]

    def smr_cell(sup, tot, rej):
        if not tot:
            return "--"
        return f"{int(sup)}/{int(tot)}" + (f" !{int(rej)}" if rej else "")

    for _, r in t.iterrows():
        L.append(
            f"| {r['gene']} | {r['trait']} | {r['tissue']} | {_f(r['pp4_sqtl'])} | "
            f"{_f(r['pp4_eqtl'])} | {int(r['n_tissue_sqtl_coloc'])}/"
            f"{int(r['n_tissue_eqtl_coloc'])} | "
            f"{smr_cell(r['smr_sqtl_supported'], r['smr_sqtl_tested'], r['smr_sqtl_heidi_rejected'])} | "
            f"{smr_cell(r['smr_eqtl_supported'], r['smr_eqtl_tested'], r['smr_eqtl_heidi_rejected'])} | "
            f"{'yes' if r['go_invisible'] else 'no'} | {r['genetics_call']} | "
            f"{_f(r['clpp'])} |")

    L += [
        "", "## Orthogonal validation", "",
        "Three assays with different failure modes. **Long-read**: ONT DLPFC BA9/46, "
        "n = 12, Bambu -- do the switching transcripts exist and trade off on an "
        "independent platform? **Junction recount**: allele-aware junction counts from "
        "BrainSEQ BAMs at ~n = 500 -- is the junction itself observed, and in how many "
        "donors? **PSI**: the LIBD event catalogue, the narrowest of the three -- "
        "`junction not in catalogue` is a gap in that resource, not a negative.", "",
        "**Minor-form usage is a descriptor, not a criterion.** BrainSEQ and GTEx are "
        "neurotypical tissue, so an isoform that matters in disease is often a minority "
        "form here precisely because disease is what raises it. `usage` is reported for "
        "interpretation; the call keys on whether both forms are observed and in how many "
        f"donors (`both` >= {MIN_DONORS_BOTH}). "
        f"For reference only, {MIN_USAGE_REFERENCE} is the pre-registered threshold the "
        "separate PSI arm applies.", "",
        "| gene | trait | LR conf | LR pairs det/switch | ASE pairs test | ASE donors both |"
        " ASE usage | PSI | switch replicates | call |",
        "| --- | --- | :---: | :---: | :---: | ---: | ---: | --- | :---: | --- |",
    ]
    for _, r in t.iterrows():
        both = int(r["ase_donors_both_forms"])
        n_d = int(r["ase_usage_donors"])
        L.append(
            f"| {r['gene']} | {r['trait']} | {'yes' if r['lr_confirmed'] else 'no'} | "
            f"{int(r['lr_pairs_detected'])}/{int(r['lr_pairs_switchlike'])} | "
            f"{int(r['ase_testable'])}/{int(r['ase_pairs'])} | "
            f"{both}/{n_d} | "
            f"{_f(r['ase_minor_usage'], 3)} | {r['psi_verdict']} | "
            f"{'yes' if r['brainseq_switch_replicates'] else 'no'} | "
            f"{r['orthogonal_call']} |")

    gc = t["genetics_call"].value_counts()
    oc = t["orthogonal_call"].value_counts()
    L += ["", "## Counts", "", "| genetics call | rows |", "| --- | ---: |"]
    L += [f"| {k} | {int(v)} |" for k, v in gc.items()]
    L += ["", "| orthogonal call | rows |", "| --- | ---: |"]
    L += [f"| {k} | {int(v)} |" for k, v in oc.items()]
    if prov:
        L += ["", "## Reproducibility", "",
              "| input | bytes | mtime (UTC) |", "| --- | ---: | --- |"]
        for f in prov["inputs"]:
            if f.get("missing"):
                L.append(f"| `{f['path']}` | missing | -- |")
            else:
                L.append(f"| `{f['path']}` | {f['bytes']} | {f['mtime']} |")
        L += ["", "Thresholds: "
              + ", ".join(f"`{k}` = {v}" for k, v in prov["thresholds"].items()) + "."]
    return "\n".join(L)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args()
    t = build()
    prov = provenance()
    od = out_dir()
    t.to_parquet(od / "anchored_gene_summary.parquet", index=False)
    t.to_csv(od / "anchored_gene_summary.tsv", sep="\t", index=False)
    (od / "provenance.json").write_text(json.dumps(prov, indent=2))
    (od / "ANCHORED_GENE_SUMMARY.md").write_text(report(t, prov))
    if prov["git_dirty"]:
        print("WARNING: working tree is dirty; these numbers are not reproducible from "
              f"{prov['git_commit'][:8]} alone", flush=True)
    print(f"{t['gene'].nunique()} genes, {len(t)} gene x trait rows -> {od}", flush=True)
    print(t["genetics_call"].value_counts().to_string())
    print()
    print(t["orthogonal_call"].value_counts().to_string())


if __name__ == "__main__":
    main()
