"""Orthogonal validation of IsoGraph switches against junction-derived splicing.

IsoGraph's switch signal is computed from the transcript quantifier (Salmon for
BrainSEQ, RSEM for GTEx). This CLI asks whether the SAME genes / switch events
also show age-associated splicing in an INDEPENDENT, split-read-based
quantification that never touches the transcript quantifier:

  * BrainSEQ -- LIBD PSI events (``inputs/processed/brainseq/<region>/psi_events.parquet``),
    percent-spliced-in values computed directly from junction reads.
  * GTEx     -- STAR junction usage. The raw counts are present
    (``inputs/raw/gtex_v11/counts/..._junctions.gct.gz``) but not yet processed
    into a per-region junction-usage table; the GTEx path therefore consumes a
    ``junction_usage.parquet`` produced by a separate ingestion step (see
    ``--cohort gtex`` below). The statistical core is shared.

Because PSI / junction usage are derived from split reads and not from the
transcript quantifier that feeds IsoGraph, concordance is genuine orthogonal
validation rather than a re-description of the same numbers.

Two analyses, both using the manuscript's df=3 natural-cubic-spline age F-test
(:func:`run_models.spline_age_association`) with the standard aging covariates,
so the model is identical to every other real-data association in the paper:

  A. gene_corroboration  (population, direction-free) -- the headline statistic.
     IsoGraph age-switch genes (switch feature-score spline FDR <= alpha) are
     compared to the background of switch-scored genes that are NOT age-switching:
     are switch genes enriched for an independent junction-level age-splicing
     signal (>= 1 PSI event with within-gene Bonferroni + BH FDR <= alpha)?
     Reported as a Fisher exact odds ratio plus a transcript-count-adjusted
     logistic odds ratio (identifiability is the obvious confounder -- genes with
     more isoforms have more PSI events and noisier switch axes).

  B. event_resolved  (targeted, per-gene) -- for the manuscript colocalization
     switches (``deep_dive_events`` rows whose implicated junction lies in a
     resolved IsoGraph switch pair) each junction is matched to the BrainSEQ PSI
     event that measures it (coordinate overlap) and that specific junction's
     PSI ~ age is tested. Each row reports the junction-PSI age p-value and
     direction next to the gene's IsoGraph switch-score age p-value and direction,
     and flags concordance (both significant, same age direction). Missing matches
     (junction not among the gene's PSI events, or gene not expressed) are kept as
     explicit NaN rows rather than dropped.

Outputs land in ``real_data/_m/switch_validation/<cohort>_<region>_<trait>/``.

Interpretation notes (from the first runs):
  * Analysis A (age) is the headline. In caudate (Phase 3, n=238) IsoGraph age-switch
    genes are ~100x enriched for an independent PSI age-splicing signal (9/28 genes,
    32% vs 0.4% background; Fisher OR ~105, transcript-count-adjusted OR ~153). The
    per-gene MARGINAL switch age-test is conservative, though: the smaller Phase 2
    regions yield too few switch-positive genes to test (hippocampus n_switch_pos=0),
    so replication is better assessed on the module-member switch-gene set (genes in
    age-associated co-switch modules, the manuscript's actual unit) rather than
    marginal per-gene calls -- a planned ``--switch-set module`` refinement.
  * Analysis B tests the disease-anchored colocalization switches. Neither age nor
    case/control Dx is the natural validation for a COLOC switch: colocalization means
    the switch is regulated by the risk VARIANT, not that it differs on average between
    cases and controls (an allele-frequency-diluted signal). Accordingly the Dx/age
    event tests are near-null; what Analysis B does establish is STRUCTURAL
    corroboration -- 9/17 GWAS-implicated switch junctions are independently detected
    as real, used PSI events. The strong genetic replication is the sQTL (genotype ->
    junction-PSI) test, which needs the controlled TOPMed genotypes (follow-on).
"""
from __future__ import annotations

import argparse
import json
import re

import numpy as np
import pandas as pd
from scipy import stats

from isograph.io.artifacts import load_dataset_bundle

from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.run_models import (
    diagnosis_association,
    spline_age_association,
)
from isograph_benchmark.real_data.sweep_leiden import (
    AGING_COVARIATE_COLS,
    GTEX_AGE_COL,
    GTEX_COVARIATE_COLS,
    SCZD_COVARIATE_COLS,
)

_META_COLS = ("feature_id", "gene_id", "feature_type", "n_transcripts")
_JUNC_RE = re.compile(r"(chr[\w]+):(\d+)-(\d+)")
_COORD_MATCH_TOL = 5  # bp tolerance when matching a switch junction to a PSI event

BRAINSEQ_AGING_REGIONS = ("caudate", "hippocampus", "dlpfc")
GTEX_REGIONS = (
    "amygdala", "anterior_cingulate_cortex_ba24", "caudate_basal_ganglia",
    "cerebellar_hemisphere", "cerebellum", "cortex", "frontal_cortex_ba9",
    "hippocampus", "hypothalamus", "nucleus_accumbens_basal_ganglia",
    "putamen_basal_ganglia", "spinal_cord_cervical_c_1", "substantia_nigra",
)


# --------------------------------------------------------------------------- #
# Cohort/trait configuration
# --------------------------------------------------------------------------- #
def _covariates(cohort: str, trait: str) -> list[str]:
    if cohort == "gtex":
        return GTEX_COVARIATE_COLS
    return SCZD_COVARIATE_COLS if trait == "dx" else AGING_COVARIATE_COLS


def _age_col(cohort: str) -> str:
    return GTEX_AGE_COL if cohort == "gtex" else "Age"


# --------------------------------------------------------------------------- #
# Path helpers
# --------------------------------------------------------------------------- #
def _artifact_dir(cohort: str, region: str, variant: str, trait: str) -> "object":
    subdir = "isograph_vae_with_abundance" if variant == "with-abundance" else "isograph_vae"
    if cohort == "gtex":
        return rel("real_data", "gtex", region, "_m", subdir)
    tree = "caudate_sczd" if trait == "dx" else region
    return rel("real_data", "brainseq", tree, "_m", subdir)


def _bundle_path(cohort: str, region: str, trait: str) -> "object":
    if cohort == "gtex":
        return rel("inputs", "bundles", "gtex_v11_brain", region)
    if trait == "dx":
        return rel("inputs", "bundles", "brainseq_sczd", "caudate")
    return rel("inputs", "bundles", "brainseq_v1", region)


def _splice_path(cohort: str, region: str) -> "object":
    if cohort == "gtex":
        return rel("inputs", "processed", "gtex_v11", region, "junction_usage.parquet")
    return rel("inputs", "processed", "brainseq", region, "psi_events.parquet")


def _out_dir(cohort: str, region: str, trait: str) -> "object":
    return ensure_dir(rel("real_data", "_m", "switch_validation", f"{cohort}_{region}_{trait}"))


def _strip_ver(s: pd.Series) -> pd.Series:
    return s.astype(str).str.split(".").str[0]


# --------------------------------------------------------------------------- #
# Loaders
# --------------------------------------------------------------------------- #
def _load_switch_scores(cohort: str, region: str, variant: str, trait: str, sample_ids: list[str]):
    """Return (meta, wide) for IsoGraph's switch channel.

    meta: gene_id, gene, n_transcripts (one row per gene)
    wide: sample_id + one column per gene (versioned gene_id), switch score.
    """
    fs = pd.read_parquet(_artifact_dir(cohort, region, variant, trait) / "feature_scores.parquet")
    sw = fs[fs["feature_type"] == "switch"].copy()
    samp = [c for c in sw.columns if c not in _META_COLS and c in sample_ids]
    meta = sw[["gene_id", "n_transcripts"]].copy()
    meta["gene"] = _strip_ver(meta["gene_id"])
    wide = sw.set_index("gene_id")[samp].T
    wide.index.name = "sample_id"
    wide = wide.reset_index()
    return meta, wide, samp


def _module_switch_genes(cohort: str, region: str, variant: str, trait: str, alpha: float) -> set[str]:
    """Unversioned gene ids in trait-associated co-switch modules.

    This is the manuscript's actual unit (module-level co-switching), and unlike the
    conservative per-gene marginal switch test it has enough power in the smaller
    Phase 2 regions. Modules are called trait-associated at eigengene FDR <= alpha:
    age via the df=3 spline F-test (age_spline.parquet:fdr_ftest), dx via the
    diagnosis t-test (diagnosis_assoc.parquet:fdr).
    """
    art = _artifact_dir(cohort, region, variant, trait)
    modules = pd.read_parquet(art / "modules.parquet")
    if trait == "dx":
        assoc = pd.read_parquet(art / "diagnosis_assoc.parquet")
        sig_col = "fdr"
    else:
        assoc = pd.read_parquet(art / "age_spline.parquet")
        sig_col = "fdr_ftest"
    sig_modules = set(assoc.loc[assoc[sig_col] <= alpha, "module_id"].unique())
    member = modules[modules["module_id"].isin(sig_modules)]
    return set(_strip_ver(member["gene_id"]))


def _load_splice(cohort: str, region: str, sample_ids: list[str], genes: set[str] | None = None):
    """Return (meta, wide) for the independent junction-derived splicing evidence.

    brainseq -> LIBD PSI events (psi_events.parquet), one row per alt-splice event.
    gtex     -> STAR within-gene junction usage (junction_usage.parquet from
                build_gtex_junction_usage), one row per junction. Both are split-read
                based and independent of the transcript quantifier feeding IsoGraph.

    meta: event_id, gene (unversioned), event_type, event_info, chrom, coords (list[int])
    wide: sample_id + one column per event_id (value in [0, 1], NaN where unquantified).
    """
    if cohort == "gtex":
        ju = pd.read_parquet(_splice_path(cohort, region))
        if genes is not None:
            ju = ju[ju["gene"].isin(genes)]
        ju = ju.reset_index(drop=True)
        ju["event_id"] = ju["gene"] + "#" + ju.index.astype(str)
        samp = [c for c in ju.columns if c in sample_ids]
        coords = [[int(a) for pair in _JUNC_RE.findall(str(j)) for a in pair[1:]]
                  for j in ju["junction"]]
        chrom = ju["junction"].astype(str).str.extract(r"(chr[\w]+):")[0].to_numpy()
        meta = pd.DataFrame({
            "event_id": ju["event_id"].to_numpy(), "gene": ju["gene"].to_numpy(),
            "event_type": "junction", "event_info": ju["junction"].to_numpy(),
            "chrom": chrom, "coords": coords,
        })
        wide = ju.set_index("event_id")[samp].T
        wide.index.name = "sample_id"
        return meta, wide.reset_index()

    psi = pd.read_parquet(_splice_path(cohort, region)).copy()
    psi["gene"] = _strip_ver(psi["gene_id"])
    if genes is not None:
        psi = psi[psi["gene"].isin(genes)]
    psi = psi.reset_index(drop=True)
    psi["event_id"] = psi["gene"] + "#" + psi.index.astype(str)

    samp = [c for c in psi.columns if c in sample_ids]
    coords = [
        [int(a) for pair in _JUNC_RE.findall(str(info)) for a in pair[1:]]
        for info in psi["event_info"]
    ]
    meta = pd.DataFrame({
        "event_id": psi["event_id"].to_numpy(),
        "gene": psi["gene"].to_numpy(),
        "event_type": psi["event_type"].to_numpy(),
        "event_info": psi["event_info"].to_numpy(),
        "chrom": psi["chr"].to_numpy(),
        "coords": coords,
    })
    wide = psi.set_index("event_id")[samp].T
    wide.index.name = "sample_id"
    wide = wide.reset_index()
    return meta, wide


# --------------------------------------------------------------------------- #
# Trait association (mirrors every other real-data association)
# --------------------------------------------------------------------------- #
_EMPTY_STATS = pd.DataFrame(columns=["unit", "p", "fdr", "direction"])


def _prune_units(
    wide: pd.DataFrame, min_finite: int = 50, min_var: float = 1e-6, min_unique: int = 10
) -> pd.DataFrame:
    """Drop unit columns that would degenerate the F-test.

    PSI events are frequently all-0/all-1, unquantified in most samples, or take
    only a handful of distinct values; such a response can be fit exactly by the
    covariate + spline design (rss_full ~ 0), giving a non-finite F-statistic and
    a NaN p-value that the downstream FDR cannot handle. Keep sample_id plus only
    columns with enough finite, variable, non-degenerate values.
    """
    units = [c for c in wide.columns if c != "sample_id"]
    vals = wide[units].apply(pd.to_numeric, errors="coerce")
    keep = [
        c for c in units
        if vals[c].notna().sum() >= min_finite
        and float(vals[c].var(skipna=True)) > min_var
        and vals[c].nunique(dropna=True) >= min_unique
    ]
    return wide[["sample_id"] + keep]


def _assoc_stats(wide: pd.DataFrame, sample_table: pd.DataFrame, trait: str,
                 covariate_cols: list[str], age_col: str) -> pd.DataFrame:
    """Per-unit trait association + a signed direction, with uniform columns.

    age -> df=3 natural-spline age F-test (direction = late p75 - early p25 effect).
    dx  -> SCZD-vs-Control diagnosis t-test (direction = case-vs-control effect).
    Returns unit, p, fdr, direction.
    """
    if trait == "dx":
        res = diagnosis_association(wide, sample_table, covariate_cols=covariate_cols)
        if res.empty:
            return _EMPTY_STATS.copy()
        return res[["module_id", "pvalue", "fdr", "effect"]].rename(
            columns={"module_id": "unit", "pvalue": "p", "effect": "direction"})
    res = spline_age_association(
        wide, sample_table, covariate_cols=covariate_cols, age_col=age_col)
    if res.empty:
        return _EMPTY_STATS.copy()
    eff = res.pivot_table(index="module_id", columns="age_label", values="effect")
    direction = eff.get("p75", 0.0) - eff.get("p25", 0.0)
    stat = res.drop_duplicates("module_id").set_index("module_id")[["pvalue_ftest", "fdr_ftest"]]
    out = stat.join(direction.rename("direction")).reset_index()
    return out.rename(columns={"module_id": "unit", "pvalue_ftest": "p", "fdr_ftest": "fdr"})


# --------------------------------------------------------------------------- #
# Analysis A: gene-level orthogonal corroboration
# --------------------------------------------------------------------------- #
def gene_corroboration(
    cohort: str, region: str, variant: str, trait: str, sample_table: pd.DataFrame,
    alpha: float, switch_set: str = "module",
) -> tuple[pd.DataFrame, dict]:
    covs, age_col = _covariates(cohort, trait), _age_col(cohort)
    sample_ids = list(sample_table["sample_id"].astype(str))
    sw_meta, sw_wide, _ = _load_switch_scores(cohort, region, variant, trait, sample_ids)

    if switch_set == "module":
        # switch_pos = member of a trait-associated co-switch module (manuscript unit;
        # powered where marginal per-gene calls are not). No per-gene switch spline needed
        # -- Analysis A is a direction-free Fisher on switch_pos x psi_pos.
        mod_genes = _module_switch_genes(cohort, region, variant, trait, alpha)
        sw_stat = sw_meta[["gene", "n_transcripts"]].drop_duplicates("gene").copy()
        sw_stat["switch_pos"] = sw_stat["gene"].isin(mod_genes)
        print(f"[{cohort}/{region}/{trait}] module switch-set: {int(sw_stat['switch_pos'].sum())} "
              f"member genes in trait-associated modules (of {len(sw_stat)} scored)", flush=True)
    else:
        print(f"[{cohort}/{region}/{trait}] marginal switch-score {trait} test over "
              f"{sw_wide.shape[1] - 1} genes ...", flush=True)
        sw_stat = _assoc_stats(sw_wide, sample_table, trait, covs, age_col)
        sw_stat = sw_stat.merge(sw_meta.rename(columns={"gene_id": "unit"}), on="unit", how="left")
        sw_stat["switch_pos"] = sw_stat["fdr"] <= alpha
        sw_stat = sw_stat[["gene", "n_transcripts", "switch_pos"]]

    universe_genes = set(sw_stat["gene"].dropna())
    psi_meta, psi_wide = _load_splice(cohort, region, sample_ids, genes=universe_genes)
    psi_wide = _prune_units(psi_wide)
    psi_meta = psi_meta[psi_meta["event_id"].isin(psi_wide.columns)]
    print(f"[{cohort}/{region}/{trait}] splice {trait} test over {psi_wide.shape[1] - 1} "
          f"testable events ({psi_meta['gene'].nunique()} genes) ...", flush=True)
    psi_stat = _assoc_stats(psi_wide, sample_table, trait, covs, age_col)
    psi_stat = psi_stat.merge(psi_meta[["event_id", "gene"]].rename(columns={"event_id": "unit"}),
                              on="unit", how="left")

    # Per gene: min event p, Bonferroni across the gene's events, then BH across genes.
    g = psi_stat.dropna(subset=["p"]).groupby("gene")
    per_gene = pd.DataFrame({
        "n_psi_events": g.size(),
        "psi_min_p": g["p"].min(),
    }).reset_index()
    per_gene["psi_gene_p"] = np.minimum(per_gene["psi_min_p"] * per_gene["n_psi_events"], 1.0)
    per_gene["psi_gene_fdr"] = stats.false_discovery_control(per_gene["psi_gene_p"], method="bh")
    per_gene["psi_pos"] = per_gene["psi_gene_fdr"] <= alpha

    genes = sw_stat.merge(per_gene, on="gene", how="inner")  # genes with switch score AND >=1 PSI event
    genes = genes.dropna(subset=["psi_gene_fdr"])

    # 2x2 Fisher (headline).
    a = int(((genes["switch_pos"]) & (genes["psi_pos"])).sum())
    b = int(((genes["switch_pos"]) & (~genes["psi_pos"])).sum())
    c = int(((~genes["switch_pos"]) & (genes["psi_pos"])).sum())
    d = int(((~genes["switch_pos"]) & (~genes["psi_pos"])).sum())
    fisher_or, fisher_p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

    # Transcript-count-adjusted logistic OR (identifiability confounder).
    adj_or = adj_p = float("nan")
    try:
        import statsmodels.api as sm
        m = genes.dropna(subset=["n_transcripts"]).copy()
        X = pd.DataFrame({
            "const": 1.0,
            "switch_pos": m["switch_pos"].astype(float),
            "log_n_tx": np.log1p(m["n_transcripts"].astype(float)),
        })
        fit = sm.Logit(m["psi_pos"].astype(float), X).fit(disp=0)
        adj_or = float(np.exp(fit.params["switch_pos"]))
        adj_p = float(fit.pvalues["switch_pos"])
    except Exception as exc:  # pragma: no cover - defensive
        print(f"[{region}] adjusted logistic skipped: {exc}", flush=True)

    rate_pos = a / max(a + b, 1)
    rate_neg = c / max(c + d, 1)
    summary = {
        "switch_set": switch_set,
        "n_genes": int(len(genes)),
        "n_switch_pos": int(genes["switch_pos"].sum()),
        "n_psi_pos": int(genes["psi_pos"].sum()),
        "contingency": {"switch_psi": a, "switch_only": b, "psi_only": c, "neither": d},
        "psi_pos_rate_in_switch_genes": rate_pos,
        "psi_pos_rate_in_background": rate_neg,
        "fisher_or": float(fisher_or),
        "fisher_p": float(fisher_p),
        "adjusted_or": adj_or,
        "adjusted_p": adj_p,
        "alpha": alpha,
    }
    return genes, summary


# --------------------------------------------------------------------------- #
# Analysis B: event-resolved concordance for the manuscript switches
# --------------------------------------------------------------------------- #
def _parse_junction(junction: str) -> tuple[str, int, int] | None:
    m = _JUNC_RE.search(str(junction))
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def _match_event(chrom: str, s: int, e: int, gene_events: pd.DataFrame) -> tuple[str | None, int]:
    """Best PSI event for a switch junction; returns (event_id, endpoints_matched 0/1/2)."""
    best, best_hits = None, 0
    for _, ev in gene_events.iterrows():
        if str(ev["chrom"]) != chrom:
            continue
        cs = ev["coords"]
        hits = sum(any(abs(x - t) <= _COORD_MATCH_TOL for x in cs) for t in (s, e))
        if hits > best_hits:
            best, best_hits = ev["event_id"], hits
    return best, best_hits


def event_resolved(
    region: str, variant: str, trait: str, sample_table: pd.DataFrame, alpha: float
) -> tuple[pd.DataFrame, dict]:
    ev = pd.read_parquet(rel("real_data", "_m", "deep_dive", "deep_dive_events.parquet"))
    ev = ev[ev["junction_in_switch_pair"] == True].copy()  # noqa: E712
    ev["gene"] = _strip_ver(ev["ens"])

    covs, age_col = _covariates("brainseq", trait), _age_col("brainseq")
    sample_ids = list(sample_table["sample_id"].astype(str))
    genes = set(ev["gene"])
    psi_meta, psi_wide = _load_splice("brainseq", region, sample_ids, genes=genes)
    sw_meta, sw_wide, _ = _load_switch_scores("brainseq", region, variant, trait, sample_ids)
    sw_meta["gene"] = _strip_ver(sw_meta["gene_id"])
    sw_wide_sub = sw_wide[["sample_id"] + [g for g in sw_meta.loc[sw_meta["gene"].isin(genes), "gene_id"]
                                          if g in sw_wide.columns]]

    # Match each resolved junction to a PSI event.
    rows = []
    for _, r in ev.iterrows():
        parsed = _parse_junction(r["junction"])
        gene_events = psi_meta[psi_meta["gene"] == r["gene"]]
        eid, hits = (None, 0) if parsed is None else _match_event(*parsed, gene_events)
        rows.append({
            "gene_name": r["gene_name"], "gene": r["gene"], "trait": r.get("trait"),
            "junction": r["junction"], "matched_event": eid, "endpoints_matched": hits,
        })
    matched = pd.DataFrame(rows)

    # Spline test the matched PSI events + the corresponding switch scores.
    want_events = [e for e in matched["matched_event"].dropna().unique()
                   if e in psi_wide.columns]
    psi_sub = _prune_units(psi_wide[["sample_id"] + want_events]) if want_events else psi_wide[["sample_id"]]
    psi_stat = (_assoc_stats(psi_sub, sample_table, trait, covs, age_col).set_index("unit")
                if psi_sub.shape[1] > 1 else pd.DataFrame())
    sw_stat = (_assoc_stats(_prune_units(sw_wide_sub), sample_table, trait, covs, age_col).set_index("unit")
               if sw_wide_sub.shape[1] > 1 else pd.DataFrame())
    gene_to_uid = dict(zip(sw_meta["gene"], sw_meta["gene_id"]))

    def _lookup(tbl, key):
        if key is not None and not tbl.empty and key in tbl.index:
            return float(tbl.loc[key, "p"]), float(tbl.loc[key, "direction"])
        return float("nan"), float("nan")

    out = []
    for _, r in matched.iterrows():
        psi_p, psi_dir = _lookup(psi_stat, r["matched_event"])
        sw_p, sw_dir = _lookup(sw_stat, gene_to_uid.get(r["gene"]))
        both_sig = np.isfinite(psi_p) and np.isfinite(sw_p) and psi_p <= alpha and sw_p <= alpha
        concordant = bool(both_sig and (np.sign(psi_dir) == np.sign(sw_dir)))
        out.append({
            **r.to_dict(), "tested_trait": trait,
            "psi_p": psi_p, "psi_direction": psi_dir,
            "switch_p": sw_p, "switch_direction": sw_dir,
            "both_significant": bool(both_sig), "concordant": concordant,
        })
    tab = pd.DataFrame(out)
    summary = {
        "n_resolved_junctions": int(len(tab)),
        "n_matched_to_psi": int(tab["matched_event"].notna().sum()),
        "n_both_significant": int(tab["both_significant"].sum()),
        "n_concordant": int(tab["concordant"].sum()),
        "alpha": alpha,
    }
    return tab, summary


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def run_region(cohort: str, region: str, variant: str, trait: str, alpha: float,
               analyses: list[str], switch_set: str = "module") -> dict:
    bundle = load_dataset_bundle(_bundle_path(cohort, region, trait))
    st = bundle.sample_table.copy()
    st["sample_id"] = st["sample_id"].astype(str)
    out = _out_dir(cohort, region, trait)
    summary: dict = {"cohort": cohort, "region": region, "variant": variant,
                     "trait": trait, "switch_set": switch_set}

    if "gene" in analyses:
        genes, s = gene_corroboration(cohort, region, variant, trait, st, alpha,
                                      switch_set=switch_set)
        genes.to_parquet(out / f"gene_corroboration_{switch_set}.parquet",
                         index=False, compression="zstd")
        summary["gene_corroboration"] = s
        print(f"[{cohort}/{region}/{trait}] GENE: switch-genes splice-{trait}-sig "
              f"{s['psi_pos_rate_in_switch_genes']:.1%} vs background "
              f"{s['psi_pos_rate_in_background']:.1%} | "
              f"Fisher OR={s['fisher_or']:.2f} p={s['fisher_p']:.2e} | "
              f"adj OR={s['adjusted_or']:.2f} p={s['adjusted_p']:.2e}", flush=True)
    if "event" in analyses:
        if cohort == "gtex":
            print(f"[{cohort}/{region}/{trait}] EVENT skipped: event-resolved analysis is "
                  "BrainSEQ-PSI-specific (deep-dive coloc switches).", flush=True)
        else:
            tab, s = event_resolved(region, variant, trait, st, alpha)
            tab.to_parquet(out / "event_resolved.parquet", index=False, compression="zstd")
            summary["event_resolved"] = s
            print(f"[{cohort}/{region}/{trait}] EVENT: {s['n_matched_to_psi']}/"
                  f"{s['n_resolved_junctions']} junctions matched a PSI event; "
                  f"{s['n_both_significant']} both-significant, {s['n_concordant']} concordant",
                  flush=True)

    (out / f"summary_{switch_set}.json").write_text(json.dumps(summary, indent=2))
    print(f"[{region}] written to {out}", flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cohort", choices=["brainseq", "gtex"], default="brainseq",
                        help="brainseq (LIBD PSI) or gtex (STAR within-gene junction usage; "
                             "requires build_gtex_junction_usage to have produced "
                             "inputs/processed/gtex_v11/<region>/junction_usage.parquet).")
    parser.add_argument("--region", action="append", dest="regions",
                        help="Region(s); repeatable. Default: all BrainSEQ aging regions "
                             "(brainseq age), caudate (brainseq dx), or all 13 GTEx regions.")
    parser.add_argument("--trait", choices=["age", "dx"], default="age",
                        help="age -> non-linear age (brainseq_v1 controls / GTEx); "
                             "dx -> SCZD-vs-Control in the caudate SCZD cohort (brainseq only). "
                             "The deep-dive coloc switches are disease-anchored, so dx is the "
                             "fair test for the event-resolved analysis.")
    parser.add_argument("--switch-set", choices=["module", "marginal"], default="module",
                        dest="switch_set",
                        help="Analysis A switch-gene set. module (default) = members of "
                             "trait-associated co-switch modules (manuscript unit, powered in "
                             "Phase 2). marginal = per-gene switch feature-score FDR (conservative; "
                             "yields 0 switch genes in the small regions).")
    parser.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    parser.add_argument("--alpha", type=float, default=0.05, help="FDR threshold (default 0.05).")
    parser.add_argument("--analysis", choices=["gene", "event", "both"], default="both")
    args = parser.parse_args()

    analyses = ["gene", "event"] if args.analysis == "both" else [args.analysis]

    if args.cohort == "gtex":
        if args.trait == "dx":
            raise SystemExit("--trait dx is BrainSEQ-only (GTEx has no case/control cohort).")
        regions = args.regions or list(GTEX_REGIONS)
        missing = [r for r in regions
                   if not _splice_path("gtex", r).exists()]
        if missing:
            raise SystemExit(
                "Missing GTEx junction_usage.parquet for: " + ", ".join(missing) + ". Run "
                "`python -m isograph_benchmark.inputs.build_gtex_junction_usage --region <r>` "
                "(SLURM: real_data/gtex/_h/... ) first.")
    elif args.trait == "dx":
        regions = args.regions or ["caudate"]
        if regions != ["caudate"]:
            raise SystemExit("--trait dx is caudate-only (the SCZD cohort is caudate).")
    else:
        regions = args.regions or list(BRAINSEQ_AGING_REGIONS)

    for region in regions:
        run_region(args.cohort, region, args.variant, args.trait, args.alpha, analyses,
                   switch_set=args.switch_set)


if __name__ == "__main__":
    main()
