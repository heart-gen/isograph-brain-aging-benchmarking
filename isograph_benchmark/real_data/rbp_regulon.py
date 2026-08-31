"""Stage 2: are IsoGraph co-switch modules RBP regulons?

The coding-consequence analysis showed the switch axis is a UTR-remodeling layer (UTR change
enriched 1.27x, 10/10 regions). 3'UTR usage is heavily RBP-controlled, which motivates the
coordination hypothesis: a co-switch module is coordinated because its member genes share a
common trans-acting RBP — the one thing cis-sQTL anchoring does not explain.

Using the per-(transcript, RBP) motif counts from `rbp_scan.py`, for each switch gene we call
whether the switch GAINS or LOSES an RBP binding site (motif present in one switch-pair isoform,
absent in the other). Then, per module, we test whether a given RBP's site-switching is
over-represented among the module's genes relative to the pooled switch-gene background
(hypergeometric) — a significant RBP marks the module as a candidate regulon for it. The
within-pair "gained/lost" call compares two isoforms of the SAME gene, controlling transcript
length/composition. Results are stratified GO-invisible vs GO-visible.

Output (07_rbp_regulation/_m/rbp/): rbp_switch_calls.parquet (per gene x RBP), rbp_regulon.parquet
(per module x RBP enrichment), RBP_REGULON.md, and a cross-region meta.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from isograph_benchmark.paths import region_store, stage_out
from isograph_benchmark.real_data.coloc_prep import load_switch_genes
from isograph_benchmark.real_data.qtl_anchoring import _bare, build_gene_sets

_RBP_DIR = stage_out("regulation", "rbp")
_COUNTS = _RBP_DIR / "rbp_counts.parquet"
_FLANK_NOTE = 100        # intronic flank window (nt); mirrors rbp_scan_intronic._FLANK
# per-scope Stage-1 count tables; "combined" unions the mature + intronic presence
_SCOPE_COUNTS = {
    "mature": [_RBP_DIR / "rbp_counts.parquet"],
    "intronic": [_RBP_DIR / "rbp_counts_intronic.parquet"],
    "combined": [_RBP_DIR / "rbp_counts.parquet", _RBP_DIR / "rbp_counts_intronic.parquet"],
}
_TREE_OF = {**{r: "brainseq" for r in ("caudate_sczd", "caudate", "hippocampus", "dlpfc")}}
_REGIONS = [
    ("brainseq", "caudate_sczd"), ("brainseq", "caudate"), ("brainseq", "hippocampus"),
    ("brainseq", "dlpfc"),
    *[("gtex", r) for r in (
        "amygdala", "anterior_cingulate_cortex_ba24", "caudate_basal_ganglia",
        "cerebellar_hemisphere", "cerebellum", "cortex", "frontal_cortex_ba9",
        "hippocampus", "hypothalamus", "nucleus_accumbens_basal_ganglia",
        "putamen_basal_ganglia", "spinal_cord_cervical_c_1", "substantia_nigra")],
]


def _presence(counts: pd.DataFrame, unit: str = "rbp") -> dict[tuple[str, str], int]:
    """(transcript_id, unit) -> hit count (missing = 0)."""
    return {(t, r): int(c) for t, r, c in
            counts[["transcript_id", unit, "count"]].itertuples(index=False)}


def _load_counts(paths: list[Path], unit: str, bg_mode: str,
                 region_class: str | None) -> pd.DataFrame:
    """Stage-1 counts, filtered to one background mode and optionally one region class.

    ``rbp_scan`` now emits a row per (transcript, unit, region, bg_mode). Collapsing over
    ``region`` reproduces the old whole-mature-transcript view; selecting one region class
    asks the question within 5'UTR / CDS / 3'UTR, so a switch that merely lengthens the
    3'UTR cannot masquerade as regulatory rewiring.
    """
    frames = []
    for p in paths:
        df = pd.read_parquet(p)
        if "bg_mode" in df.columns:
            df = df[df["bg_mode"] == bg_mode]
        if region_class is not None:
            if "region" not in df.columns:
                raise SystemExit(f"{p} has no region column; re-run rbp_scan.py for the "
                                 f"5'UTR/CDS/3'UTR partition.")
            df = df[df["region"] == region_class]
        frames.append(df)
    counts = pd.concat(frames, ignore_index=True)
    return counts.groupby(["transcript_id", unit], as_index=False)["count"].sum()


def _gene_covariates(pool_genes: pd.DataFrame, sp: pd.DataFrame,
                     opportunity: pd.DataFrame | None) -> pd.DataFrame:
    """Per-gene opportunity covariates for the adjusted model.

    A longer, more AU-rich, more UTR-heavy switch interval offers more chances for any motif
    to be gained or lost, so an unadjusted hypergeometric confuses "this module is an RBP
    regulon" with "this module's genes have long UTRs". These are the covariates that make
    the comparison length-, composition- and region-matched without having to construct
    matched null sets.
    """
    n_tx = sp.groupby("gene")[["transcript_id_1", "transcript_id_2"]].apply(
        lambda g: len(set(g["transcript_id_1"]) | set(g["transcript_id_2"])))
    cov = pd.DataFrame({"gene": n_tx.index, "n_transcripts": n_tx.to_numpy()})
    cov["log_n_transcripts"] = np.log1p(cov["n_transcripts"])

    if opportunity is None or opportunity.empty:
        for c in ("log_length", "gc", "frac_5utr", "frac_cds", "frac_3utr"):
            cov[c] = 0.0
        return cov

    tx_long = pd.concat([
        sp[["gene", "transcript_id_1"]].rename(columns={"transcript_id_1": "transcript_id"}),
        sp[["gene", "transcript_id_2"]].rename(columns={"transcript_id_2": "transcript_id"}),
    ], ignore_index=True).drop_duplicates()
    opp = tx_long.merge(opportunity, on="transcript_id", how="inner")
    if opp.empty:
        for c in ("log_length", "gc", "frac_5utr", "frac_cds", "frac_3utr"):
            cov[c] = 0.0
        return cov

    total = opp.groupby("gene")["length"].sum()
    gc = (opp.assign(w=opp["gc"] * opp["length"]).groupby("gene")["w"].sum()
          / total.replace(0, np.nan))
    frac = (opp.pivot_table(index="gene", columns="region", values="length",
                            aggfunc="sum", fill_value=0)
            .div(total.replace(0, np.nan), axis=0))
    prof = pd.DataFrame({"gene": total.index, "log_length": np.log1p(total.to_numpy()),
                         "gc": gc.reindex(total.index).to_numpy()})
    for cls, col in (("5utr", "frac_5utr"), ("cds", "frac_cds"), ("3utr", "frac_3utr")):
        prof[col] = (frac[cls].reindex(total.index).to_numpy()
                     if cls in frac.columns else 0.0)
    out = cov.merge(prof, on="gene", how="left")
    for c in ("log_length", "gc", "frac_5utr", "frac_cds", "frac_3utr"):
        out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0.0)
    return out


_GLM_COVARIATES = ("log_length", "gc", "frac_5utr", "frac_cds", "frac_3utr",
                   "log_n_transcripts")
# Estimability guards for the adjusted GLM (see _regulon_glm).  The observed fits are
# sharply bimodal -- every estimable cell lands at SE < 10 and every separated one at
# SE > 1e4, with the band between them empty -- so _MAX_SE sits in that gap and the exact
# value is immaterial.  _MAX_ABS_COEF = 30 is a log-odds ratio of ~1e13, i.e. far outside
# anything a real effect can produce; it only catches numeric blow-up.
_SEP_TOL = 1e-8          # fitted probability treated as pinned at 0 or 1
_MAX_SE = 100.0
_MAX_ABS_COEF = 30.0


def _gene_tags(art, fdr: float) -> tuple[pd.DataFrame, str]:
    """(gene, module_id, go_invisible) tags + provenance.

    Primary pool = phenotype-associated switch genes (`load_switch_genes`). When that is empty
    (regions with no FDR-significant phenotype/age association, e.g. several GTEx tissues), fall
    back to the region's FULL module-gene pool so the RBP-regulon question — "is this co-switch
    module coordinated by an RBP?" — is still asked; the GO-invisible tag is preserved from
    `build_gene_sets`, and a `pool_source` flag records which universe was used.
    """
    sg = load_switch_genes(art, fdr)
    if not sg.empty:
        return sg[["gene", "module_id", "go_invisible"]].drop_duplicates(), "switch_genes"
    enrich_path = art.parent / "module_enrichment" / "isograph_modules.parquet"
    # _gene_tags is called with the artifact dir only, so name the region from the path
    # (…/<cohort>/<region>/_m/isograph_vae) rather than from unavailable arguments.
    inv = build_gene_sets(art, enrich_path, fdr,
                          context=f"rbp_regulon {art.parent.parent.name}"
                          ).get("go_invisible_modules", set())
    mods = pd.read_parquet(art / "modules.parquet")
    mods["gene"] = _bare(mods["gene_id"])
    mods["module_id"] = mods["module_id"].astype(str)
    mods["go_invisible"] = mods["gene"].isin(inv)
    return mods[["gene", "module_id", "go_invisible"]].drop_duplicates(), "module_genes"


def _switch_calls(region_tree: str, region: str, pres: dict, rbps: list[str],
                  fdr: float) -> pd.DataFrame:
    """Per (gene, RBP): does the switch gain/lose the motif (present in one isoform only)?"""
    art = region_store(region_tree, region, "isograph_vae")
    sp_path = art / "module_interpret" / "structure_switch_pairs.parquet"
    if not sp_path.exists():
        return pd.DataFrame()
    tag, pool_source = _gene_tags(art, fdr)
    if tag.empty:
        return pd.DataFrame()
    sp = pd.read_parquet(sp_path)
    sp["gene"] = _bare(sp["gene_id"])
    tag = tag.groupby("gene").agg(go_invisible=("go_invisible", "max"),
                                  module_id=("module_id", "first")).reset_index()
    sp = sp.merge(tag, on="gene", how="inner")
    if sp.empty:
        return pd.DataFrame()

    rows = []
    for gene, sub in sp.groupby("gene"):
        go_inv = bool(sub["go_invisible"].iloc[0])
        mod = sub["module_id"].iloc[0]
        for rbp in rbps:
            # gained/lost across ANY of the gene's switch pairs
            switched = False
            for r in sub.itertuples():
                c1 = pres.get((r.transcript_id_1, rbp), 0)
                c2 = pres.get((r.transcript_id_2, rbp), 0)
                if (c1 > 0) != (c2 > 0):
                    switched = True
                    break
            rows.append((region, gene, mod, go_inv, rbp, switched, pool_source))
    return pd.DataFrame(rows, columns=["region", "gene", "module_id", "go_invisible",
                                       "rbp", "switched", "pool_source"])


def _regulon_enrich(calls: pd.DataFrame) -> pd.DataFrame:
    """Per (region, module, RBP): hypergeometric over-representation of RBP site-switching
    among the module's genes vs the region's switch-gene pool."""
    out = []
    for region, rc in calls.groupby("region"):
        genes = rc["gene"].unique()
        N = len(genes)
        pool_source = rc["pool_source"].iloc[0] if "pool_source" in rc.columns else "switch_genes"
        # pool: number of genes with this RBP switched (region-wide)
        pool = rc.groupby("rbp")["switched"].sum().to_dict()
        gene_inv = rc.drop_duplicates("gene").set_index("gene")["go_invisible"]
        for module, mg in rc.groupby("module_id"):
            module_genes = mg["gene"].unique()
            n = len(module_genes)
            if n < 3:
                continue
            go_inv = bool(gene_inv.reindex(module_genes).mode().iloc[0]) \
                if len(module_genes) else False
            per_rbp = mg.groupby("rbp")["switched"].sum()
            for rbp, k in per_rbp.items():
                K = int(pool.get(rbp, 0))
                if K == 0 or k == 0:
                    continue
                p = float(hypergeom.sf(int(k) - 1, N, K, n))
                exp = n * K / N
                out.append({
                    "region": region, "module_id": module, "go_invisible": go_inv,
                    "rbp": rbp, "module_size": n, "n_switched": int(k),
                    "pool_switched": K, "pool_size": N,
                    "expected": exp, "enrichment": (k / exp if exp > 0 else np.nan),
                    "p": p, "pool_source": pool_source})
    res = pd.DataFrame(out)
    if not res.empty:
        res["q"] = multipletests(res["p"], method="fdr_bh")[1]
    return res


def _safe_exp(x: float) -> float:
    return float(np.exp(np.clip(x, -700, 700)))


def _regulon_glm(calls: pd.DataFrame, covariates: pd.DataFrame,
                 unit: str, min_module_genes: int = 3) -> pd.DataFrame:
    """Covariate-adjusted binomial GLM per (region, module, unit).

    ``switched ~ in_module + log(interval length) + GC + region composition +
    log(n_transcripts)``, fit over the region's whole gene pool.  This is the tractable form
    of the reviewer's "background matched for length, GC, region and UTR-vs-CDS context":
    adjust in the model rather than construct matched null sets, which the repo has no
    utility for and which would be far more expensive.  The hypergeometric p/q are retained
    by the caller as legacy columns.

    **Estimability gate.**  Only cells with an existing MLE for the ``in_module`` coefficient
    enter the BH family.  Two layers:

    1. *Complete separation*, detected from the (in_module x switched) 2x2 before fitting.
       Any empty cell makes ``in_module`` a perfect predictor and drives its coefficient to
       +/-inf.  The dominant case is ``module_switched == 0`` -- exactly the cell the
       hypergeometric arm already declines to test (``_regulon_enrich`` skips ``k == 0``), so
       this also keeps the two multiplicity families over the same set of module x unit cells.
    2. *Quasi-complete separation* and numeric degeneracy, detected after fitting from
       non-convergence, fitted probabilities pinned at 0/1, or an implausible estimate.
       These arise when a module contributes a single switched gene that the covariates can
       isolate; IRLS converges to a nominal answer whose Wald statistic is meaningless.

    Excluded cells keep ``model_status`` and their 2x2 counts for audit but carry no p-value,
    so ``multipletests`` never sees them.  Without the gate they inflate the family by ~7%
    and a handful surface as q = 0 with odds ratios of ~1e304.
    """
    import warnings

    import statsmodels.api as sm
    from statsmodels.tools.sm_exceptions import (ConvergenceWarning,
                                                 PerfectSeparationError,
                                                 PerfectSeparationWarning)

    cov_cols = [c for c in _GLM_COVARIATES if c in covariates.columns]
    rows: list[dict] = []
    for region, rc in calls.groupby("region"):
        rc = rc.merge(covariates, on="gene", how="left")
        for c in cov_cols:
            rc[c] = pd.to_numeric(rc[c], errors="coerce").fillna(0.0)
        gene_cov = rc.drop_duplicates("gene").set_index("gene")
        for unit_value, uc in rc.groupby(unit):
            genes = uc.drop_duplicates("gene")
            outcome = genes["switched"].astype(int)
            module_of = genes["module_id"].astype(str)
            design_cov = gene_cov.loc[genes["gene"], cov_cols].to_numpy(float)
            y = outcome.to_numpy()
            n_universe, n_universe_switched = int(len(genes)), int(y.sum())
            for module in module_of.unique():
                in_module = (module_of == module).astype(float).to_numpy()
                n_mod = int(in_module.sum())
                n_mod_switched = int(y[in_module > 0].sum())
                base = {"region": region, "module_id": module, unit: unit_value,
                        "module_genes": n_mod, "universe_genes": n_universe,
                        "module_switched": n_mod_switched,
                        "universe_switched": n_universe_switched}
                if n_mod < min_module_genes:
                    rows.append({**base, "model_status": "module_too_small"})
                    continue
                if outcome.nunique() < 2 or len(np.unique(in_module)) < 2:
                    rows.append({**base, "model_status": "outcome_or_predictor_constant"})
                    continue
                # (in_module x switched) 2x2; an empty cell => complete separation.
                cells = (n_mod_switched, n_mod - n_mod_switched,
                         n_universe_switched - n_mod_switched,
                         (n_universe - n_mod) - (n_universe_switched - n_mod_switched))
                if min(cells) == 0:
                    rows.append({**base, "model_status": "separated_zero_cell"})
                    continue
                X = np.column_stack([np.ones(len(genes)), in_module, design_cov])
                # Drop covariate columns that are constant in this region; they add nothing
                # and make the design singular.
                keep = np.ones(X.shape[1], dtype=bool)
                keep[2:] = X[:, 2:].std(axis=0) > 1e-12
                try:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter("always")
                        fit = sm.GLM(y, X[:, keep],
                                     family=sm.families.Binomial()).fit(maxiter=200)
                    # statsmodels downgraded perfect separation from an error to a warning;
                    # blanket-suppressing warnings here is what let these through before.
                    if any(issubclass(w.category,
                                      (PerfectSeparationWarning, ConvergenceWarning))
                           for w in caught):
                        rows.append({**base, "model_status": "quasi_separated"})
                        continue
                    coef = float(fit.params[1])
                    se = float(fit.bse[1])
                    p = float(fit.pvalues[1])
                    lo, hi = (float(v) for v in fit.conf_int()[1])
                    if not all(np.isfinite([coef, se, p, lo, hi])):
                        raise ValueError("non-finite logistic estimate")
                    mu = np.asarray(fit.fittedvalues, dtype=float)
                    degenerate = (not getattr(fit, "converged", True)
                                  or mu.max() > 1.0 - _SEP_TOL or mu.min() < _SEP_TOL
                                  or se > _MAX_SE or abs(coef) > _MAX_ABS_COEF)
                    if degenerate:
                        rows.append({**base, "coefficient": coef, "standard_error": se,
                                     "model_status": "quasi_separated"})
                        continue
                    rows.append({**base, "coefficient": coef, "standard_error": se,
                                 "odds_ratio": _safe_exp(coef),
                                 "ci_low": _safe_exp(lo), "ci_high": _safe_exp(hi),
                                 "pvalue_glm": p, "model_status": "fit"})
                except (PerfectSeparationError, ValueError, np.linalg.LinAlgError) as err:
                    rows.append({**base, "model_status": f"fit_failed:{type(err).__name__}"})

    res = pd.DataFrame(rows)
    if res.empty:
        return res
    # Stable schema: a scope in which nothing is estimable must still carry the estimate
    # columns, so callers never have to test for their existence.
    for col in ("coefficient", "standard_error", "odds_ratio", "ci_low", "ci_high",
                "pvalue_glm", "qvalue_glm"):
        if col not in res.columns:
            res[col] = np.nan
    # Only estimable cells (model_status == "fit") carry a p-value, so the BH family is
    # exactly the set of tests that were actually performed.
    ok = res["pvalue_glm"].notna()
    if ok.any():   # multipletests divides by the test count and cannot take an empty family
        res.loc[ok, "qvalue_glm"] = multipletests(
            res.loc[ok, "pvalue_glm"], method="fdr_bh")[1]
    return res


_UNIT_COUNTS = {
    "rbp": {"mature": [_RBP_DIR / "rbp_counts.parquet"],
            "intronic": [_RBP_DIR / "rbp_counts_intronic.parquet"],
            "combined": [_RBP_DIR / "rbp_counts.parquet",
                         _RBP_DIR / "rbp_counts_intronic.parquet"]},
    "family_id": {"mature": [_RBP_DIR / "rbp_family_counts.parquet"],
                  "intronic": [_RBP_DIR / "rbp_family_counts_intronic.parquet"],
                  "combined": [_RBP_DIR / "rbp_family_counts.parquet",
                               _RBP_DIR / "rbp_family_counts_intronic.parquet"]},
}


def run(fdr: float, scope: str = "mature", unit: str = "rbp", bg_mode: str = "composition",
        region_class: str | None = None) -> None:
    paths = _UNIT_COUNTS[unit][scope]
    missing = [p for p in paths if not p.exists()]
    if missing:
        hint = "rbp_scan.py" if scope == "mature" else "rbp_scan_intronic.py"
        raise SystemExit(f"{missing[0]} missing; run {hint} (motif env) first "
                         f"(scope={scope} needs {[p.name for p in paths]}).")
    counts = _load_counts(paths, unit, bg_mode, region_class)
    pres = _presence(counts, unit)
    rbps = sorted(counts[unit].unique())
    opp_path = _RBP_DIR / "rbp_scan_opportunity.parquet"
    opportunity = pd.read_parquet(opp_path) if opp_path.exists() else None
    if opportunity is None:
        print("  NOTE: rbp_scan_opportunity.parquet missing; the adjusted GLM will run "
              "without the length/GC/region covariates. Re-run rbp_scan.py to produce it.")
    print(f"[{scope}/{unit}/bg={bg_mode}"
          f"{'/region=' + region_class if region_class else ''}] "
          f"loaded motif counts: {len(rbps)} {unit}s, "
          f"{counts['transcript_id'].nunique():,} transcripts")

    calls = []
    for tree, region in _REGIONS:
        c = _switch_calls(tree, region, pres, rbps, fdr)
        if not c.empty:
            calls.append(c)
            print(f"  {region}: {c['gene'].nunique()} switch genes")
    calls = pd.concat(calls, ignore_index=True)
    calls = calls.rename(columns={"rbp": unit}) if unit != "rbp" else calls
    _RBP_DIR.mkdir(parents=True, exist_ok=True)
    suffix = "" if scope == "mature" else f"_{scope}"
    if unit != "rbp":
        suffix += "_family"
    if region_class:
        suffix += f"_{region_class}"
    calls.to_parquet(_RBP_DIR / f"rbp_switch_calls{suffix}.parquet", index=False)

    # Per-gene opportunity covariates, pooled across regions (a gene's switch interval is a
    # property of its transcripts, not of the region it was called in).
    sp_frames = []
    for tree, region in _REGIONS:
        f = region_store(tree, region, "isograph_vae", "module_interpret",
                "structure_switch_pairs.parquet")
        if f.exists():
            d = pd.read_parquet(f, columns=["gene_id", "transcript_id_1", "transcript_id_2"])
            d["gene"] = _bare(d["gene_id"])
            sp_frames.append(d)
    sp_all = pd.concat(sp_frames, ignore_index=True).drop_duplicates()
    covariates = _gene_covariates(calls, sp_all, opportunity)

    reg = _regulon_enrich(calls.rename(columns={unit: "rbp"}) if unit != "rbp" else calls)
    if unit != "rbp":
        reg = reg.rename(columns={"rbp": unit})
    glm = _regulon_glm(calls, covariates, unit)
    if not glm.empty:
        # Keep the GLM's 2x2 counts: for module x unit cells the hypergeometric arm skips
        # (k == 0) they are the only record of why the cell was not testable.
        reg = reg.merge(glm, on=["region", "module_id", unit], how="outer")
    reg.to_parquet(_RBP_DIR / f"rbp_regulon{suffix}.parquet", index=False)
    _write_report(reg, calls, _RBP_DIR / f"RBP_REGULON{suffix}.md", scope, unit)
    n_sig = int((reg["q"] < 0.05).sum()) if not reg.empty else 0
    n_glm = (int((reg["qvalue_glm"] < 0.05).sum())
             if not reg.empty and "qvalue_glm" in reg.columns else 0)
    n_family = (int(reg["pvalue_glm"].notna().sum())
                if not reg.empty and "pvalue_glm" in reg.columns else 0)
    print(f"[{scope}/{unit}] candidate regulons: hypergeometric q<0.05 {n_sig}; "
          f"covariate-adjusted GLM q<0.05 {n_glm} of {n_family} estimable tests across "
          f"{reg['region'].nunique() if not reg.empty else 0} regions")
    if not reg.empty and "model_status" in reg.columns:
        print("  GLM cell status: " + ", ".join(
            f"{k} {v}" for k, v in reg["model_status"].value_counts().items()))


def _write_report(reg: pd.DataFrame, calls: pd.DataFrame,
                  out_path: Path | None = None, scope: str = "mature",
                  unit: str = "rbp") -> None:
    out_path = out_path or (_RBP_DIR / "RBP_REGULON.md")
    scope_note = {
        "mature": "Motifs are scanned on the **mature transcript** (exonic + UTR) sequence.",
        "intronic": "Motifs are scanned on **intronic splice-site flanks** "
                    f"(pre-mRNA sense; up to {_FLANK_NOTE} nt into each intron), the binding "
                    "niche for splicing-regulatory RBPs invisible to the mature-transcript scan.",
        "combined": "Motif presence unions the **mature-transcript** and **intronic "
                    "splice-site flank** scans.",
    }[scope]
    unit_label = "RBP" if unit == "rbp" else "motif family"
    sig = reg[reg["q"] < 0.05].sort_values("q") if not reg.empty else reg
    n_sig = len(sig)
    n_goinv = int(sig["go_invisible"].sum()) if n_sig else 0

    def _num(frame, col):
        return (pd.to_numeric(frame[col], errors="coerce") if col in frame.columns
                else pd.Series(np.nan, index=frame.index))

    has_glm = (not reg.empty) and "qvalue_glm" in reg.columns
    lines = [
        f"# Candidate RBP regulons among IsoGraph co-switch modules ({scope} scope)", "",
        scope_note, "",
        f"For each module and {unit_label}, whether binding-site switching (motif gained/lost "
        "between the switch-pair isoforms) is over-represented among the module's genes vs "
        f"the region's switch-gene pool. Two arms are reported for every module x {unit_label} "
        "cell:", "",
        "1. **Hypergeometric** (legacy) -- over-representation against the switch-gene pool, "
        "BH-corrected across all cells. It conditions on nothing else.",
        "2. **Covariate-adjusted binomial GLM** (reviewer item 6c) -- the same contrast with "
        "transcript length, GC content, 5'UTR/CDS/3'UTR composition and `n_transcripts` as "
        "covariates, so a module cannot score simply because its genes are long, GC-rich or "
        "UTR-heavy and therefore offer more motif *opportunity*. BH is applied over the "
        "estimable cells only.", "",
        "The adjusted arm is the stricter reading and the two can disagree in both "
        "directions; where they do, the disagreement is reported below rather than resolved "
        "in favour of the larger number.", "",
        f"- switch genes tested: **{calls['gene'].nunique()}** over "
        f"{calls['region'].nunique()} regions",
        f"- module x {unit_label} cells: **{len(reg)}**",
        f"- hypergeometric q<0.05: **{n_sig}** (GO-invisible {n_goinv})",
    ]

    if has_glm:
        status = reg["model_status"].value_counts() if "model_status" in reg.columns \
            else pd.Series(dtype=int)
        n_est = int(_num(reg, "pvalue_glm").notna().sum())
        n_glm_sig = int((_num(reg, "qvalue_glm") < 0.05).sum())
        s_or, s_lo, s_hi = _num(sig, "odds_ratio"), _num(sig, "ci_low"), _num(sig, "ci_high")
        n_sig_est = int(s_or.notna().sum())
        n_sig_glm = int((_num(sig, "qvalue_glm") < 0.05).sum())
        n_up = int((s_lo > 1).sum())
        contra = sig[(s_hi < 1).fillna(False)]
        lines += [
            f"- estimable GLM cells (BH family): **{n_est}** of {len(reg)}",
            f"- covariate-adjusted q<0.05: **{n_glm_sig}**",
            f"- of the {n_sig} hypergeometric hits: {n_sig_est} estimable, **{n_sig_glm}** "
            f"also adjusted-q<0.05, {n_up} with an adjusted CI entirely above 1, "
            f"**{len(contra)}** entirely *below* it (enriched on raw counts, depleted once "
            "opportunity is adjusted for)", "",
            "### Estimability", "",
            "A cell whose 2x2 has an empty margin admits no maximum-likelihood fit; including "
            "such cells would put an arbitrary point estimate into the BH family and dilute "
            "every real test. They are excluded from the correction and reported here instead "
            "of silently carrying a p-value.", "",
            "| GLM cell status | n | meaning |",
            "|---|---|---|",
        ]
        meaning = {
            "fit": "estimable; carries a p-value and enters the BH family",
            "separated_zero_cell": "an empty cell in the (in-module x switched) 2x2 -- "
                                   "complete separation, no finite MLE",
            "quasi_separated": "fitted probabilities pinned at 0/1, or a degenerate "
                               "coefficient/standard error",
            "outcome_or_predictor_constant": "no variation to model in the region",
        }
        for k, v in status.items():
            lines.append(f"| `{k}` | {v} | {meaning.get(k, '')} |")
        lines.append("")

    lines += [
        "### Top candidate regulons (by hypergeometric q)", "",
        f"| region | module | {unit_label} | module size | switched | enrichment | q | "
        "adj. OR (95% CI) | adj. q | GO-inv |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in sig.head(30).itertuples():
        name = getattr(r, unit)
        o = pd.to_numeric(getattr(r, "odds_ratio", np.nan), errors="coerce")
        lo = pd.to_numeric(getattr(r, "ci_low", np.nan), errors="coerce")
        hi = pd.to_numeric(getattr(r, "ci_high", np.nan), errors="coerce")
        qg = pd.to_numeric(getattr(r, "qvalue_glm", np.nan), errors="coerce")
        if np.isfinite(o) and np.isfinite(lo) and np.isfinite(hi):
            cell = f"{o:.2f} ({lo:.2f}-{hi:.2f})"
            if hi < 1:
                cell = f"**{cell}**"
        else:
            cell = f"not estimable ({getattr(r, 'model_status', 'n/a')})"
        qg_cell = f"{qg:.2e}" if np.isfinite(qg) else "-"
        lines.append(f"| {r.region} | {r.module_id} | {name} | {r.module_size} | "
                     f"{r.n_switched} | {r.enrichment:.2f} | {r.q:.2e} | {cell} | "
                     f"{qg_cell} | {'yes' if r.go_invisible else 'no'} |")

    if has_glm and len(contra):
        lines += [
            "", "### Contradicted by opportunity adjustment", "",
            f"These {len(contra)} cells are hypergeometric-significant yet have an adjusted "
            "95% CI entirely below 1 -- the raw enrichment is explained, and then some, by "
            "the motif opportunity their genes carry. They must not be described as "
            f"candidate regulons.", "",
            f"| region | module | {unit_label} | enrichment | q | adj. OR (95% CI) |",
            "|---|---|---|---|---|---|",
        ]
        for r in contra.sort_values("q").head(20).itertuples():
            o = pd.to_numeric(r.odds_ratio, errors="coerce")
            lo = pd.to_numeric(r.ci_low, errors="coerce")
            hi = pd.to_numeric(r.ci_high, errors="coerce")
            lines.append(f"| {r.region} | {r.module_id} | {getattr(r, unit)} | "
                         f"{r.enrichment:.2f} | {r.q:.2e} | {o:.2f} ({lo:.2f}-{hi:.2f}) |")
        lines.append("")

    out_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="Per-module RBP-regulon enrichment (stage 2).")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--scope", choices=("mature", "intronic", "combined"), default="mature",
                   help="which Stage-1 count table(s) to consume (default: mature, canonical)")
    p.add_argument("--unit", choices=("rbp", "family_id"), default="rbp",
                   help="count motifs per RBP (legacy) or per motif-similarity family "
                        "(deduplicated; see rbp_motif_families.py)")
    p.add_argument("--bg-mode", choices=("composition", "flat"), default="composition",
                   help="composition-matched background (default) or the legacy flat 0.25")
    p.add_argument("--region-class", choices=("5utr", "cds", "3utr", "noncoding"),
                   default=None, help="restrict motif hits to one transcript region class "
                                      "(default: whole mature transcript)")
    args = p.parse_args()
    run(args.fdr, args.scope, args.unit, args.bg_mode, args.region_class)


if __name__ == "__main__":
    main()
