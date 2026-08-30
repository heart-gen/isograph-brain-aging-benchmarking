"""Empirical null for the cross-cohort aging-replication count (reviewer item 2).

``module_trust replication`` reports how many trusted BrainSEQ modules have a GTEx
gene-overlap match whose age association is nominally significant in both cohorts and points
the same way.  That count is descriptive: nothing establishes what the complete
matching-and-testing procedure would produce by chance.  This module supplies the null.

Two nulls are computed, because they answer different questions and only the pair of them is
honest:

``age`` (Freedman-Lane)
    Keep module construction, the trusted sets, the BrainSEQ→GTEx gene-Jaccard matching and
    every BrainSEQ statistic **fixed**; regenerate only the GTEx age association under the
    age-null.  Residuals from the covariate-only model are permuted and added back to the
    covariate-only fit, so the covariate structure is preserved rather than broken by naively
    shuffling age.  Answers: "is there more cross-cohort age concordance than chance?"

``matching``
    Keep every age statistic fixed and permute **which** GTEx module each BrainSEQ module is
    matched to (within the observed match pool, and within a size-matched pool).  Answers the
    question a reviewer asks next: "is the *specific* gene-overlap match doing any work?"

A significant ``age`` null with a null ``matching`` null means the cohorts agree that module-
level age association exists, but not on which modules — which is a materially weaker claim
than "replication" and must be written as such.

No model is re-fitted.  Real-data fits do not persist an eigengene table, but the eigengene
is the plain mean of ``feature_scores`` rows over a module's genes
(``isograph.models.base.compute_trait_associations``), and both ``feature_scores.parquet``
and ``modules.parquet`` are on disk — so it is reconstructed exactly and checked against the
published ``age_linear.parquet`` before any permutation runs.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from patsy import dmatrix
from scipy import stats

from isograph_benchmark.paths import rel
from isograph_benchmark.real_data.module_trust import (
    METHOD_DIRS,
    PROD_ROOTS,
    REGION_PAIRS,
    _crosscohort_rows,
    _out_dir,
)
from isograph_benchmark.real_data.run_models import (
    BRAINSEQ_COVARIATES,
    BRAINSEQ_DISCOVERY_COVARIATES,
    GTEX_COVARIATES,
    GTEX_DISCOVERY_COVARIATES,
)

SEED = 13
N_PERMUTATIONS = 10_000
SIG = 0.05

# GTEx cohort side of each pair; the BrainSEQ side is never permuted.
GTEX_BUNDLES = {
    "caudate_basal_ganglia": ("inputs", "bundles", "gtex_v11_brain", "caudate_basal_ganglia"),
    "hippocampus": ("inputs", "bundles", "gtex_v11_brain", "hippocampus"),
    "frontal_cortex_ba9": ("inputs", "bundles", "gtex_v11_brain", "frontal_cortex_ba9"),
}
GTEX_AGE_COL = "AGE"

# BrainSEQ discovery side.  The statistic must be computed the SAME way on both cohorts:
# scoring BrainSEQ with the published Pearson while scoring GTEx with a spline would make
# `both_sig` a hybrid of two different age models, and its permutation p would calibrate a
# statistic the manuscript never reports.
BRAINSEQ_BUNDLES = {
    "caudate": ("inputs", "bundles", "brainseq_v1", "caudate"),
    "hippocampus": ("inputs", "bundles", "brainseq_v1", "hippocampus"),
    "dlpfc": ("inputs", "bundles", "brainseq_v1", "dlpfc"),
}
BRAINSEQ_AGE_COL = "Age"

STATISTICS = ("pearson", "partial_linear", "spline_f")

# How much covariate adjustment the *age model* should apply on top of whatever the fit
# already did.  IsoGraph residualizes the DISCOVERY covariates inside the fit when building
# the network, so adjusting for those again downstream would be a double adjustment — but
# only if the persisted feature_scores carry the residualized values.  Measure before
# choosing: regress the discovery covariates on the reconstructed eigengenes and look at R².
# Near 0 means the features are residualized and `none`/`complement` is right; substantially
# above 0 means they are raw and `full` is the single, correct adjustment.
COVARIATE_MODES = ("full", "complement", "none")


def _covariates_for(cohort: str, mode: str) -> list[str]:
    """The covariate list the age model should use, given what the fit already removed."""
    full = BRAINSEQ_COVARIATES if cohort == "brainseq" else GTEX_COVARIATES
    discovery = (BRAINSEQ_DISCOVERY_COVARIATES if cohort == "brainseq"
                 else GTEX_DISCOVERY_COVARIATES)
    if mode == "full":
        return list(full)
    if mode == "none":
        return []
    if mode == "complement":
        # Only the covariates the fit did NOT residualize, so each is adjusted exactly once.
        return [c for c in full if c not in set(discovery)]
    raise SystemExit(f"unknown covariate mode {mode!r}; choose from {COVARIATE_MODES}")


def covariate_leakage(cohort: str, region: str, method: str, sample_table: pd.DataFrame,
                      age_col: str) -> float:
    """Median R² of the *discovery* covariates on the eigengenes.

    The diagnostic that decides which covariate mode is correct: if the fit's residualization
    reached the persisted feature_scores, these covariates explain ~no eigengene variance.
    """
    discovery = (BRAINSEQ_DISCOVERY_COVARIATES if cohort == "brainseq"
                 else GTEX_DISCOVERY_COVARIATES)
    eig = load_eigengenes(cohort, region, method)
    samples = [s for s in eig.columns if s in set(sample_table["sample_id"].astype(str))]
    X0, age, _ = _design(sample_table, samples, age_col, discovery)
    keep = np.isfinite(age) & np.all(np.isfinite(X0), axis=1)
    Y, X = eig[samples].to_numpy(float).T[keep], X0[keep]
    coef, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ coef
    ss_tot = ((Y - Y.mean(axis=0)) ** 2).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = 1 - (resid ** 2).sum(axis=0) / np.where(ss_tot > 0, ss_tot, np.nan)
    return float(np.nanmedian(r2))


# --------------------------------------------------------------------------- #
# Eigengene reconstruction
# --------------------------------------------------------------------------- #
def load_eigengenes(cohort: str, region: str, method: str) -> pd.DataFrame:
    """Module × sample eigengene matrix, from persisted artifacts.

    For the IsoGraph backends this is a reconstruction: ``compute_trait_associations``
    defines the eigengene as the unweighted mean of every ``feature_scores`` row (both the
    switch and abundance channels) belonging to the module's genes, and both inputs are on
    disk.  WGCNA's eigengene is a first principal component and is *not* recoverable that
    way, so its R runner persists ``eigengenes.parquet`` directly and we read it.  Either
    way ``verify_reconstruction`` checks the matrix against the published
    ``age_linear.parquet`` before any permutation runs.
    """
    root = PROD_ROOTS[(cohort, region)]
    art = rel(*root, METHOD_DIRS[method])

    stored = art / "eigengenes.parquet"
    if stored.exists():
        eig = pd.read_parquet(stored)
        if "sample_id" not in eig.columns:
            raise SystemExit(f"{stored} has no sample_id column")
        eig = eig.set_index(eig["sample_id"].astype(str)).drop(columns=["sample_id"])
        return eig.T.sort_index()          # modules × samples, matching the IsoGraph branch

    fs_path = art / "feature_scores.parquet"
    if not fs_path.exists():
        raise SystemExit(
            f"{art} carries neither eigengenes.parquet nor feature_scores.parquet, so its "
            f"eigengenes cannot be recovered. For the classical WGCNA baseline, re-run "
            f"real_data/gtex/_h/02.wgcna_gene.sh (or brainseq 01.wgcna_gene_aging.sh), "
            f"which now persists eigengenes.parquet."
        )
    fs = pd.read_parquet(fs_path)
    modules = pd.read_parquet(art / "modules.parquet")

    meta = {"feature_id", "gene_id", "feature_type", "n_transcripts"}
    samples = [c for c in fs.columns if c not in meta]
    fs["gene_id"] = fs["gene_id"].astype(str)
    modules["gene_id"] = modules["gene_id"].astype(str)
    modules["module_id"] = modules["module_id"].astype(str)

    gene_to_module = dict(zip(modules["gene_id"], modules["module_id"]))
    fs = fs[fs["gene_id"].isin(gene_to_module)].copy()
    fs["module_id"] = fs["gene_id"].map(gene_to_module)
    eig = fs.groupby("module_id")[samples].mean()
    return eig.sort_index()


def verify_reconstruction(cohort: str, region: str, method: str, sample_table: pd.DataFrame,
                          age_col: str, tol: float = 1e-6) -> float:
    """Assert the reconstructed eigengenes reproduce the published ``age_linear.parquet``.

    Returns the max absolute deviation in the Pearson effect.  A failure here means the
    permutation would be calibrating a different quantity than the published statistic.
    """
    root = PROD_ROOTS[(cohort, region)]
    published = pd.read_parquet(rel(*root, METHOD_DIRS[method], "age_linear.parquet"))
    published = published.set_index(published["module_id"].astype(str))

    eig = load_eigengenes(cohort, region, method)
    st = sample_table.set_index(sample_table["sample_id"].astype(str))
    shared = [s for s in eig.columns if s in st.index]
    age = pd.to_numeric(st.loc[shared, age_col], errors="coerce").to_numpy(float)

    worst = 0.0
    for module_id, row in eig[shared].iterrows():
        if module_id not in published.index:
            continue
        y = row.to_numpy(float)
        mask = np.isfinite(y) & np.isfinite(age)
        if mask.sum() < 10:
            continue
        r, _ = stats.pearsonr(age[mask], y[mask])
        worst = max(worst, abs(r - float(published.loc[module_id, "effect"])))
    if worst > tol:
        raise SystemExit(
            f"eigengene reconstruction for {cohort}/{region} deviates from the published "
            f"age_linear.parquet by {worst:.3g} (> {tol:g}); the permutation would not be "
            f"calibrating the published statistic."
        )
    return worst


# --------------------------------------------------------------------------- #
# Age models (vectorised over modules; one design, many responses)
# --------------------------------------------------------------------------- #
def _design(sample_table: pd.DataFrame, samples: list[str], age_col: str,
            covariates: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (X0 covariate-only, age vector, spline basis) aligned to ``samples``."""
    st = sample_table.set_index(sample_table["sample_id"].astype(str)).loc[samples]
    age = pd.to_numeric(st[age_col], errors="coerce").to_numpy(float)

    cols = [c for c in covariates if c in st.columns]
    # An empty covariate list is legitimate (covariates="none"): the reduced model is then
    # the intercept alone, and get_dummies of no columns cannot be concatenated.
    X0 = np.ones((len(st), 1))
    if cols:
        cov = pd.get_dummies(st[cols].copy(), drop_first=True).to_numpy(float)
        X0 = np.hstack([X0, cov])

    age_z = (age - np.nanmean(age)) / np.nanstd(age)
    knots = np.nanquantile(age_z, [0.5])
    basis = np.asarray(dmatrix(
        f"cr(age_z, knots={list(knots)}, lower_bound={np.nanmin(age_z)}, "
        f"upper_bound={np.nanmax(age_z)}) - 1",
        {"age_z": age_z}, return_type="matrix"))
    return X0, age, basis


def _ols(X: np.ndarray, Y: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """Least squares of every column of Y on X.  Returns (coef, residuals, rank)."""
    coef, _, rank, _ = np.linalg.lstsq(X, Y, rcond=None)
    return coef, Y - X @ coef, int(rank)


def age_statistics(Y: np.ndarray, X0: np.ndarray, age: np.ndarray, basis: np.ndarray,
                   statistic: str) -> tuple[np.ndarray, np.ndarray]:
    """Per-module (p_value, signed effect) for one of the three age models.

    ``Y`` is samples × modules.  ``pearson`` reproduces the published statistic exactly
    (raw correlation, no covariates); ``partial_linear`` is the covariate-adjusted signed
    age term; ``spline_f`` is the df=3 spline block F-test with the sign taken from the
    covariate-adjusted linear term (the F-test is unsigned, so it cannot define concordance
    on its own).
    """
    n = Y.shape[0]
    if statistic == "pearson":
        a = age - age.mean()
        Yc = Y - Y.mean(axis=0, keepdims=True)
        denom = np.sqrt((a @ a) * (Yc * Yc).sum(axis=0))
        r = np.where(denom > 0, (a @ Yc) / np.where(denom > 0, denom, 1.0), 0.0)
        r = np.clip(r, -1 + 1e-12, 1 - 1e-12)
        dof = n - 2
        t = r * np.sqrt(dof / (1 - r**2))
        return 2 * stats.t.sf(np.abs(t), dof), r

    Xfull = np.hstack([X0, age[:, None]])
    coef, resid, rank = _ols(Xfull, Y)
    dof = max(n - rank, 1)
    sigma2 = (resid * resid).sum(axis=0) / dof
    var_beta = np.linalg.pinv(Xfull.T @ Xfull)[-1, -1]
    beta = coef[-1]
    se = np.sqrt(np.maximum(sigma2 * var_beta, 1e-300))
    t = beta / se
    p_linear = 2 * stats.t.sf(np.abs(t), dof)
    if statistic == "partial_linear":
        return p_linear, beta

    # spline_f: joint F-test of the spline block against the covariate-only model.
    Xs = np.hstack([X0, basis])
    _, resid_s, rank_s = _ols(Xs, Y)
    _, resid_0, rank_0 = _ols(X0, Y)
    rss_f = (resid_s * resid_s).sum(axis=0)
    rss_r = (resid_0 * resid_0).sum(axis=0)
    dfn = max(rank_s - rank_0, 1)
    dfd = max(n - rank_s, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        F = ((rss_r - rss_f) / dfn) / (rss_f / dfd)
    p = np.where(np.isfinite(F), stats.f.sf(F, dfn, dfd), np.nan)
    return p, beta


# --------------------------------------------------------------------------- #
# The replication statistic
# --------------------------------------------------------------------------- #
def _observed_rows(method: str, k: int, sig: float) -> pd.DataFrame:
    """The fixed 130-row matching table, straight from module_trust (never permuted)."""
    frames = [pd.DataFrame(_crosscohort_rows(pair, method, k, sig)) for pair in REGION_PAIRS]
    return pd.concat(frames, ignore_index=True)


def _replication_count(bs_eff, bs_p, gt_eff, gt_p, sig: float) -> int:
    ok = np.isfinite(bs_eff) & np.isfinite(gt_eff) & np.isfinite(bs_p) & np.isfinite(gt_p)
    return int((ok & (bs_p < sig) & (gt_p < sig) & (np.sign(bs_eff) == np.sign(gt_eff))).sum())


def run(method: str, statistic: str, null: str, n_perm: int, seed: int, sig: float,
        min_jaccard: float, k: int, covariates: str = "full") -> None:
    from isograph.io.artifacts import load_dataset_bundle

    rows = _observed_rows(method, k, sig)
    rows = rows[rows["gene_jaccard"] >= min_jaccard].reset_index(drop=True)
    print(f"[{method}/{statistic}/{null}] matched rows at jaccard>={min_jaccard}: {len(rows)}",
          flush=True)

    # ---- GTEx side: reconstruct eigengenes, verify, build designs -------------------
    gtex_state: dict[str, dict] = {}
    for pair, (_, (gc, gr)) in REGION_PAIRS.items():
        bundle = load_dataset_bundle(rel(*GTEX_BUNDLES[gr]))
        st = bundle.sample_table
        dev = verify_reconstruction(gc, gr, method, st, GTEX_AGE_COL)
        eig = load_eigengenes(gc, gr, method)
        samples = [s for s in eig.columns if s in set(st["sample_id"].astype(str))]
        X0, age, basis = _design(st, samples, GTEX_AGE_COL,
                                 _covariates_for("gtex", covariates))
        keep = _keep_mask(statistic, age, X0)
        Y = eig[samples].to_numpy(float).T[keep]
        gtex_state[pair] = {
            "modules": list(eig.index), "Y": Y, "X0": X0[keep],
            "age": age[keep], "basis": basis[keep],
        }
        print(f"  [{pair}] {gc}/{gr}: {Y.shape[0]} samples × {Y.shape[1]} modules "
              f"(eigengene check max|Δr| = {dev:.2e})", flush=True)

    # ---- BrainSEQ side: same statistic, same code path, never permuted ---------------
    bs_eff, bs_p = _brainseq_vectors(rows, method, statistic, covariates)
    rows["age_effect_bs"] = bs_eff
    rows["age_p_bs"] = bs_p

    gt_eff_obs, gt_p_obs = _gtex_vectors(rows, gtex_state, statistic, permuted=None)
    t_obs = _replication_count(bs_eff, bs_p, gt_eff_obs, gt_p_obs, sig)
    per_pair_obs = {
        pair: _replication_count(*_subset(rows, pair, bs_eff, bs_p, gt_eff_obs, gt_p_obs), sig)
        for pair in REGION_PAIRS
    }
    print(f"  T_obs = {t_obs}  {per_pair_obs}", flush=True)

    # ---- null ------------------------------------------------------------------------
    rng = np.random.default_rng(seed)
    t_null = np.empty(n_perm, dtype=int)
    per_pair_null = {pair: np.empty(n_perm, dtype=int) for pair in REGION_PAIRS}

    if null == "age":
        # Freedman-Lane: residuals of the covariate-only (age-null) model are permuted and
        # added back to that model's fit, so covariate structure survives the permutation.
        for pair, state in gtex_state.items():
            _, resid, _ = _ols(state["X0"], state["Y"])
            state["fitted0"] = state["Y"] - resid
            state["resid0"] = resid
        for b in range(n_perm):
            permuted = {}
            for pair, state in gtex_state.items():
                idx = rng.permutation(state["resid0"].shape[0])
                permuted[pair] = state["fitted0"] + state["resid0"][idx]
            gt_eff, gt_p = _gtex_vectors(rows, gtex_state, statistic, permuted=permuted)
            t_null[b] = _replication_count(bs_eff, bs_p, gt_eff, gt_p, sig)
            for pair in REGION_PAIRS:
                per_pair_null[pair][b] = _replication_count(
                    *_subset(rows, pair, bs_eff, bs_p, gt_eff, gt_p), sig)
            _progress(b, n_perm)
    elif null == "matching":
        # Age statistics stay fixed; the BrainSEQ->GTEx assignment is reshuffled within each
        # pair's observed GTEx module pool.
        pool = {pair: np.flatnonzero((rows["pair"] == pair).to_numpy()) for pair in REGION_PAIRS}
        gtex_by_pair = {
            pair: _pair_gtex_table(gtex_state[pair], statistic) for pair in REGION_PAIRS
        }
        for b in range(n_perm):
            gt_eff = np.full(len(rows), np.nan)
            gt_p = np.full(len(rows), np.nan)
            for pair, idx in pool.items():
                eff_pool, p_pool = gtex_by_pair[pair]
                pick = rng.integers(0, len(eff_pool), size=len(idx))
                gt_eff[idx] = eff_pool[pick]
                gt_p[idx] = p_pool[pick]
            t_null[b] = _replication_count(bs_eff, bs_p, gt_eff, gt_p, sig)
            for pair in REGION_PAIRS:
                per_pair_null[pair][b] = _replication_count(
                    *_subset(rows, pair, bs_eff, bs_p, gt_eff, gt_p), sig)
            _progress(b, n_perm)
    else:
        raise SystemExit(f"unknown null: {null!r}")

    p_emp = float((1 + int((t_null >= t_obs).sum())) / (n_perm + 1))
    null_mean, null_sd = float(t_null.mean()), float(t_null.std(ddof=1))
    stats_out = {
        "method": method, "statistic": statistic, "null": null,
        "covariates": covariates,
        "n_matched_rows": int(len(rows)), "min_jaccard": min_jaccard,
        "T_obs": int(t_obs), "per_pair_obs": per_pair_obs,
        "B": int(n_perm), "seed": int(seed), "sig": sig,
        "p_emp": p_emp, "null_mean": null_mean, "null_sd": null_sd,
        "z": float((t_obs - null_mean) / null_sd) if null_sd > 0 else None,
        "null_max": int(t_null.max()), "null_q95": float(np.quantile(t_null, 0.95)),
    }
    print(f"  T_obs={t_obs} null={null_mean:.2f}±{null_sd:.2f} p_emp={p_emp:.4g}", flush=True)

    out_dir = _out_dir()
    tag = f"{method}__{statistic}__{null}"
    if covariates != "full":
        tag += f"__cov-{covariates}"
    draws = pd.DataFrame({"b": np.arange(n_perm), "T": t_null,
                          **{f"T_{p}": v for p, v in per_pair_null.items()}})
    draws.to_parquet(out_dir / f"replication_permutation__{tag}.parquet", index=False,
                     compression="zstd")
    (out_dir / f"replication_permutation__{tag}__stats.json").write_text(
        json.dumps(stats_out, indent=2))
    print(f"  wrote replication_permutation__{tag}.*", flush=True)


def _progress(b: int, n_perm: int) -> None:
    if (b + 1) % 2000 == 0:
        print(f"    {b + 1:,}/{n_perm:,}", flush=True)


def _subset(rows, pair, bs_eff, bs_p, gt_eff, gt_p):
    m = (rows["pair"] == pair).to_numpy()
    return bs_eff[m], bs_p[m], gt_eff[m], gt_p[m]


def _keep_mask(statistic: str, age: np.ndarray, X0: np.ndarray) -> np.ndarray:
    """Which samples enter the fit.

    ``pearson`` is the *published* covariate-free statistic, so it must be computed on the
    samples the published fit used — everything with a finite age.  Dropping samples with an
    incomplete covariate row there would silently change the reported number.  The
    covariate-adjusted statistics genuinely cannot use those samples, so they filter.
    """
    ok_age = np.isfinite(age)
    if statistic == "pearson":
        return ok_age
    return ok_age & np.all(np.isfinite(X0), axis=1)


def _brainseq_vectors(rows: pd.DataFrame, method: str, statistic: str,
                      covariates: str) -> tuple[np.ndarray, np.ndarray]:
    """Per-row BrainSEQ age effect and p under the *same* statistic used on GTEx.

    Computed from BrainSEQ's own eigengenes rather than read from ``age_linear.parquet``,
    so that ``spline_f`` really is a spline test on both sides.  The BrainSEQ side is never
    permuted — it is the fixed discovery half of the concordance statistic.
    """
    from isograph.io.artifacts import load_dataset_bundle

    eff = np.full(len(rows), np.nan)
    pval = np.full(len(rows), np.nan)
    for pair, ((bc, br), _) in REGION_PAIRS.items():
        bundle = load_dataset_bundle(rel(*BRAINSEQ_BUNDLES[br]))
        st = bundle.sample_table
        verify_reconstruction(bc, br, method, st, BRAINSEQ_AGE_COL)
        eig = load_eigengenes(bc, br, method)
        samples = [s for s in eig.columns if s in set(st["sample_id"].astype(str))]
        X0, age, basis = _design(st, samples, BRAINSEQ_AGE_COL,
                                 _covariates_for("brainseq", covariates))
        keep = _keep_mask(statistic, age, X0)
        Y = eig[samples].to_numpy(float).T[keep]
        p, e = age_statistics(Y, X0[keep], age[keep], basis[keep], statistic)
        index = {m: i for i, m in enumerate(eig.index)}
        for row_i in np.flatnonzero((rows["pair"] == pair).to_numpy()):
            j = index.get(str(rows["bs_module"].iloc[row_i]))
            if j is not None:
                eff[row_i] = e[j]
                pval[row_i] = p[j]
    return eff, pval


def _pair_gtex_table(state: dict, statistic: str) -> tuple[np.ndarray, np.ndarray]:
    p, eff = age_statistics(state["Y"], state["X0"], state["age"], state["basis"], statistic)
    return eff, p


def _gtex_vectors(rows: pd.DataFrame, gtex_state: dict, statistic: str,
                  permuted: dict | None) -> tuple[np.ndarray, np.ndarray]:
    """Map each matched row's GTEx module onto its age effect and p under the current draw."""
    eff = np.full(len(rows), np.nan)
    pval = np.full(len(rows), np.nan)
    for pair, state in gtex_state.items():
        Y = state["Y"] if permuted is None else permuted[pair]
        p, e = age_statistics(Y, state["X0"], state["age"], state["basis"], statistic)
        index = {m: i for i, m in enumerate(state["modules"])}
        mask = (rows["pair"] == pair).to_numpy()
        for row_i in np.flatnonzero(mask):
            j = index.get(str(rows["gtex_match"].iloc[row_i]))
            if j is not None:
                eff[row_i] = e[j]
                pval[row_i] = p[j]
    return eff, pval


_STATISTIC_LABEL = {
    "pearson": "Pearson r, no covariates (the published statistic)",
    "partial_linear": "linear age term, covariate-adjusted",
    "spline_f": "df=3 natural cubic spline block F-test, covariate-adjusted",
}


def report() -> None:
    """Collect every stats json into REPLICATION_PERMUTATION.md."""
    out_dir = _out_dir()
    rows = []
    for p in sorted(out_dir.iterdir()):
        if p.name.startswith("replication_permutation__") and p.name.endswith("__stats.json"):
            rows.append(json.loads(p.read_text()))
    if not rows:
        raise SystemExit(f"no replication_permutation__*__stats.json in {out_dir}")
    df = pd.DataFrame(rows)
    if "covariates" not in df.columns:
        df["covariates"] = "full"
    df["covariates"] = df["covariates"].fillna("full")

    lines = [
        "# Empirical null for the cross-cohort aging-replication count", "",
        "`T_obs` counts matched BrainSEQ<->GTEx module pairs whose age association is "
        "nominally significant (p<0.05) in **both** cohorts with a concordant sign. The two "
        "nulls answer different questions and both are reported: **age** (Freedman-Lane) "
        "regenerates the GTEx age association under the age-null while holding module "
        "construction, the trusted sets and the gene-Jaccard matching fixed; **matching** "
        "holds every age statistic fixed and permutes which GTEx module each BrainSEQ module "
        "is matched to. `matching` is the stricter of the two.", "",
        "The statistic is computed **the same way on both cohorts**. Scoring BrainSEQ with "
        "the published covariate-free Pearson while scoring GTEx with a covariate-adjusted "
        "model would make `both_sig` a hybrid of two age models and its p-value would "
        "calibrate a statistic that is never reported.", "",
    ]
    lines += [
        "Covariate modes: `full` adds every covariate; `complement` adds only those the fit "
        "did not already residualize, so each is adjusted exactly once; `none` adds nothing. "
        "IsoGraph residualizes its DISCOVERY covariates inside the fit, so `complement` is "
        "the non-double-adjusting choice there — but only if the persisted feature_scores "
        "carry residualized values, which `covariate_leakage()` measures. WGCNA's eigengenes "
        "are not residualized at all, so it needs `full`.", "",
    ]
    for method, g in df.groupby("method"):
        lines += [f"## `{method}`", "",
                  "| covariates | statistic | T_obs | n pairs | null | null mean ± sd | p_emp |",
                  "|---|---|---|---|---|---|---|"]
        for cov in ("complement", "none", "full"):
            for stat in STATISTICS:
                for null in ("age", "matching"):
                    r = g[(g["statistic"] == stat) & (g["null"] == null)
                          & (g["covariates"] == cov)]
                    if r.empty:
                        continue
                    r = r.iloc[0]
                    lines.append(
                        f"| {cov} | {stat} | {int(r['T_obs'])} | {int(r['n_matched_rows'])} | "
                        f"{null} | {r['null_mean']:.2f} ± {r['null_sd']:.2f} | "
                        f"{r['p_emp']:.4g} |")
        lines.append("")
        lines += ["Statistics: " + "; ".join(f"`{k}` = {v}" for k, v in _STATISTIC_LABEL.items()),
                  ""]

    lines += ["## How to write this up", ""]
    iso = df[df["method"] == "isograph"]
    prim = iso[(iso["statistic"] == "spline_f") & (iso["null"] == "matching")
               & (iso["covariates"] == "complement")]
    pear = iso[(iso["statistic"] == "pearson") & (iso["null"] == "matching")
               & (iso["covariates"] == "complement")]
    if not prim.empty and not pear.empty:
        prim, pear = prim.iloc[0], pear.iloc[0]
        ok = prim["p_emp"] < 0.05
        lines.append(
            f"The primary (covariate-adjusted spline) count is "
            f"**{int(prim['T_obs'])}/{int(prim['n_matched_rows'])}**, p_emp="
            f"{prim['p_emp']:.4g} against the matching null. The published covariate-free "
            f"Pearson count is **{int(pear['T_obs'])}/{int(pear['n_matched_rows'])}**, "
            f"p_emp={pear['p_emp']:.4g}.")
        lines += ["", (
            "The count survives the matching null under the primary model, so a subset of "
            "matched modules may be described as showing greater cross-cohort age "
            "concordance than expected by chance."
        ) if ok else (
            "**The count does not survive under the primary model.** Per the pre-registered "
            "decision rule the word \"replication\" must not be used; the honest wording is "
            "\"matched modules with concordant age effects\", reported alongside the "
            "covariate-free sensitivity analysis and the note that the two disagree."
        )]
    (out_dir / "REPLICATION_PERMUTATION.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {out_dir / 'REPLICATION_PERMUTATION.md'} ({len(df)} cells)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", action="store_true",
                    help="collect existing stats json files into the markdown report and exit")
    ap.add_argument("--method", default="isograph", choices=["isograph", "wgcna"])
    ap.add_argument("--statistic", default="partial_linear", choices=STATISTICS)
    ap.add_argument("--null", default="age", choices=["age", "matching"])
    ap.add_argument("--n-perm", type=int, default=N_PERMUTATIONS)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--sig", type=float, default=SIG)
    ap.add_argument("--min-jaccard", type=float, default=0.0,
                    help="sensitivity: restrict to matches at or above this gene Jaccard")
    ap.add_argument("--k", type=int, default=5, help="top-k driver transcripts (matching is "
                                                     "by gene set; k only affects reporting)")
    ap.add_argument("--covariates", default="full", choices=list(COVARIATE_MODES),
                    help="covariate adjustment the age model applies on top of the fit: "
                         "full = every covariate; complement = only those the fit did not "
                         "residualize (use when feature_scores carry residualized values, so "
                         "each covariate is adjusted exactly once); none = no adjustment")
    args = ap.parse_args()
    if args.report:
        report()
        return
    run(args.method, args.statistic, args.null, args.n_perm, args.seed, args.sig,
        args.min_jaccard, args.k, args.covariates)


if __name__ == "__main__":
    main()
