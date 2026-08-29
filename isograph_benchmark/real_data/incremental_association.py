"""Switch-vs-abundance incremental association — strengths/limitations of IsoGraph.

Two complementary, reproducible analyses that ask whether isoform COMPOSITION
(IsoGraph's switch channel) carries phenotype signal beyond total ABUNDANCE,
using IsoGraph's own per-sample switch and abundance features (identical
preprocessing, so the contrast is apples-to-apples):

1. gene_level_deconfounded  — the rigorous, de-confounded test. Per gene:
     switch_g    ~ phenotype + abundance_g + covariates   (composition | abundance)
     abundance_g ~ phenotype + switch_g    + covariates   (abundance | composition)
   No clustering / eigengenes, so it removes the module-definition circularity.
   "composition-unique" genes (switch significant, abundance not) are genuine
   DTU-without-DGE: isoform-switching regulation a DGE/WGCNA pipeline cannot see.

2. module_level_incremental — the eigengene-level version (switch eigengene vs
   abundance eigengene of the SAME module genes). Reported for completeness/
   contrast: it INFLATES composition-uniqueness because modules are defined by
   switch co-expression (coherent switch eigengene, incoherent abundance
   eigengene). Keep it next to the gene-level result so the circularity is visible.

phenotype = SCZD diagnosis (brainseq-sczd) or non-linear age (brainseq-aging,
df=3 spline F-test). Outputs land in <artifact_dir>/incremental_association/.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from patsy import dmatrix
from scipy import stats

from isograph.io.artifacts import load_dataset_bundle

from isograph_benchmark.paths import ensure_dir
from isograph_benchmark.real_data.run_models import (
    GTEX_REGIONS,
    diagnosis_association,
    spline_age_association,
)
from isograph_benchmark.real_data.sweep_leiden import (
    AGING_COVARIATE_COLS,
    GTEX_AGE_COL,
    GTEX_COVARIATE_COLS,
    SCZD_COVARIATE_COLS,
    _artifact_dir,
    _bundle_path,
    _pivot_eigengenes,
)

_META_COLS = ("feature_id", "gene_id", "feature_type", "n_transcripts")


def _covariate_cols(analysis: str) -> list[str]:
    if analysis == "brainseq-sczd":
        return SCZD_COVARIATE_COLS
    if analysis == "gtex-aging":
        return GTEX_COVARIATE_COLS
    return AGING_COVARIATE_COLS


def _age_col(analysis: str) -> str:
    return GTEX_AGE_COL if analysis == "gtex-aging" else "Age"


def _load(analysis: str, region: str | None, variant: str):
    artifact_dir = _artifact_dir(analysis, region, variant)
    fs = pd.read_parquet(artifact_dir / "feature_scores.parquet")
    modules = pd.read_parquet(artifact_dir / "modules.parquet")
    bundle = load_dataset_bundle(_bundle_path(analysis, region))
    return artifact_dir, fs, modules, bundle


def _sample_cols(fs: pd.DataFrame, sample_ids: set[str]) -> list[str]:
    return [c for c in fs.columns if c not in _META_COLS and c in sample_ids]


def _channel_matrix(fs: pd.DataFrame, feature_type: str, samp: list[str]) -> pd.DataFrame:
    """gene_id-indexed (genes × samples) matrix for one feature channel."""
    sub = fs[fs["feature_type"] == feature_type].set_index("gene_id")
    return sub[samp]


def _eigengene_table(channel: pd.DataFrame, modules: pd.DataFrame, samp: list[str]) -> pd.DataFrame:
    """module_id × samples eigengene = mean of the channel over module genes
    (mirrors compute_trait_associations' switch eigengene construction)."""
    rows = {}
    for mid in sorted(modules["module_id"].unique()):
        genes = modules.loc[modules["module_id"] == mid, "gene_id"]
        present = channel.index.intersection(genes)
        if len(present) == 0:
            continue
        rows[mid] = channel.loc[present].to_numpy(float).mean(axis=0)
    eg = pd.DataFrame(rows, index=samp).T
    eg.index.name = "module_id"
    return eg.reset_index()


# --------------------------------------------------------------------------- #
# 1. Gene-level de-confounded test
# --------------------------------------------------------------------------- #
_OLS_KEYS = ("p", "beta", "se", "t", "dof", "partial_r2")
_FTEST_KEYS = ("p", "F", "dfn", "dfd", "partial_r2", "cohens_f2", "lrt")


def _ols_coef(X: np.ndarray, y: np.ndarray, idx: int) -> dict[str, float]:
    """Two-sided t-test on coefficient ``idx``, plus its effect size.

    ``partial_r2 = t² / (t² + dof)`` is the squared partial correlation of that single term
    — the share of residual variance in ``y`` it explains once every other column of ``X``
    is accounted for.  The p-value arithmetic is byte-for-byte what ``_ols_coef_p`` did.
    """
    out = dict.fromkeys(_OLS_KEYS, np.nan)
    if not (np.all(np.isfinite(X)) and np.all(np.isfinite(y))):
        return out
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    dof = len(y) - rank
    if dof <= 0:
        return out
    s2 = float((y - X @ beta) @ (y - X @ beta)) / dof
    se = np.sqrt(max(s2 * np.linalg.pinv(X.T @ X)[idx, idx], 1e-30))
    t = beta[idx] / se
    out.update(
        p=float(2 * stats.t.sf(abs(t), dof)),
        beta=float(beta[idx]), se=float(se), t=float(t), dof=float(dof),
        partial_r2=float(t * t / (t * t + dof)),
    )
    return out


def _ols_coef_p(X: np.ndarray, y: np.ndarray, idx: int) -> float:
    return _ols_coef(X, y, idx)["p"]


def _ftest_block(Xfull: np.ndarray, Xred: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """Nested F-test of the block in ``Xfull`` but not ``Xred``, plus its effect sizes.

    ``partial_r2 = (rss_r − rss_f) / rss_r`` is the fraction of the reduced model's residual
    variance the block explains; ``cohens_f2 = (rss_r − rss_f) / rss_f``; ``lrt`` is the
    Gaussian likelihood-ratio statistic ``n·log(rss_r / rss_f)``, asymptotically χ²(dfn).
    All three are read off quantities the F-test already computes, so the returned ``p`` is
    identical to what ``_ftest_block_p`` returned.
    """
    out = dict.fromkeys(_FTEST_KEYS, np.nan)
    if not (np.all(np.isfinite(Xfull)) and np.all(np.isfinite(y))):
        return out
    bf, _, rf, _ = np.linalg.lstsq(Xfull, y, rcond=None)
    br, _, rr, _ = np.linalg.lstsq(Xred, y, rcond=None)
    rss_f = float(((y - Xfull @ bf) ** 2).sum())
    rss_r = float(((y - Xred @ br) ** 2).sum())
    dfn, dfd = max(rf - rr, 1), max(len(y) - rf, 1)
    F = ((rss_r - rss_f) / dfn) / (rss_f / dfd)
    if not np.isfinite(F):
        return out
    out.update(
        p=float(stats.f.sf(F, dfn, dfd)), F=float(F), dfn=float(dfn), dfd=float(dfd),
        partial_r2=float((rss_r - rss_f) / rss_r) if rss_r > 0 else np.nan,
        cohens_f2=float((rss_r - rss_f) / rss_f) if rss_f > 0 else np.nan,
        lrt=float(len(y) * np.log(rss_r / rss_f)) if rss_f > 0 and rss_r > 0 else np.nan,
    )
    return out


def _ftest_block_p(Xfull: np.ndarray, Xred: np.ndarray, y: np.ndarray) -> float:
    return _ftest_block(Xfull, Xred, y)["p"]


def gene_level_deconfounded(
    analysis: str,
    fs: pd.DataFrame,
    bundle,
    fdr_alpha: float = 0.10,
    min_samples: int = 30,
    comp_cols: tuple[str, ...] = (),
) -> pd.DataFrame:
    st = bundle.sample_table
    sample_ids = set(st["sample_id"].astype(str))
    samp = _sample_cols(fs, sample_ids)
    st = st.set_index(st["sample_id"].astype(str)).loc[samp]
    covs = _covariate_cols(analysis) + list(comp_cols)

    SW = _channel_matrix(fs, "switch", samp)
    AB = _channel_matrix(fs, "abundance", samp)
    genes = SW.index.intersection(AB.index)
    SWm = SW.loc[genes].to_numpy(float).T   # samples × genes
    ABm = AB.loc[genes].to_numpy(float).T

    covdf = pd.get_dummies(st[[c for c in covs if c in st.columns]].copy(), drop_first=True)
    C = np.hstack([np.ones((len(st), 1)), covdf.to_numpy(dtype=float)])
    C_ok = np.all(np.isfinite(C), axis=1)

    pA = np.full(len(genes), np.nan)
    pB = np.full(len(genes), np.nan)
    # Effect sizes for the same two tests.  Binary FDR counts are sample-size dependent, so
    # the magnitude of the conditional signal is reported alongside them.
    effect_keys = _OLS_KEYS if analysis == "brainseq-sczd" else _FTEST_KEYS
    effA = {k: np.full(len(genes), np.nan) for k in effect_keys}
    effB = {k: np.full(len(genes), np.nan) for k in effect_keys}

    def _store(target: dict, j: int, res: dict) -> None:
        for key, value in res.items():
            target[key][j] = value

    if analysis == "brainseq-sczd":
        pheno = (st["Dx"].astype(str).to_numpy() == "SCZD").astype(float)[:, None]
        dx_idx = C.shape[1]
        for j in range(len(genes)):
            ys, ya = SWm[:, j], ABm[:, j]
            m = np.isfinite(ys) & np.isfinite(ya) & C_ok
            if m.sum() < min_samples:
                continue
            ra = _ols_coef(np.hstack([C[m], pheno[m], ya[m, None]]), ys[m], dx_idx)
            rb = _ols_coef(np.hstack([C[m], pheno[m], ys[m, None]]), ya[m], dx_idx)
            pA[j], pB[j] = ra["p"], rb["p"]
            _store(effA, j, ra)
            _store(effB, j, rb)
    else:
        age = st[_age_col(analysis)].to_numpy(float)
        age_z = (age - age.mean()) / age.std()
        knots = np.quantile(age_z, [0.5])
        bd = [age_z.min(), age_z.max()]
        B = np.asarray(dmatrix(
            f"cr(age_z, knots={list(knots)}, lower_bound={bd[0]}, upper_bound={bd[1]}) - 1",
            {"age_z": age_z}, return_type="matrix"))
        for j in range(len(genes)):
            ys, ya = SWm[:, j], ABm[:, j]
            m = np.isfinite(ys) & np.isfinite(ya) & C_ok
            if m.sum() < min_samples:
                continue
            ra = _ftest_block(np.hstack([C[m], B[m], ya[m, None]]), np.hstack([C[m], ya[m, None]]), ys[m])
            rb = _ftest_block(np.hstack([C[m], B[m], ys[m, None]]), np.hstack([C[m], ys[m, None]]), ya[m])
            pA[j], pB[j] = ra["p"], rb["p"]
            _store(effA, j, ra)
            _store(effB, j, rb)

    ok = np.isfinite(pA) & np.isfinite(pB)
    res = pd.DataFrame({
        "gene_id": genes[ok],
        "p_switch_given_abund": pA[ok],
        "p_abund_given_switch": pB[ok],
    })
    for key in effect_keys:
        if key == "p":
            continue
        res[f"{key}_switch_given_abund"] = effA[key][ok]
        res[f"{key}_abund_given_switch"] = effB[key][ok]
    res["fdr_switch_given_abund"] = stats.false_discovery_control(res["p_switch_given_abund"], method="bh")
    res["fdr_abund_given_switch"] = stats.false_discovery_control(res["p_abund_given_switch"], method="bh")
    sw = res["fdr_switch_given_abund"] <= fdr_alpha
    ab = res["fdr_abund_given_switch"] <= fdr_alpha
    res["category"] = np.select(
        [sw & ~ab, ab & ~sw, sw & ab], ["composition_unique", "abundance_unique", "both"], "neither"
    )
    return res.sort_values("p_switch_given_abund").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# 2. Module-level incremental (eigengene switch vs abundance) — for contrast
# --------------------------------------------------------------------------- #
def module_level_incremental(
    analysis: str,
    fs: pd.DataFrame,
    modules: pd.DataFrame,
    bundle,
    fdr_alpha: float = 0.10,
    comp_cols: tuple[str, ...] = (),
) -> pd.DataFrame:
    st = bundle.sample_table
    sample_ids = set(st["sample_id"].astype(str))
    samp = _sample_cols(fs, sample_ids)
    covs = _covariate_cols(analysis) + list(comp_cols)

    eg_switch = _pivot_eigengenes(_eigengene_table(_channel_matrix(fs, "switch", samp), modules, samp))
    eg_abund = _pivot_eigengenes(_eigengene_table(_channel_matrix(fs, "abundance", samp), modules, samp))

    def assoc(egp):
        if analysis == "brainseq-sczd":
            d = diagnosis_association(egp, st, covariate_cols=covs)
            return d.set_index("module_id")["fdr"]
        d = spline_age_association(egp, st, covariate_cols=covs, age_col=_age_col(analysis))
        return d.drop_duplicates("module_id").set_index("module_id")["fdr_ftest"]

    res = pd.DataFrame({"switch_fdr": assoc(eg_switch), "abund_fdr": assoc(eg_abund)}).dropna()
    res["n_genes"] = [int((modules["module_id"] == m).sum()) for m in res.index]
    sw = res["switch_fdr"] <= fdr_alpha
    ab = res["abund_fdr"] <= fdr_alpha
    res["category"] = np.select(
        [sw & ~ab, ab & ~sw, sw & ab], ["composition_unique", "abundance_unique", "both"], "neither"
    )
    return res.reset_index().rename(columns={"index": "module_id"})


_QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90]


def _effect_distribution(gene_level: pd.DataFrame) -> dict:
    """Quantiles of the conditional effect sizes, overall and per FDR category.

    An FDR count answers "how many genes cross a threshold", which is strongly
    sample-size dependent; the distribution answers "how much conditional signal is
    there".  ``category == "neither"`` is the internal noise floor these should be read
    against.
    """
    out: dict = {}
    for col in [c for c in gene_level.columns if c.startswith("partial_r2_")]:
        vals = pd.to_numeric(gene_level[col], errors="coerce")
        entry = {
            "n": int(vals.notna().sum()),
            "quantiles": {str(q): _q(vals, q) for q in _QUANTILES},
            "by_category": {},
        }
        for cat, grp in gene_level.groupby("category"):
            v = pd.to_numeric(grp[col], errors="coerce")
            entry["by_category"][str(cat)] = {
                "n": int(v.notna().sum()),
                "median": _q(v, 0.5),
                "q90": _q(v, 0.9),
            }
        out[col] = entry
    return out


def _q(values: pd.Series, q: float) -> float | None:
    v = values.dropna()
    return float(np.quantile(v, q)) if len(v) else None


def _summary(gene_level: pd.DataFrame, module_level: pd.DataFrame) -> dict:
    gc = gene_level["category"].value_counts().to_dict()
    mc = module_level["category"].value_counts().to_dict()
    return {
        "gene_level": {"n_tested": int(len(gene_level)), **{k: int(gc.get(k, 0)) for k in
                       ["composition_unique", "abundance_unique", "both", "neither"]}},
        "module_level": {"n_tested": int(len(module_level)), **{k: int(mc.get(k, 0)) for k in
                         ["composition_unique", "abundance_unique", "both", "neither"]}},
        "gene_level_effect_sizes": _effect_distribution(gene_level),
    }


# --------------------------------------------------------------------------- #
# 1b. Explicit nested ladder with the phenotype as response
# --------------------------------------------------------------------------- #
def gene_level_ladder(
    analysis: str,
    fs: pd.DataFrame,
    bundle,
    min_samples: int = 30,
    comp_cols: tuple[str, ...] = (),
) -> pd.DataFrame:
    """Incremental variance explained, with the phenotype as the response variable.

    Complements ``gene_level_deconfounded`` rather than replacing it.  That test puts the
    switch channel on the left-hand side and asks whether the phenotype term survives
    conditioning on abundance; this one is the literal ladder

        M0  : pheno ~ covariates
        M1a : pheno ~ covariates + abundance_g
        M1b : pheno ~ covariates + switch_g
        M2  : pheno ~ covariates + abundance_g + switch_g

    so ``delta_r2_switch_given_abund = R²(M2) − R²(M1a)`` is exactly "the additional
    phenotype-associated variation the switch coordinate contributes, conditional on
    abundance", and the reverse term is its abundance-side mirror.  Continuous traits use
    OLS R² and a nested F-test; the binary SCZD arm uses logistic McFadden pseudo-R² and a
    likelihood-ratio χ².
    """
    st = bundle.sample_table
    sample_ids = set(st["sample_id"].astype(str))
    samp = _sample_cols(fs, sample_ids)
    st = st.set_index(st["sample_id"].astype(str)).loc[samp]
    covs = _covariate_cols(analysis) + list(comp_cols)

    SW = _channel_matrix(fs, "switch", samp)
    AB = _channel_matrix(fs, "abundance", samp)
    genes = SW.index.intersection(AB.index)
    SWm = SW.loc[genes].to_numpy(float).T
    ABm = AB.loc[genes].to_numpy(float).T

    covdf = pd.get_dummies(st[[c for c in covs if c in st.columns]].copy(), drop_first=True)
    C = np.hstack([np.ones((len(st), 1)), covdf.to_numpy(dtype=float)])
    C_ok = np.all(np.isfinite(C), axis=1)

    binary = analysis == "brainseq-sczd"
    if binary:
        y_all = (st["Dx"].astype(str).to_numpy() == "SCZD").astype(float)
    else:
        y_all = st[_age_col(analysis)].to_numpy(float)

    rows = []
    for j in range(len(genes)):
        sw, ab = SWm[:, j], ABm[:, j]
        m = np.isfinite(sw) & np.isfinite(ab) & np.isfinite(y_all) & C_ok
        if m.sum() < min_samples:
            continue
        y = y_all[m]
        X0 = C[m]
        Xa = np.hstack([X0, ab[m, None]])
        Xs = np.hstack([X0, sw[m, None]])
        X2 = np.hstack([X0, ab[m, None], sw[m, None]])
        fit = _logistic_ladder if binary else _ols_ladder
        rows.append({"gene_id": genes[j], **fit(y, X0, Xa, Xs, X2)})

    res = pd.DataFrame(rows)
    if res.empty:
        return res
    for col in ("p_switch_given_abund", "p_abund_given_switch"):
        res[col.replace("p_", "fdr_", 1)] = stats.false_discovery_control(res[col], method="bh")
    return res.sort_values("p_switch_given_abund").reset_index(drop=True)


def _ols_ladder(y, X0, Xa, Xs, X2) -> dict[str, float]:
    tss = float(((y - y.mean()) ** 2).sum())

    def rss(X):
        b, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
        return float(((y - X @ b) ** 2).sum()), int(rank)

    r0, k0 = rss(X0)
    ra, ka = rss(Xa)
    rs, ks = rss(Xs)
    r2, k2 = rss(X2)
    n = len(y)

    def f_p(rss_r, rank_r, rss_f, rank_f):
        dfn, dfd = max(rank_f - rank_r, 1), max(n - rank_f, 1)
        F = ((rss_r - rss_f) / dfn) / (rss_f / dfd)
        return (float(stats.f.sf(F, dfn, dfd)) if np.isfinite(F) else np.nan,
                float(F) if np.isfinite(F) else np.nan)

    p_sw, f_sw = f_p(ra, ka, r2, k2)
    p_ab, f_ab = f_p(rs, ks, r2, k2)
    return {
        "r2_m0": 1 - r0 / tss if tss else np.nan,
        "r2_m1a_abund": 1 - ra / tss if tss else np.nan,
        "r2_m1b_switch": 1 - rs / tss if tss else np.nan,
        "r2_m2": 1 - r2 / tss if tss else np.nan,
        "delta_r2_switch_given_abund": (ra - r2) / tss if tss else np.nan,
        "delta_r2_abund_given_switch": (rs - r2) / tss if tss else np.nan,
        "lrt_switch_given_abund": float(n * np.log(ra / r2)) if r2 > 0 else np.nan,
        "lrt_abund_given_switch": float(n * np.log(rs / r2)) if r2 > 0 else np.nan,
        "f_switch_given_abund": f_sw, "f_abund_given_switch": f_ab,
        "p_switch_given_abund": p_sw, "p_abund_given_switch": p_ab,
        "n": float(n), "model": "ols",
    }


def _logistic_ladder(y, X0, Xa, Xs, X2) -> dict[str, float]:
    import statsmodels.api as sm

    def loglike(X):
        try:
            fit = sm.GLM(y, X, family=sm.families.Binomial()).fit()
            return float(fit.llf), int(np.linalg.matrix_rank(X))
        except Exception:
            return np.nan, 0

    ll_null, _ = loglike(np.ones((len(y), 1)))
    ll0, k0 = loglike(X0)
    lla, ka = loglike(Xa)
    lls, ks = loglike(Xs)
    ll2, k2 = loglike(X2)

    def mcfadden(ll):
        return 1 - ll / ll_null if np.isfinite(ll) and ll_null else np.nan

    def lrt(ll_r, rank_r, ll_f, rank_f):
        stat = 2 * (ll_f - ll_r)
        df = max(rank_f - rank_r, 1)
        return (float(stat), float(stats.chi2.sf(stat, df))) if np.isfinite(stat) else (np.nan, np.nan)

    lrt_sw, p_sw = lrt(lla, ka, ll2, k2)
    lrt_ab, p_ab = lrt(lls, ks, ll2, k2)
    return {
        "r2_m0": mcfadden(ll0), "r2_m1a_abund": mcfadden(lla),
        "r2_m1b_switch": mcfadden(lls), "r2_m2": mcfadden(ll2),
        "delta_r2_switch_given_abund": mcfadden(ll2) - mcfadden(lla),
        "delta_r2_abund_given_switch": mcfadden(ll2) - mcfadden(lls),
        "lrt_switch_given_abund": lrt_sw, "lrt_abund_given_switch": lrt_ab,
        "f_switch_given_abund": np.nan, "f_abund_given_switch": np.nan,
        "p_switch_given_abund": p_sw, "p_abund_given_switch": p_ab,
        "n": float(len(y)), "model": "logistic",
    }


def _attach_composition(artifact_dir, bundle) -> list[str]:
    """Merge written cell-type fractions into bundle.sample_table by sample_id and
    return the composition covariate columns. Requires celltype_composition to have run."""
    from isograph_benchmark.real_data.celltype_composition import load_composition_covariates
    frac, comp_cols = load_composition_covariates(artifact_dir)
    fr = frac[comp_cols].reset_index()
    fr.columns = ["sample_id", *comp_cols]
    fr["sample_id"] = fr["sample_id"].astype(str)
    st = bundle.sample_table.copy()
    st["sample_id"] = st["sample_id"].astype(str)
    bundle.sample_table = st.merge(fr, on="sample_id", how="left")
    return comp_cols


def run_analysis(analysis: str, region: str | None, variant: str, fdr_alpha: float,
                 composition: bool = False, ladder: bool = True) -> dict:
    label = f"{analysis}/{region}" if region else analysis
    artifact_dir, fs, modules, bundle = _load(analysis, region, variant)
    comp_cols = tuple(_attach_composition(artifact_dir, bundle)) if composition else ()
    if composition:
        print(f"[{label}] composition-adjusted; covariates += {list(comp_cols)}", flush=True)
    print(f"[{label}] gene-level de-confounded test ...", flush=True)
    gene_level = gene_level_deconfounded(analysis, fs, bundle, fdr_alpha=fdr_alpha,
                                         comp_cols=comp_cols)
    print(f"[{label}] module-level incremental ...", flush=True)
    module_level = module_level_incremental(analysis, fs, modules, bundle, fdr_alpha=fdr_alpha,
                                            comp_cols=comp_cols)

    subdir = "incremental_association_composition" if composition else "incremental_association"
    out = ensure_dir(artifact_dir / subdir)
    gene_level.to_parquet(out / "gene_level.parquet", index=False, compression="zstd")
    module_level.to_parquet(out / "module_level.parquet", index=False, compression="zstd")
    summary = _summary(gene_level, module_level)
    summary.update({"analysis": analysis, "region": region, "variant": variant,
                    "fdr_alpha": fdr_alpha, "composition_adjusted": composition,
                    "composition_covariates": list(comp_cols)})

    if ladder:
        print(f"[{label}] nested ladder (phenotype as response) ...", flush=True)
        lad = gene_level_ladder(analysis, fs, bundle, comp_cols=comp_cols)
        if not lad.empty:
            lad.to_parquet(out / "gene_level_ladder.parquet", index=False, compression="zstd")
            summary["gene_level_ladder"] = {
                "n_tested": int(len(lad)),
                "model": str(lad["model"].iloc[0]),
                **{
                    col: {"quantiles": {str(q): _q(lad[col], q) for q in _QUANTILES}}
                    for col in ("delta_r2_switch_given_abund", "delta_r2_abund_given_switch")
                },
            }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    g = summary["gene_level"]
    m = summary["module_level"]
    print(f"[{label}] GENE-LEVEL  (n={g['n_tested']}): composition_unique={g['composition_unique']} "
          f"abundance_unique={g['abundance_unique']} both={g['both']}")
    print(f"[{label}] MODULE-LEVEL(n={m['n_tested']}): composition_unique={m['composition_unique']} "
          f"abundance_unique={m['abundance_unique']} both={m['both']}")
    print(f"[{label}] written to {out}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis", choices=["brainseq-sczd", "brainseq-aging", "gtex-aging"])
    parser.add_argument("--region", action="append", dest="regions",
                        help="For brainseq-aging/gtex-aging: region(s). Repeatable. "
                             "Default: all brainseq (3) or all gtex (13) regions.")
    parser.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    parser.add_argument("--fdr", type=float, default=0.10)
    parser.add_argument("--composition", action="store_true",
                        help="add estimated cell-type fractions as inference covariates "
                             "(requires celltype_composition to have run; brainseq only). "
                             "Writes to incremental_association_composition/.")
    parser.add_argument("--no-ladder", dest="ladder", action="store_false",
                        help="skip the phenotype-as-response nested ladder "
                             "(gene_level_ladder.parquet)")
    args = parser.parse_args()

    if args.analysis == "brainseq-sczd":
        run_analysis("brainseq-sczd", None, args.variant, args.fdr, args.composition, args.ladder)
    elif args.analysis == "gtex-aging":
        for region in (args.regions or GTEX_REGIONS):
            run_analysis("gtex-aging", region, args.variant, args.fdr, args.composition,
                         args.ladder)
    else:
        for region in (args.regions or ["caudate", "hippocampus", "dlpfc"]):
            run_analysis("brainseq-aging", region, args.variant, args.fdr, args.composition,
                         args.ladder)


if __name__ == "__main__":
    main()
