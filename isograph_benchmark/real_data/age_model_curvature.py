"""Is the module-level Age trajectory nonlinear? A direct curvature test.

Why this exists
---------------
Two age models are fitted for every module: a covariate-free linear (Pearson) arm and a
covariate-adjusted df=3 natural cubic spline.  They disagree about how many cross-cohort
aging replications survive (`module_trust replication`), and the manuscript has to say
which one it quotes.  That choice is only defensible if we know *why* they disagree.

There are two possible reasons, and they call for opposite decisions:

1. **The spline sees curvature the line misses.**  Then the spline is the better model,
   the linear arm is misspecified, and the lower replication count is the honest one.
2. **There is no curvature, and the extra 2 df are spent on noise.**  Then the spline is
   over-parameterized for these data, its lower count is a power tax rather than a
   correction, and the linear arm is the one to quote -- with the spline shown as a
   declared sensitivity rather than a competing result.

Reporting the two counts side by side without resolving this leaves a reviewer to guess,
and the natural guess is (1) -- the unfavourable one.  So test it directly.

The test
--------
The spline fit already writes, for every module, its projected Age effect at k age
percentiles together with the full k x k covariance of those projections
(``cov_age_ij``).  A set of points is collinear exactly when every second difference
vanishes, so nonlinearity is the joint contrast

    C = second-difference rows, e.g. [1, -2, 1] for k = 3
    W = (Cb)' (C V C')^-1 (Cb)  ~  chi2 with k - 2 df

evaluated on b = the projected effects in ascending age order.  This is a Wald test
*inside the spline's own fit*: it asks whether that fit's own trajectory is
distinguishable from a straight line.  Nothing is refitted, so it cannot drift from the
numbers the age tables carry.  For k = 3 the statistic is a single z, reported as such.

Two properties make this the right statistic rather than a weak proxy:

* It is **more powerful than comparing the two models' p-values.**  The projections share
  almost all of their uncertainty (the covariance matrix is nearly constant), and because
  every contrast row sums to zero that shared part cancels out of ``C V C'``.  Curvature
  is therefore estimated far more precisely than the level of the effect is -- a null
  result here is informative, not merely underpowered.
* It is **orientation- and scale-free** in the way that matters: it asks only whether the
  effects are collinear, so it does not care about the sign of the aging effect or the
  units of the eigenswitch.

The full covariance is **required**, not optional: dropping the off-diagonal terms would
discard exactly the shared uncertainty the contrast is designed to cancel, inflating the
variance and manufacturing a null.  A fit whose age table carries only per-point standard
errors is therefore refused rather than tested on the diagonal.  (This is why the WGCNA
baseline cannot be run here -- its ``age_spline.parquet`` has no ``cov_age_ij`` columns.
The decision this test informs is about which IsoGraph age model to quote, so that is not
a gap in the argument.)

Caveat, stated in the output: the percentiles are equally spaced in *quantile*, not in
years, so the contrast measures deviation-from-linearity across the age quantiles rather
than a second derivative in years.  For a monotone age variable the two vanish together,
which is all the decision needs.

Usage
-----
    python -m isograph_benchmark.real_data.age_model_curvature
    python -m isograph_benchmark.real_data.age_model_curvature --method wgcna
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from isograph_benchmark.paths import ensure_dir, rel, stage_out
from isograph_benchmark.real_data.module_trust import METHOD_DIRS, PROD_ROOTS


def _out_dir():
    return ensure_dir(stage_out("trust.stability", "module_trust"))


def second_difference_contrasts(k: int) -> np.ndarray:
    """(k-2) x k matrix of second differences; every row sums to zero.

    ``C @ b == 0`` for all b exactly when the k points are collinear, so these rows span
    the departures-from-a-straight-line that the projections can express.
    """
    if k < 3:
        raise ValueError(f"need at least 3 age projections to test curvature, got {k}")
    c = np.zeros((k - 2, k))
    for i in range(k - 2):
        c[i, i:i + 3] = (1.0, -2.0, 1.0)
    return c


def _cov_matrix(g: pd.DataFrame, k: int) -> np.ndarray | None:
    """The k x k covariance of the projected effects, or None if the fit did not write it."""
    if not all(f"cov_age_{i}{j}" in g.columns for i in range(1, k + 1) for j in range(1, k + 1)):
        return None
    return np.array(
        [[float(g[f"cov_age_{i}{j}"].iloc[0]) for j in range(1, k + 1)]
         for i in range(1, k + 1)]
    )


def _curvature_rows(cohort: str, region: str, method: str) -> list[dict]:
    """One row per module: curvature contrast, its Wald z, and the two models' p-values."""
    root = PROD_ROOTS[(cohort, region)]
    spline_path = rel(*root, METHOD_DIRS[method], "age_spline.parquet")
    linear_path = rel(*root, METHOD_DIRS[method], "age_linear.parquet")
    if not spline_path.exists() or not linear_path.exists():
        raise SystemExit(
            f"age tables missing for {cohort}/{region}/{method}: "
            f"{spline_path if not spline_path.exists() else linear_path}"
        )
    spl = pd.read_parquet(spline_path)
    lin = pd.read_parquet(linear_path)
    lin["module_id"] = lin["module_id"].astype(str)
    lin_p = lin.set_index("module_id")["pvalue"]

    rows: list[dict] = []
    for mid, g in spl.groupby("module_id"):
        g = g.sort_values("age_prob")
        k = len(g)
        b = g["effect"].to_numpy(dtype=float)
        cov = _cov_matrix(g, k)
        if cov is None:
            # Testing on the diagonal alone would inflate the variance by exactly the
            # shared term the contrast cancels, and so would fabricate a null result.
            raise SystemExit(
                f"{cohort}/{region}/{method}: age_spline.parquet carries no cov_age_ij "
                f"columns, so the curvature contrast has no valid variance. Refusing to "
                f"substitute the diagonal (it would inflate the variance and manufacture "
                f"a null). This test is only available for fits that write the full "
                f"projection covariance."
            )
        c = second_difference_contrasts(k)
        d = c @ b
        vd = c @ cov @ c.T
        try:
            stat = float(d @ np.linalg.solve(vd, d))
        except np.linalg.LinAlgError:
            stat = np.nan
        df_curv = k - 2
        if not np.isfinite(stat) or stat < 0:
            stat = p = z = np.nan
        else:
            p = float(stats.chi2.sf(stat, df_curv))
            # for a single contrast the chi2 is z^2; keep the signed z, which is readable
            z = float(np.sign(d[0]) * np.sqrt(stat)) if df_curv == 1 else np.nan
        mid_s = str(mid)
        rows.append({
            "cohort": cohort,
            "region": region,
            "method": method,
            "module_id": mid_s,
            "n_age_points": k,
            "effect_first": b[0],
            "effect_last": b[-1],
            "curvature": float(d[0]) if df_curv >= 1 else np.nan,
            "chi2": stat,
            "df": df_curv,
            "z": z,
            "pvalue": p,
            "p_spline_ftest": float(g["pvalue_ftest"].iloc[0]),
            "p_linear": float(lin_p.get(mid_s, np.nan)),
        })
    return rows


def run(method: str = "isograph", fdr: float = 0.05) -> pd.DataFrame:
    rows: list[dict] = []
    for cohort, region in PROD_ROOTS:
        rows.extend(_curvature_rows(cohort, region, method))
    df = pd.DataFrame(rows)
    ok = df["pvalue"].notna()
    df["qvalue"] = np.nan
    if ok.any():
        df.loc[ok, "qvalue"] = multipletests(df.loc[ok, "pvalue"], method="fdr_bh")[1]

    n = int(ok.sum())
    n_nom = int((df["pvalue"] < fdr).sum())
    n_fdr = int((df["qvalue"] < fdr).sum())
    # modules the spline calls significant that the linear model does not -- the other
    # way the spline could earn its keep even with no measurable curvature
    spline_only = int(((df["p_spline_ftest"] < 0.05) & (df["p_linear"] >= 0.05)).sum())
    lin_only = int(((df["p_linear"] < 0.05) & (df["p_spline_ftest"] >= 0.05)).sum())

    out_dir = _out_dir()
    path = out_dir / f"age_model_curvature__{method}.parquet"
    df.to_parquet(path, index=False, compression="zstd")

    summary = {
        "method": method,
        "n_modules": n,
        "fdr": fdr,
        "n_curved_nominal": n_nom,
        "n_curved_nominal_expected_by_chance": round(fdr * n, 1),
        "n_curved_fdr": n_fdr,
        "min_pvalue": float(df["pvalue"].min()) if n else None,
        "median_chi2": float(df["chi2"].median()) if n else None,
        "spline_only_significant": spline_only,
        "linear_only_significant": lin_only,
        "verdict": (
            "no detectable curvature; the spline's extra df are a power tax"
            if n_fdr == 0 and n_nom <= fdr * n
            else "curvature present in at least some modules; the spline is not redundant"
        ),
    }
    (out_dir / f"age_model_curvature__{method}.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print(f"modules tested: {n}", flush=True)
    print(
        f"curvature nominal p<{fdr}: {n_nom} "
        f"(expected by chance {fdr * n:.1f}); BH q<{fdr}: {n_fdr}",
        flush=True,
    )
    if n:
        print(
            f"smallest curvature p: {df['pvalue'].min():.3g}; "
            f"median chi2: {df['chi2'].median():.3f} on df={int(df['df'].max())}",
            flush=True,
        )
    print(
        f"spline-only significant modules: {spline_only}; "
        f"linear-only: {lin_only}",
        flush=True,
    )
    per_fit = (
        df.assign(curved=df["pvalue"] < fdr)
        .groupby(["cohort", "region"])
        .agg(n=("module_id", "size"), curved=("curved", "sum"),
             median_chi2=("chi2", lambda s: round(float(s.median()), 3)))
    )
    print("\n" + per_fit.to_string(), flush=True)
    print(f"\nverdict: {summary['verdict']}", flush=True)
    print(f"wrote {path}", flush=True)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--method", default="isograph", choices=list(METHOD_DIRS))
    ap.add_argument("--fdr", type=float, default=0.05)
    args = ap.parse_args()
    run(method=args.method, fdr=args.fdr)


if __name__ == "__main__":
    main()
