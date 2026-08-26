"""Inverse-variance meta-analysis with DerSimonian-Laird random effects.

Extracted verbatim from ``real_data/qtl_anchoring_meta.py`` so more than one analysis can
pool per-analysis log effects without re-deriving the arithmetic.  Fisher's method combines
*p*-values only: it grows more significant as analyses are added regardless of how small the
effects are, which is exactly the failure mode that makes a 1.04-fold enrichment look
impressive.  Pooling on the effect scale with a heterogeneity estimate does not.

Callers supply per-analysis ``beta`` (a log effect) and ``se``; the returned point estimates
are exponentiated, so they read as ratios / odds ratios on the original scale.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm

Z = 1.959963984540054  # two-sided 95%


def meta_keys(effect_name: str = "or", count_col: str | None = "n_foreground") -> tuple[str, ...]:
    """Column names ``meta`` produces for a given effect label."""
    keys = ["k"]
    if count_col:
        keys.append("n_fg_total")
    keys += [
        f"{effect_name}_fe", f"{effect_name}_fe_low", f"{effect_name}_fe_high", "p_fe",
        f"{effect_name}_re", f"{effect_name}_re_low", f"{effect_name}_re_high", "p_re",
        "Q", "I2",
    ]
    return tuple(keys)


def meta(
    group: pd.DataFrame,
    effect_name: str = "or",
    count_col: str | None = "n_foreground",
) -> pd.Series:
    """Fixed-effect (IVW) and DerSimonian-Laird random-effects pooling of ``beta``/``se``.

    Also returns Cochran's ``Q`` and ``I2`` so between-analysis heterogeneity is reported
    rather than hidden inside a pooled point estimate.  Designed for
    ``df.groupby(...).apply(meta, include_groups=False)``.
    """
    keys = meta_keys(effect_name, count_col)
    g = group[np.isfinite(group["beta"]) & np.isfinite(group["se"]) & (group["se"] > 0)]
    k = len(g)
    if k == 0:
        return pd.Series({key: (0 if key in ("k", "n_fg_total") else np.nan) for key in keys})

    beta, se = g["beta"].to_numpy(), g["se"].to_numpy()
    w = 1.0 / se**2
    beta_fe = float(np.sum(w * beta) / np.sum(w))
    se_fe = float(np.sqrt(1.0 / np.sum(w)))
    q = float(np.sum(w * (beta - beta_fe) ** 2))
    i2 = float(max(0.0, (q - (k - 1)) / q)) if k > 1 and q > 0 else 0.0
    tau2 = max(0.0, (q - (k - 1)) / (np.sum(w) - np.sum(w**2) / np.sum(w))) if k > 1 else 0.0
    wr = 1.0 / (se**2 + tau2)
    beta_re = float(np.sum(wr * beta) / np.sum(wr))
    se_re = float(np.sqrt(1.0 / np.sum(wr)))

    out = {"k": k}
    if count_col:
        out["n_fg_total"] = int(g[count_col].sum()) if count_col in g else 0
    out.update({
        f"{effect_name}_fe": np.exp(beta_fe),
        f"{effect_name}_fe_low": np.exp(beta_fe - Z * se_fe),
        f"{effect_name}_fe_high": np.exp(beta_fe + Z * se_fe),
        "p_fe": float(2 * norm.sf(abs(beta_fe / se_fe))),
        f"{effect_name}_re": np.exp(beta_re),
        f"{effect_name}_re_low": np.exp(beta_re - Z * se_re),
        f"{effect_name}_re_high": np.exp(beta_re + Z * se_re),
        "p_re": float(2 * norm.sf(abs(beta_re / se_re))),
        "Q": q, "I2": i2,
    })
    return pd.Series(out)
