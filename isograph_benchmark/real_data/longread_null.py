"""Matched nulls and table helper shared by the long-read confirmation arms.

Kept free of IsoGraph-package imports so ``switch_orthogonal_confirm`` (02b) and
``longread_driver_confirm`` (02f) score their focal sets against *identical* null code, and
so 02f can run on a machine with only the analysis requirements installed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Abundance strata for matching. Pairs are binned on the mean isoform fraction of their
# better-expressed member, which is what actually governs whether a usage correlation is
# estimable at n=12.
_N_STRATA = 10

# The stricter joint null (Methods, "Long-read transcript support"). Abundance alone does
# not equalize transcript detection, within-gene isoform count or usage variability, and
# each of those moves the chance rate of a negative usage correlation. This null matches
# jointly on the abundance quintile of the better-expressed member, the tertile of the
# pair's mean per-transcript detection rate, the tertile of the gene's long-read isoform
# count, and the tertile of the anchored transcript's across-sample usage SD.
JOINT_KEYS = ("max_mean_if", "detect_rate", "n_iso", "usage_sd")
JOINT_BINS = (5, 3, 3, 3)


def _md_table(frame: pd.DataFrame, floats: int = 4) -> str:
    """Pipe table without pulling in ``tabulate`` for one call."""
    def cell(v):
        if isinstance(v, float):
            return "" if pd.isna(v) else f"{v:.{floats}g}"
        return "" if v is None or (isinstance(v, float) and pd.isna(v)) else str(v)

    cols = list(frame.columns)
    head = "| " + " | ".join(cols) + " |"
    rule = "|" + "|".join("---" for _ in cols) + "|"
    body = [
        "| " + " | ".join(cell(v) for v in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([head, rule, *body])


def assign_strata(values: pd.Series, edges: np.ndarray | None = None,
                  n_strata: int = _N_STRATA):
    """Quantile strata on ``values`` (deciles by default); returns (labels, edges) so the
    background's edges can be reused for the focal set. Tied quantiles merge strata."""
    finite = values[np.isfinite(values)]
    if edges is None:
        qs = np.linspace(0, 1, n_strata + 1)
        edges = np.unique(np.nanquantile(finite, qs))
    labels = np.digitize(values.to_numpy(dtype=float), edges[1:-1], right=True)
    labels = np.where(np.isfinite(values.to_numpy(dtype=float)), labels, -1)
    return labels, edges


def matched_null(
    focal: pd.DataFrame,
    background: pd.DataFrame,
    statistic: str,
    n_draws: int,
    seed: int,
    keys: tuple[str, ...] = ("max_mean_if",),
    bins: tuple[int, ...] | None = None,
) -> dict:
    """Draw stratum-matched background sets and compare ``statistic`` against the focal set.

    Matching is on the joint stratum of ``keys``, each cut into ``bins`` quantile strata on
    the background (default: abundance deciles alone), with replacement *within* stratum.
    Pairs with a non-finite key are excluded from both sides. A focal pair whose stratum is
    unrepresented in the background is dropped and counted, so the comparison is never made
    across a stratum gap. ``JOINT_KEYS``/``JOINT_BINS`` give the stricter joint null.
    """
    keys = tuple(keys)
    bins = tuple(bins) if bins is not None else (_N_STRATA,) * len(keys)
    if len(bins) != len(keys):
        raise ValueError("one bin count per matching key")
    finite_f = np.logical_and.reduce([np.isfinite(focal[k].to_numpy(float)) for k in keys])
    finite_b = np.logical_and.reduce([np.isfinite(background[k].to_numpy(float)) for k in keys])
    focal = focal[finite_f].copy()
    bg = background[finite_b].copy()
    # Joint stratum as one integer, key by key in mixed radix. With a single key the label
    # is that key's stratum itself, so the default draws are those of the abundance null.
    bg["stratum"] = 0
    focal["stratum"] = 0
    radix = 1
    for key, n in zip(keys, bins):
        lab_b, edges = assign_strata(bg[key], n_strata=n)
        lab_f, _ = assign_strata(focal[key], edges)
        bg["stratum"] += lab_b * radix
        focal["stratum"] += lab_f * radix
        radix *= n

    pools = {s: g.index.to_numpy() for s, g in bg.groupby("stratum")}
    usable = focal[focal["stratum"].isin(pools)].copy()
    dropped = int(len(focal) - len(usable))
    if usable.empty:
        return {
            "statistic": statistic,
            "n_focal": 0,
            "n_focal_dropped_no_stratum": dropped,
            "note": "no focal pair shares a matching stratum with the background",
        }

    def stat_of(frame: pd.DataFrame) -> float:
        if statistic == "switch_like_rate":
            return float(frame["switch_like"].mean())
        if statistic == "mean_usage_spearman":
            return float(frame["usage_spearman"].mean(skipna=True))
        raise ValueError(f"unknown statistic {statistic}")

    obs = stat_of(usable)
    rng = np.random.default_rng(seed)
    want = usable["stratum"].to_numpy()
    draws = np.empty(n_draws, dtype=float)
    for b in range(n_draws):
        pick = [rng.choice(pools[s]) for s in want]
        draws[b] = stat_of(bg.loc[pick])

    finite = draws[np.isfinite(draws)]
    # Two-sided empirical p with the +1 correction (a p of exactly 0 is not defensible).
    centred = np.abs(finite - np.nanmean(finite))
    p_emp = float((np.sum(centred >= abs(obs - np.nanmean(finite))) + 1) / (len(finite) + 1))
    # The matching keys are recorded only for a non-default null, so abundance-only
    # summaries keep their committed shape.
    matching = {} if keys == ("max_mean_if",) and bins == (_N_STRATA,) else {
        "match_keys": list(keys),
        "match_bins": list(bins),
        "n_strata_focal": int(usable["stratum"].nunique()),
    }
    return {
        "statistic": statistic,
        "n_focal": int(len(usable)),
        "n_focal_dropped_no_stratum": dropped,
        "n_background": int(len(bg)),
        "observed": obs,
        "null_mean": float(np.nanmean(finite)),
        "null_sd": float(np.nanstd(finite)),
        "null_q025": float(np.nanquantile(finite, 0.025)),
        "null_q975": float(np.nanquantile(finite, 0.975)),
        "p_empirical_two_sided": p_emp,
        "n_draws": int(len(finite)),
        **matching,
    }
