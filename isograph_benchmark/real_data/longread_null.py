"""Abundance-matched null and table helper shared by the long-read confirmation arms.

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


def assign_strata(values: pd.Series, edges: np.ndarray | None = None):
    """Decile strata on ``values``; returns (labels, edges) so the background reuses them."""
    finite = values[np.isfinite(values)]
    if edges is None:
        qs = np.linspace(0, 1, _N_STRATA + 1)
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
) -> dict:
    """Draw stratum-matched background sets and compare ``statistic`` against the focal set.

    Matching is on abundance stratum with replacement *within* stratum. A focal pair whose
    stratum is unrepresented in the background is dropped from both sides and counted, so
    the comparison is never made across a stratum gap.
    """
    focal = focal[np.isfinite(focal["max_mean_if"])].copy()
    bg = background[np.isfinite(background["max_mean_if"])].copy()
    strata, edges = assign_strata(bg["max_mean_if"])
    bg["stratum"] = strata
    focal["stratum"], _ = assign_strata(focal["max_mean_if"], edges)

    pools = {s: g.index.to_numpy() for s, g in bg.groupby("stratum")}
    usable = focal[focal["stratum"].isin(pools)].copy()
    dropped = int(len(focal) - len(usable))
    if usable.empty:
        return {
            "statistic": statistic,
            "n_focal": 0,
            "n_focal_dropped_no_stratum": dropped,
            "note": "no focal pair shares an abundance stratum with the background",
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
    }
