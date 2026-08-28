"""Per-sweep-point breakdown of the synthetic benchmark metrics.

`summarize.py` collapses each (scenario, method, metric) to a single mean/median over
every parameter point in that scenario's sweep.  For scenarios whose sweep spans a
difficulty gradient that summary is misleading in both directions: the median is set by
how many easy points the sweep happens to contain, and the mean blends regimes where the
method is exact with regimes where it has broken down.  `feature_space_interactions` is
the clearest case — IsoGraph's ARI runs from 1.000 at `interaction_strength <= 0.5` to
0.020 at `interaction_strength >= 2` with `interaction_fraction = 0.75`, and the
scenario-level median of 0.973 shows neither number.

This module re-groups the long table by the sweep parameters that actually vary *within*
each scenario, so the dose-response is reported instead of averaged away.  It changes no
published value: `summarize.py` remains the source of the scenario-level table and of the
paired tests, which are within-simulation contrasts and therefore unaffected by sweep
composition.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.stats.hypothesis_tests import _paired_diff_ci, _paired_wilcoxon
from isograph_benchmark.stats.summarize import (
    _load_stats_config,
    benjamini_hochberg,
    bootstrap_ci,
)

# Run-level identifiers that index replicates rather than sweep position.  Grouping on
# these would return one run per cell and defeat the purpose.
_NOT_SWEEP = frozenset({"run_run_id", "run_replicate", "run_seed"})

DEFAULT_METRICS = ("metrics_ari_planted", "metrics_module_recovery")


def sweep_axes(frame: pd.DataFrame) -> list[str]:
    """Sweep parameters that take more than one value in `frame`.

    Determined per scenario rather than globally: every scenario carries the full set of
    `run_*` columns, but holds all but its own knobs fixed (usually at NaN).
    """
    cols = [c for c in frame.columns if c.startswith("run_") and c not in _NOT_SWEEP]
    return [c for c in cols if frame[c].nunique(dropna=False) > 1]


def breakdown(
    long: pd.DataFrame,
    metrics: tuple[str, ...],
    n_iter: int,
    alpha: float,
    seed: int = 0,
) -> pd.DataFrame:
    """One row per (scenario, method, metric, sweep point)."""
    rows = []
    for (scenario, metric), scen in long[long["metric"].isin(metrics)].groupby(
        ["scenario", "metric"], dropna=False
    ):
        axes = sweep_axes(scen)
        # A scenario with no varying knob has a single sweep point; report it as such so
        # the table is complete rather than silently dropping the scenario.
        group_keys = ["method"] + axes
        for keys, cell in scen.groupby(group_keys, dropna=False):
            keys = keys if isinstance(keys, tuple) else (keys,)
            vals = cell["value"].to_numpy(dtype=float)
            lo, hi = bootstrap_ci(vals, n_iter=n_iter, alpha=alpha, seed=seed)
            row = {
                "scenario": scenario,
                "metric": metric,
                "method": keys[0],
                "n": int(np.isfinite(vals).sum()),
                "mean": float(np.nanmean(vals)) if np.isfinite(vals).any() else np.nan,
                "median": float(np.nanmedian(vals)) if np.isfinite(vals).any() else np.nan,
                "ci_low": lo,
                "ci_high": hi,
                "sweep_axes": "|".join(axes),
            }
            for axis, value in zip(axes, keys[1:]):
                row[axis] = value
            rows.append(row)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    lead = ["scenario", "metric", "method", "sweep_axes", "n", "mean", "median",
            "ci_low", "ci_high"]
    return out[lead + [c for c in out.columns if c not in lead]]


def paired_by_sweep_point(
    long: pd.DataFrame,
    metric: str,
    method: str,
    method_ref: str,
    n_iter: int,
    alpha: float,
    seed: int = 13,
) -> pd.DataFrame:
    """Paired method-minus-reference contrast computed *within* each sweep point.

    `hypothesis_tests.paired_tests` pairs by simulation but pools every sweep point in a
    scenario into one contrast.  That answers "does the method win on this scenario as
    sampled", which is not the same question as "does the method win at each difficulty
    the scenario spans" — a scenario-level win can be carried entirely by its easy points
    and can conceal a reversal at its hardest one.  BH is applied across the sweep points
    of all scenarios for one metric, so the per-point claims carry their own multiplicity
    correction rather than borrowing the scenario-level one.
    """
    sub = long[(long["metric"] == metric) & (long["method"].isin([method, method_ref]))]
    rows = []
    for scenario, scen in sub.groupby("scenario", dropna=False):
        axes = sweep_axes(scen)
        # Pair on the simulation identity: sweep position plus the replicate/seed that
        # index which draw it is.  Both arms must have run the same dataset to pair.
        pair_keys = axes + [c for c in ("run_replicate", "run_seed", "run_n_genes",
                                        "run_n_samples") if c in scen.columns]
        pair_keys = list(dict.fromkeys(pair_keys))
        wide = scen.pivot_table(index=pair_keys, columns="method", values="value")
        if method not in wide.columns or method_ref not in wide.columns:
            continue
        wide = wide.dropna(subset=[method, method_ref])
        if wide.empty:
            continue
        group_axes = axes if axes else None
        groups = wide.groupby(level=axes) if group_axes else [((), wide)]
        for keys, cell in groups:
            keys = keys if isinstance(keys, tuple) else (keys,)
            a = cell[method].to_numpy(dtype=float)
            b = cell[method_ref].to_numpy(dtype=float)
            lo, hi = _paired_diff_ci(a, b, n_iter=n_iter, alpha=alpha, seed=seed)
            stat, pval = _paired_wilcoxon(a, b)
            row = {
                "scenario": scenario,
                "metric": metric,
                "method": method,
                "method_ref": method_ref,
                "n_pairs": int(len(a)),
                "mean_diff": float(np.mean(a - b)),
                "diff_ci_low": lo,
                "diff_ci_high": hi,
                "statistic": stat,
                "p_value": pval,
                "sweep_axes": "|".join(axes),
            }
            for axis, value in zip(axes, keys):
                row[axis] = value
            rows.append(row)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["p_adj"] = benjamini_hochberg(out["p_value"])
    # A point is only called either way when the CI excludes zero; "no advantage" here
    # means "not demonstrated at this point", not "proven equal".
    out["direction"] = np.where(
        out["diff_ci_low"] > 0, "method_better",
        np.where(out["diff_ci_high"] < 0, "ref_better", "indistinguishable"),
    )
    lead = ["scenario", "metric", "method", "method_ref", "sweep_axes", "n_pairs",
            "mean_diff", "diff_ci_low", "diff_ci_high", "p_value", "p_adj", "direction"]
    return out[lead + [c for c in out.columns if c not in lead]]


def sweep_range(table: pd.DataFrame) -> pd.DataFrame:
    """Best and worst sweep point per (scenario, method, metric).

    The spread between them is the quantity the scenario-level median conceals, so this
    is what the report leads with.
    """
    rows = []
    for keys, group in table.groupby(["scenario", "metric", "method"], dropna=False):
        finite = group[np.isfinite(group["mean"])]
        if finite.empty:
            continue
        best = finite.loc[finite["mean"].idxmax()]
        worst = finite.loc[finite["mean"].idxmin()]
        rows.append(
            {
                "scenario": keys[0],
                "metric": keys[1],
                "method": keys[2],
                "n_sweep_points": len(finite),
                "best_mean": best["mean"],
                "worst_mean": worst["mean"],
                "spread": best["mean"] - worst["mean"],
                "sweep_axes": best["sweep_axes"],
            }
        )
    return pd.DataFrame(rows).sort_values("spread", ascending=False).reset_index(drop=True)


def _fmt_axis_value(value) -> str:
    if isinstance(value, float) and np.isfinite(value):
        return f"{value:g}"
    return str(value)


def render_markdown(
    table: pd.DataFrame,
    scenarios: tuple[str, ...],
    methods: tuple[str, ...],
    metric: str,
) -> str:
    """Markdown fragment: one table per scenario, methods as columns."""
    lines: list[str] = []
    sub = table[(table["metric"] == metric) & (table["method"].isin(methods))]
    for scenario in scenarios:
        scen = sub[sub["scenario"] == scenario]
        if scen.empty:
            continue
        axes = [a for a in scen["sweep_axes"].iloc[0].split("|") if a]
        lines.append(f"**`{scenario}`** by {', '.join(f'`{a}`' for a in axes)}:")
        lines.append("")
        header = [a.replace("run_", "") for a in axes] + [f"`{m}`" for m in methods] + ["n"]
        lines.append("| " + " | ".join(header) + " |")
        lines.append("|" + "---|" * len(header))
        wide = scen.pivot_table(index=axes, columns="method", values="mean", dropna=False)
        counts = scen.pivot_table(index=axes, columns="method", values="n", dropna=False)
        for idx, row in wide.iterrows():
            idx = idx if isinstance(idx, tuple) else (idx,)
            cells = [_fmt_axis_value(v) for v in idx]
            cells += [f"{row[m]:.3f}" if m in row and np.isfinite(row[m]) else "—"
                      for m in methods]
            n_cell = counts.loc[idx if len(idx) > 1 else idx[0]]
            cells.append(f"{int(n_cell.max())}")
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metrics", nargs="+", default=list(DEFAULT_METRICS))
    ap.add_argument(
        "--markdown-metric",
        default="metrics_ari_planted",
        help="metric rendered into the markdown fragment",
    )
    ap.add_argument(
        "--markdown-scenarios",
        nargs="+",
        default=["feature_space_interactions", "unequal_isoform_abundance",
                 "non_switching_background"],
    )
    ap.add_argument("--markdown-methods", nargs="+",
                    default=["isograph_vae", "wgcna_gene"])
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    long_path = rel("benchmark", "03_metrics", "_m", "synthetic_metric_long.parquet")
    if not long_path.exists():
        raise SystemExit(f"Missing {long_path}\nRun isograph_benchmark.stats.summarize first.")

    cfg = _load_stats_config("configs/synthetic_grid.yaml")
    n_iter = int(cfg["bootstrap_iterations"])
    alpha = 1 - float(cfg["confidence_level"])

    long = pd.read_parquet(long_path)
    table = breakdown(long, tuple(args.metrics), n_iter=n_iter, alpha=alpha, seed=args.seed)
    if table.empty:
        raise SystemExit("No rows produced — check --metrics against the long table.")

    out = rel("benchmark", "03_metrics", "_m", "synthetic_sweep_breakdown.parquet")
    ensure_dir(out.parent)
    table.to_parquet(out, index=False, compression="zstd")
    print(f"Wrote {len(table):,} sweep-point rows to {out.name}")

    ranges = sweep_range(table)
    out_range = rel("benchmark", "03_metrics", "_m", "synthetic_sweep_range.parquet")
    ranges.to_parquet(out_range, index=False, compression="zstd")
    print(f"Wrote {len(ranges):,} rows to {out_range.name}")

    paired = paired_by_sweep_point(
        long,
        metric=args.markdown_metric,
        method=args.markdown_methods[0],
        method_ref=args.markdown_methods[-1],
        n_iter=n_iter,
        alpha=alpha,
        seed=13,
    )
    if not paired.empty:
        out_paired = rel("benchmark", "03_metrics", "_m", "synthetic_sweep_paired.parquet")
        paired.to_parquet(out_paired, index=False, compression="zstd")
        n_rev = int((paired["direction"] == "ref_better").sum())
        n_ind = int((paired["direction"] == "indistinguishable").sum())
        print(f"Wrote {len(paired):,} per-sweep-point paired contrasts to {out_paired.name} "
              f"({n_rev} reference-better, {n_ind} indistinguishable)")

    md = render_markdown(
        table,
        tuple(args.markdown_scenarios),
        tuple(args.markdown_methods),
        args.markdown_metric,
    )
    out_md = rel("benchmark", "03_metrics", "_m", "SWEEP_BREAKDOWN.md")
    out_md.write_text(
        "# Per-sweep-point breakdown\n\n"
        "Generated by `python -m isograph_benchmark.stats.sweep_breakdown`. "
        f"Metric: `{args.markdown_metric}`; cell values are means over replicates.\n\n"
        + md
    )
    print(f"Wrote {out_md.name}")

    top = ranges[ranges["metric"] == args.markdown_metric].head(8)
    print("\nLargest within-scenario sweep spreads:")
    print(top.to_string(index=False))

    if not paired.empty:
        weak = paired[paired["direction"] != "method_better"].sort_values("mean_diff")
        print("\nSweep points where the advantage is not demonstrated:")
        cols = ["scenario", "sweep_axes", "n_pairs", "mean_diff", "diff_ci_low",
                "diff_ci_high", "p_adj", "direction"]
        print(weak[cols].to_string(index=False) if not weak.empty else "  (none)")


if __name__ == "__main__":
    main()
