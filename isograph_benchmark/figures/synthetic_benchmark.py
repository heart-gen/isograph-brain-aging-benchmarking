"""Publication-quality figures for the IsoGraph synthetic benchmark."""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Design constants
# ---------------------------------------------------------------------------

METHOD_ORDER = [
    "isograph_baseline",
    "isograph_latent",
    "isograph_graph",
    "isograph_cpu_latent",
    "isograph_vae",
    "wgcna_gene",
]

# Methods available in the scale scenario (subset of METHOD_ORDER)
SCALE_METHOD_ORDER = [
    "isograph_vae",
    "wgcna_gene",
]

METHOD_LABELS: dict[str, str] = {
    "isograph_baseline":   "IsoGraph\nBaseline",
    "isograph_latent":     "IsoGraph\nLatent",
    "isograph_graph":      "IsoGraph\nGraph",
    "isograph_cpu_latent": "IsoGraph\nCPU Latent",
    "isograph_vae":        "IsoGraph\nVAE",
    "wgcna_gene":          "WGCNA",
}

METHOD_LABELS_SHORT: dict[str, str] = {
    "isograph_baseline":   "Baseline",
    "isograph_latent":     "Latent",
    "isograph_graph":      "Graph",
    "isograph_cpu_latent": "CPU Latent",
    "isograph_vae":        "VAE",
    "wgcna_gene":          "WGCNA",
}

# Okabe-Ito colorblind-safe palette (7 distinct)
METHOD_COLORS: dict[str, str] = {
    "isograph_baseline":   "#0072B2",  # blue
    "isograph_latent":     "#E69F00",  # orange
    "isograph_graph":      "#009E73",  # bluish green
    "isograph_cpu_latent": "#56B4E9",  # sky blue
    "isograph_vae":        "#D55E00",  # vermillion
    "wgcna_gene":          "#CC79A7",  # reddish purple
}

SCENARIO_ORDER = [
    "idealized_switching",
    "noise_stress",
    "feature_space_interactions",
    "non_switching_background",
    "unequal_isoform_abundance",
    "scale_realistic",
]

SCENARIO_LABELS: dict[str, str] = {
    "idealized_switching":        "Idealized\nSwitching",
    "noise_stress":               "Noise\nStress",
    "feature_space_interactions": "Feature\nInteractions",
    "non_switching_background":   "Non-Switching\nBackground",
    "unequal_isoform_abundance":  "Unequal\nAbundance",
    "scale_realistic":            "BrainSEQ\nScale (16k)",
}

# Matplotlib global style
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.75,
    "xtick.major.width": 0.75,
    "ytick.major.width": 0.75,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

PANEL_LABEL_KW = dict(fontsize=10, fontweight="bold", va="top", ha="right")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_raw(path: Path, method_order: list[str] = METHOD_ORDER) -> pd.DataFrame:
    df = pd.read_parquet(path)
    df = df[df["status"] == "completed"].copy()
    present = [m for m in method_order if m in df["run_method"].unique()]
    return df[df["run_method"].isin(present)].copy()


def load_summary(path: Path, method_order: list[str] = METHOD_ORDER) -> pd.DataFrame:
    df = pd.read_parquet(path)
    present = [m for m in method_order if m in df["method"].unique()]
    return df[df["method"].isin(present)].copy()


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _active_methods(df: pd.DataFrame, col: str = "run_method") -> list[str]:
    present = set(df[col].unique()) if col in df.columns else set(df["method"].unique())
    return [m for m in METHOD_ORDER if m in present]


def _draw_dots_ci(
    ax: plt.Axes,
    summary: pd.DataFrame,
    scenario: str,
    metric: str,
    methods: list[str],
    y_label: str | None = None,
    y_lim: tuple[float, float] | None = (0, 1),
    show_x_labels: bool = True,
) -> None:
    sub = summary[(summary["scenario"] == scenario) & (summary["metric"] == metric)]
    for i, method in enumerate(methods):
        row = sub[sub["method"] == method]
        if row.empty:
            continue
        row = row.iloc[0]
        color = METHOD_COLORS[method]
        ax.errorbar(
            i, row["mean"],
            yerr=[[row["mean"] - row["ci_low"]], [row["ci_high"] - row["mean"]]],
            fmt="o", color=color, markerfacecolor=color,
            markersize=4, capsize=3, linewidth=1, capthick=1, zorder=3,
        )

    ax.set_xlim(-0.7, len(methods) - 0.3)
    if y_lim is not None:
        ax.set_ylim(*y_lim)
    if show_x_labels:
        ax.set_xticks(range(len(methods)))
        ax.set_xticklabels(
            [METHOD_LABELS_SHORT[m] for m in methods], rotation=45, ha="right",
        )
    else:
        ax.set_xticks([])
    if y_label:
        ax.set_ylabel(y_label, labelpad=4)
    ax.axhline(0, color="black", linewidth=0.5, linestyle="--", alpha=0.4)
    ax.grid(axis="y", linewidth=0.4, alpha=0.4)


def _line_ci_from_raw(
    ax: plt.Axes,
    sub_raw: pd.DataFrame,
    param_col: str,
    metric: str,
    methods: list[str],
    param_vals: list,
) -> None:
    for method in methods:
        msub = sub_raw[sub_raw["run_method"] == method]
        color = METHOD_COLORS[method]
        xs, means, lo, hi = [], [], [], []
        for pv in param_vals:
            vals = msub[msub[param_col] == pv][metric].dropna().values
            if len(vals) == 0:
                continue
            m = float(np.mean(vals))
            se = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0
            xs.append(pv)
            means.append(m)
            lo.append(m - se)
            hi.append(m + se)
        if xs:
            ax.plot(xs, means, "o-", color=color, linewidth=1.5, markersize=4,
                    label=METHOD_LABELS_SHORT[method])
            ax.fill_between(xs, lo, hi, color=color, alpha=0.15)


def _heatmap_panel(
    ax: plt.Axes,
    sub: pd.DataFrame,
    method: str,
    metric: str,
    row_col: str,
    col_col: str,
    row_vals: list,
    col_vals: list,
    vmin: float = 0.0,
    vmax: float = 1.0,
    cmap: str = "viridis",
) -> object:
    msub = sub[sub["run_method"] == method]
    pivot = (
        msub.groupby([row_col, col_col])[metric]
        .mean()
        .unstack(col_col)
        .reindex(index=row_vals, columns=col_vals)
    )
    im = ax.imshow(
        pivot.values, aspect="auto", cmap=cmap,
        vmin=vmin, vmax=vmax, origin="lower",
    )
    for i, rv in enumerate(row_vals):
        for j, cv in enumerate(col_vals):
            val = pivot.iloc[i, j] if (i < len(row_vals) and j < len(col_vals)) else np.nan
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=5.5, color="white" if val < 0.55 else "black")
    return im


def _heatmap_figure(
    raw: pd.DataFrame,
    scenario: str,
    metric: str,
    row_col_candidates: list[str],
    col_col_candidates: list[str],
    row_label: str,
    col_label: str,
    row_fmt: str = ".2f",
    col_fmt: str = ".2f",
    cbar_label: str = "Module recovery (mean)",
) -> plt.Figure | None:
    sub = raw[
        (raw["run_scenario"] == scenario) &
        raw["run_method"].isin(METHOD_ORDER) &
        raw[metric].notna()
    ].copy()
    if sub.empty:
        return None

    row_col = _find_col(sub, row_col_candidates)
    col_col = _find_col(sub, col_col_candidates)
    if not row_col or not col_col:
        return None

    methods = _active_methods(sub)
    n_methods = len(methods)
    ncols = min(n_methods, 3)
    nrows = (n_methods + ncols - 1) // ncols

    row_vals = sorted(sub[row_col].dropna().unique())
    col_vals = sorted(sub[col_col].dropna().unique())

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.8 * nrows))
    axes = np.array(axes).reshape(nrows, ncols)
    im = None

    for idx, method in enumerate(methods):
        r, c = divmod(idx, ncols)
        ax = axes[r, c]
        im = _heatmap_panel(
            ax, sub, method=method, metric=metric,
            row_col=row_col, col_col=col_col,
            row_vals=row_vals, col_vals=col_vals,
        )
        ax.set_xticks(range(len(col_vals)))
        ax.set_xticklabels(
            [format(v, col_fmt) for v in col_vals], rotation=45, ha="right",
        )
        ax.set_yticks(range(len(row_vals)))
        ax.set_yticklabels([format(v, row_fmt) for v in row_vals])
        if r == nrows - 1:
            ax.set_xlabel(col_label, labelpad=4)
        if c == 0:
            ax.set_ylabel(row_label, labelpad=4)
        ax.set_title(METHOD_LABELS_SHORT[method], pad=4)

    for idx in range(n_methods, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r, c].set_visible(False)

    if im is not None:
        cbar = fig.colorbar(im, ax=axes, orientation="vertical", shrink=0.6, pad=0.02)
        cbar.set_label(cbar_label, labelpad=4)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def make_fig1(raw: pd.DataFrame, summary: pd.DataFrame, out_dir: Path) -> None:
    scenarios = [s for s in SCENARIO_ORDER if s in summary["scenario"].unique()]
    methods = _active_methods(summary, col="method")
    n_methods = len(methods)
    n_scenarios = len(scenarios)

    # Width: ~2.8 in per scenario; height: ~3 in per metric row + 2.5 for violin
    fig = plt.figure(figsize=(2.8 * n_scenarios, 3.0 * 3 + 2.5))
    outer = gridspec.GridSpec(
        4, 1, figure=fig,
        height_ratios=[3, 3, 3, 2.5],
        hspace=0.6,
    )

    panel_metrics = [
        ("metrics_module_recovery",            "Module recovery [0–1]",            (0, 1)),
        ("metrics_switch_gene_detection_rate",  "Switch gene detection rate [0–1]", (0, 1)),
        ("metrics_nonswitch_gene_module_rate",  "Non-switching gene\nmodule rate [0–1]", (0, 1)),
    ]
    panel_labels = ["A", "B", "C", "D"]

    # Panels A, B, C
    for row_idx, (metric, ylabel, ylim) in enumerate(panel_metrics):
        inner = gridspec.GridSpecFromSubplotSpec(
            1, n_scenarios, subplot_spec=outer[row_idx], wspace=0.4,
        )
        for col_idx, scenario in enumerate(scenarios):
            ax = fig.add_subplot(inner[col_idx])
            show_x = (row_idx == len(panel_metrics) - 1)
            _draw_dots_ci(
                ax, summary, scenario=scenario, metric=metric,
                methods=methods, y_lim=ylim, show_x_labels=show_x,
            )
            if col_idx == 0:
                ax.set_ylabel(ylabel, labelpad=4)
                ax.text(-0.28, 1.10, panel_labels[row_idx],
                        transform=ax.transAxes, **PANEL_LABEL_KW)
            else:
                ax.set_ylabel("")
            if row_idx == 0:
                ax.set_title(SCENARIO_LABELS[scenario], pad=4)

    # Panel D — Runtime violin
    inner_d = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[3])
    ax_d = fig.add_subplot(inner_d[0])
    runtime_col = "measurement_elapsed_sec"

    if runtime_col in raw.columns:
        for i, method in enumerate(methods):
            vals = raw[raw["run_method"] == method][runtime_col].dropna().values
            if len(vals) == 0:
                continue
            vp = ax_d.violinplot(vals, positions=[i], showmedians=True, showextrema=False)
            color = METHOD_COLORS[method]
            for pc in vp["bodies"]:
                pc.set_facecolor(color)
                pc.set_alpha(0.7)
                pc.set_edgecolor("none")
            vp["cmedians"].set_color("white")
            vp["cmedians"].set_linewidth(1.5)

        ax_d.set_yscale("log")
        ax_d.set_xticks(range(n_methods))
        ax_d.set_xticklabels([METHOD_LABELS_SHORT[m] for m in methods],
                              rotation=30, ha="right")
        ax_d.set_ylabel("Runtime (s, log scale)", labelpad=4)
        ax_d.set_xlim(-0.7, n_methods - 0.3)
        ax_d.grid(axis="y", linewidth=0.4, alpha=0.4)
        ax_d.text(-0.06, 1.12, panel_labels[3],
                  transform=ax_d.transAxes, **PANEL_LABEL_KW)
    else:
        ax_d.text(0.5, 0.5, "Runtime data not available",
                  ha="center", va="center", transform=ax_d.transAxes)

    # Shared legend
    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=METHOD_COLORS[m], markersize=6,
                   label=METHOD_LABELS_SHORT[m])
        for m in methods
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.01, 0.99),
               frameon=False, ncol=1)

    _save(fig, out_dir, "fig1_benchmark_overview")
    plt.close(fig)
    print("  fig1_benchmark_overview saved")


# ---------------------------------------------------------------------------
# Supplementary S1–S3: Heatmaps
# ---------------------------------------------------------------------------

def make_figS1(raw: pd.DataFrame, out_dir: Path) -> None:
    fig = _heatmap_figure(
        raw, scenario="idealized_switching", metric="metrics_module_recovery",
        row_col_candidates=["run_switching_fraction", "switching_fraction"],
        col_col_candidates=["run_noise_sd", "noise_sd"],
        row_label="Switching fraction", col_label="Noise SD",
    )
    if fig is None:
        print("  figS1: no data, skipping"); return
    _save(fig, out_dir, "figS1_idealized_switching")
    plt.close(fig)
    print("  figS1_idealized_switching saved")


def make_figS2(raw: pd.DataFrame, out_dir: Path) -> None:
    fig = _heatmap_figure(
        raw, scenario="noise_stress", metric="metrics_module_recovery",
        row_col_candidates=["run_count_dispersion", "count_dispersion"],
        col_col_candidates=["run_noise_sd", "noise_sd"],
        row_label="Count dispersion", col_label="Noise SD",
        row_fmt=".0f",
    )
    if fig is None:
        print("  figS2: no data, skipping"); return
    _save(fig, out_dir, "figS2_noise_stress")
    plt.close(fig)
    print("  figS2_noise_stress saved")


def make_figS3(raw: pd.DataFrame, out_dir: Path) -> None:
    fig = _heatmap_figure(
        raw, scenario="feature_space_interactions", metric="metrics_module_recovery",
        row_col_candidates=["run_interaction_strength", "interaction_strength"],
        col_col_candidates=["run_interaction_fraction", "interaction_fraction"],
        row_label="Interaction strength", col_label="Interaction fraction",
        row_fmt=".1f",
    )
    if fig is None:
        print("  figS3: no data, skipping"); return
    _save(fig, out_dir, "figS3_feature_interactions")
    plt.close(fig)
    print("  figS3_feature_interactions saved")


# ---------------------------------------------------------------------------
# Supplementary S4: Non-switching background
# ---------------------------------------------------------------------------

def make_figS4(raw: pd.DataFrame, out_dir: Path) -> None:
    scenario = "non_switching_background"
    sub = raw[raw["run_scenario"] == scenario].copy()
    if sub.empty:
        print("  figS4: no data, skipping"); return

    sw_col = _find_col(sub, ["run_switching_fraction", "switching_fraction"])
    if not sw_col:
        print("  figS4: missing switching_fraction column, skipping"); return

    methods = _active_methods(sub)
    sw_vals = sorted(sub[sw_col].dropna().unique())

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    for ax_idx, (metric, ylabel) in enumerate([
        ("metrics_module_recovery",           "Module recovery [0–1]"),
        ("metrics_nonswitch_gene_module_rate", "Non-switching gene\nmodule rate [0–1]"),
    ]):
        ax = axes[ax_idx]
        _line_ci_from_raw(ax, sub, sw_col, metric, methods, sw_vals)
        ax.set_xlabel("Background switching fraction", labelpad=4)
        ax.set_ylabel(ylabel, labelpad=4)
        ax.set_ylim(0, 1)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
        if ax_idx == 1:
            ax.legend(frameon=False, loc="upper left")

    fig.tight_layout()
    _save(fig, out_dir, "figS4_nonswitching_background")
    plt.close(fig)
    print("  figS4_nonswitching_background saved")


# ---------------------------------------------------------------------------
# Supplementary S5: Unequal isoform abundance
# ---------------------------------------------------------------------------

def make_figS5(raw: pd.DataFrame, out_dir: Path) -> None:
    scenario = "unequal_isoform_abundance"
    sub = raw[raw["run_scenario"] == scenario].copy()
    if sub.empty:
        print("  figS5: no data, skipping"); return

    ab_col = _find_col(sub, ["run_abundance_imbalance", "abundance_imbalance"])
    if not ab_col:
        print("  figS5: missing abundance_imbalance column, skipping"); return

    methods = _active_methods(sub)
    ab_vals = sorted(sub[ab_col].dropna().unique())

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    for ax_idx, (metric, ylabel) in enumerate([
        ("metrics_module_recovery",            "Module recovery [0–1]"),
        ("metrics_switch_gene_detection_rate",  "Switch gene detection\nrate [0–1]"),
    ]):
        ax = axes[ax_idx]
        _line_ci_from_raw(ax, sub, ab_col, metric, methods, ab_vals)
        ax.set_xlabel("Abundance imbalance ratio", labelpad=4)
        ax.set_ylabel(ylabel, labelpad=4)
        ax.set_ylim(0, 1)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
        if ax_idx == 1:
            ax.legend(frameon=False, loc="upper right")

    fig.tight_layout()
    _save(fig, out_dir, "figS5_unequal_abundance")
    plt.close(fig)
    print("  figS5_unequal_abundance saved")


# ---------------------------------------------------------------------------
# Supplementary S6: Scaling (n_genes)
# ---------------------------------------------------------------------------

def make_figS6(raw: pd.DataFrame, out_dir: Path) -> None:
    # Combine scale and scale_realistic so the scaling curve extends to 16k genes.
    sub = raw[raw["run_scenario"].isin({"scale", "scale_realistic"})].copy()
    if sub.empty:
        print("  figS6: no scale data, skipping"); return

    gene_col = _find_col(sub, ["run_n_genes", "n_genes"])
    if not gene_col:
        print("  figS6: missing n_genes column, skipping"); return

    methods = [m for m in SCALE_METHOD_ORDER if m in sub["run_method"].unique()]
    if not methods:
        print("  figS6: no scale methods present, skipping"); return

    gene_vals = sorted(sub[gene_col].dropna().unique())

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    for ax_idx, (metric, ylabel, yscale) in enumerate([
        ("metrics_module_recovery",  "Module recovery [0–1]", "linear"),
        ("measurement_elapsed_sec",  "Runtime (s, log scale)",  "log"),
    ]):
        ax = axes[ax_idx]
        _line_ci_from_raw(ax, sub, gene_col, metric, methods, gene_vals)
        ax.set_xlabel("Number of genes", labelpad=4)
        ax.set_ylabel(ylabel, labelpad=4)
        ax.set_yscale(yscale)
        if yscale == "linear":
            ax.set_ylim(0, 1)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
        if ax_idx == 0:
            ax.legend(frameon=False, loc="lower left")

    fig.tight_layout()
    _save(fig, out_dir, "figS6_scale")
    plt.close(fig)
    print("  figS6_scale saved")


# ---------------------------------------------------------------------------
# Table 1
# ---------------------------------------------------------------------------

def make_table1(
    summary: pd.DataFrame,
    out_dir: Path,
    tests: pd.DataFrame | None = None,
) -> None:
    """Write Table 1.

    When *tests* (pairwise Wilcoxon results from hypothesis_tests.py) is
    provided, each metric cell is annotated with a significance marker:
      **   BH-adjusted p < 0.05
      *    BH-adjusted p < 0.10
    A blank marker means the comparison was not significant or not available.
    Direction is encoded as ↑ (method better) or ↓ (method worse vs. WGCNA).
    """
    scenarios = [s for s in SCENARIO_ORDER if s in summary["scenario"].unique()]
    methods = _active_methods(summary, col="method")

    metric_cols = {
        "metrics_module_recovery":            "Module Recovery",
        "metrics_switch_gene_detection_rate":  "Switch Detection",
        "metrics_nonswitch_gene_module_rate":  "False Positive Rate",
    }

    # Build a lookup: (scenario, metric, method) → (p_adj, direction)
    test_lookup: dict[tuple, tuple[float, str]] = {}
    if tests is not None and not tests.empty:
        for _, tr in tests.iterrows():
            test_lookup[(tr["scenario"], tr["metric"], tr["method"])] = (
                float(tr["p_adj"]) if pd.notna(tr["p_adj"]) else np.nan,
                str(tr["direction"]),
            )

    def _sig_marker(p_adj: float, direction: str) -> str:
        arrow = "↑" if direction == "method_better" else ("↓" if direction == "ref_better" else "")
        if np.isnan(p_adj):
            return ""
        if p_adj < 0.05:
            return f"**{arrow}"
        if p_adj < 0.10:
            return f"*{arrow}"
        return ""

    rows = []
    for scenario in scenarios:
        for method in methods:
            row: dict = {
                "Scenario": scenario.replace("_", " ").title(),
                "Method": METHOD_LABELS_SHORT[method],
            }
            for metric, col_label in metric_cols.items():
                sub = summary[
                    (summary["scenario"] == scenario) &
                    (summary["method"] == method) &
                    (summary["metric"] == metric)
                ]
                if sub.empty:
                    row[f"{col_label} (mean [95% CI])"] = "—"
                    row[f"{col_label} sig."] = ""
                else:
                    r = sub.iloc[0]
                    row["N"] = int(r["n"])
                    row[f"{col_label} (mean [95% CI])"] = (
                        f"{r['mean']:.3f} [{r['ci_low']:.3f}, {r['ci_high']:.3f}]"
                    )
                    p_adj, direction = test_lookup.get(
                        (scenario, metric, method), (np.nan, "")
                    )
                    row[f"{col_label} sig."] = _sig_marker(p_adj, direction)

            rt_sub = summary[
                (summary["scenario"] == scenario) &
                (summary["method"] == method) &
                (summary["metric"] == "measurement_elapsed_sec")
            ]
            row["Runtime s (median)"] = (
                f"{rt_sub.iloc[0]['median']:.1f}" if not rt_sub.empty else "—"
            )
            rows.append(row)

    table = pd.DataFrame(rows)
    out_path = out_dir / "table1_benchmark_summary.csv"
    table.to_csv(out_path, index=False)
    print(f"  table1_benchmark_summary.csv: {len(table)} rows")
    if tests is not None and not tests.empty:
        n_sig = (tests["p_adj"] < 0.05).sum()
        print(f"  Significance annotations applied ({n_sig} FDR<0.05 tests marked)")


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _save(fig: plt.Figure, out_dir: Path, name: str) -> None:
    ensure_dir(out_dir)
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / f"{name}.{ext}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    raw_path     = rel("benchmark", "01_synthetic", "_m", "synthetic_results.parquet")
    summary_path = rel("benchmark", "02_metrics",   "_m", "synthetic_metric_summary.parquet")
    tests_path   = rel("benchmark", "02_metrics",   "_m", "synthetic_pairwise_tests.parquet")
    out_dir      = rel("benchmark", "02_metrics", "figures")
    table_dir    = rel("benchmark", "02_metrics", "_m")

    if not raw_path.exists():
        raise SystemExit(f"Missing: {raw_path}\nRun step_1_collect.sh first.")
    if not summary_path.exists():
        raise SystemExit(f"Missing: {summary_path}\nRun step_2_summarize.sh first.")

    print("Loading data...")
    raw = load_raw(raw_path)
    summary = load_summary(summary_path)
    tests = pd.read_parquet(tests_path) if tests_path.exists() else None
    if tests is None:
        print("  Note: pairwise tests file not found — Table 1 will have no significance markers.")
    print(f"  Completed runs: {len(raw):,}  |  Summary rows: {len(summary):,}")

    print("Generating figures...")
    make_fig1(raw, summary, out_dir)
    make_figS1(raw, out_dir)
    make_figS2(raw, out_dir)
    make_figS3(raw, out_dir)
    make_figS4(raw, out_dir)
    make_figS5(raw, out_dir)
    make_figS6(raw, out_dir)
    make_table1(summary, table_dir, tests=tests)

    print(f"Done. Outputs in {out_dir}")


if __name__ == "__main__":
    main()
