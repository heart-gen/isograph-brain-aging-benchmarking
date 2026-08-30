"""Residual-cost ablation: does residualization cost anything when there is no confound?

The benchmark shows residualization *repairs* confounded data (``cell_composition`` ARI
0.185 -> 0.643, ``library_depth`` 0.116 -> 0.647).  It does not yet show the other half:
that regressing covariates out is free when there is nothing to remove.  Without that, the
claim has to stay "repairs the confound losses" with no at-no-cost clause, because
residualizing removes real signal whenever signal and covariate are correlated.

Why this is a separate grid
---------------------------
The obvious route -- run ``isograph_vae_residual`` over the existing unconfounded datasets
-- is closed twice over:

1. 1,730 cached datasets predate the covariate columns, so ``build_design_matrix`` receives
   nothing and residualization is a silent no-op (it would "cost nothing" because it never
   ran).
2. Of those, 1,204 -- every core switch scenario -- no longer regenerate from the current
   generator, so their sample tables cannot be repaired.  See
   ``01_synthetic_benchmark/01_synthetic/_m/SAMPLE_TABLE_REFRESH.md``.

The ablation is a **paired within-dataset contrast**: both arms see byte-identical input,
and only the residualization flag differs.  That makes commensurability with the archived
grid unnecessary -- it needs fresh data, not *the* data.  So this builds its own grid, with
``ablation`` in the dataset hash and its own dataset root, and cannot touch the archive.

Everything downstream is the production code path: ``run_one.run`` fits the models and
``run_one.compute_metrics`` scores them, so the numbers are comparable in kind to the main
benchmark even though the datasets are fresh.

Scope
-----
In unconfounded scenarios ``RIN``, ``neuron_frac`` and ``batch`` are constants that
``build_design_matrix`` drops, so this measures the cost of residualizing **library depth
only**.  Report it that way; it is not a test of residualization in general.

Usage
-----
    python -m isograph_benchmark.benchmark.residual_cost grid
    python -m isograph_benchmark.benchmark.residual_cost materialize
    sbatch 01_synthetic_benchmark/01_synthetic/_h/run_residual_cost.sh
    python -m isograph_benchmark.benchmark.residual_cost summarize
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.benchmark.run_synthetic import RESOURCE_DEFAULTS, stable_id
from isograph_benchmark.config import load_yaml
from isograph_benchmark.paths import ensure_dir, rel

CONFIG = "configs/residual_cost_ablation.yaml"
GRID_OUT = ("benchmark", "00_design", "_m", "residual_cost_grid.parquet")
PIN_OUT = ("benchmark", "00_design", "_m", "residual_cost_generator_pin.json")
DATASET_ROOT = ("benchmark", "01_synthetic", "_m", "datasets_residual_cost")
RUN_ROOT = ("benchmark", "01_synthetic", "_o", "runs_residual_cost")
OUT_ROOT = ("benchmark", "03_metrics", "_m")

BASELINE = "isograph_vae"
TREATMENT = "isograph_vae_residual"

# Metrics the ablation is judged on. A cost shows up as the residual arm scoring LOWER.
# The best-match Jaccard is `module_recovery` -- there is no `module_recovery_score` key on
# any run, old or new, so do not go looking for one.
METRICS = (
    "metrics_ari_planted",
    "metrics_ari_assigned",
    "metrics_module_recovery",
    "metrics_homogeneity_planted",
    "metrics_completeness_planted",
    "metrics_frac_planted_assigned",
)

# Every field that feeds the dataset hash. `ablation` is what guarantees these ids can
# never collide with the archived grid, even at identical parameters.
_HASH_KEYS = [
    "ablation", "scenario", "n_genes", "n_samples", "switching_fraction", "noise_sd",
    "abundance_imbalance", "count_dispersion", "interaction_strength",
    "interaction_fraction", "seed",
]


# --------------------------------------------------------------------------- #
# Grid
# --------------------------------------------------------------------------- #
def build_grid() -> pd.DataFrame:
    cfg = load_yaml(CONFIG)
    ablation = str(cfg["ablation"])
    rows: list[dict] = []
    for scenario, params in cfg["scenarios"].items():
        for seed_idx in range(int(cfg["seed_count"])):
            dataset = {
                **params,
                "ablation": ablation,
                "scenario": scenario,
                "seed": int(cfg["base_seed"]) + seed_idx,
                "replicate": seed_idx,
            }
            dataset_id = stable_id(dataset, _HASH_KEYS)
            for method in cfg["methods"]:
                row = dict(dataset)
                row.update({
                    "method": method,
                    "dataset_id": dataset_id,
                    "run_id": stable_id({**dataset, "method": method},
                                        _HASH_KEYS + ["method"]),
                    "resource_class": "vae",
                    **RESOURCE_DEFAULTS["vae"],
                })
                rows.append(row)
    grid = pd.DataFrame(rows)
    first = ["run_id", "dataset_id", "scenario", "method", "resource_class", "seed",
             "replicate", "n_genes", "n_samples", "switching_fraction", "noise_sd"]
    ordered = [c for c in first if c in grid.columns]
    ordered += [c for c in grid.columns if c not in ordered]
    return (grid.loc[:, ordered]
            .sort_values(["scenario", "replicate", "dataset_id", "method"])
            .reset_index(drop=True))


def _generator_pin() -> dict:
    """Fingerprint the generator so a future RNG-stream change is detectable.

    The archived grid diverged silently because nothing recorded which generator produced
    it; regenerating today yields a different (equally valid) draw. These datasets are
    pinned at creation so the same failure is a loud mismatch next time.
    """
    src = rel("isograph_benchmark", "benchmark", "synthetic_data.py")
    pin = {
        "synthetic_data_sha256": hashlib.sha256(src.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(rel(CONFIG).read_bytes()).hexdigest(),
    }
    try:
        pin["git_commit"] = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=str(rel(".")), check=True).stdout.strip()
        pin["git_dirty"] = bool(subprocess.run(
            ["git", "status", "--porcelain", str(src)], capture_output=True, text=True,
            cwd=str(rel(".")), check=True).stdout.strip())
    except (subprocess.CalledProcessError, FileNotFoundError):
        pin["git_commit"] = None
        pin["git_dirty"] = None
    return pin


def cmd_grid() -> int:
    grid = build_grid()
    out = rel(*GRID_OUT)
    ensure_dir(out.parent)
    grid.to_parquet(out, index=False, compression="zstd")
    pin = _generator_pin()
    rel(*PIN_OUT).write_text(json.dumps(pin, indent=2))
    n_ds = grid["dataset_id"].nunique()
    print(f"Wrote {len(grid)} runs over {n_ds} datasets -> {out}")
    print(f"generator pin: commit {pin['git_commit']} dirty={pin['git_dirty']} "
          f"synthetic_data {pin['synthetic_data_sha256'][:12]}")
    if pin["git_dirty"]:
        print("WARNING: synthetic_data.py has uncommitted changes; the pin records a "
              "commit that does not describe the generator actually used.")
    return 0


# --------------------------------------------------------------------------- #
# Summarize
# --------------------------------------------------------------------------- #
def _load_runs() -> pd.DataFrame:
    grid = pd.read_parquet(rel(*GRID_OUT))
    root = rel(*RUN_ROOT)
    recs = []
    for _, row in grid.iterrows():
        done = root / str(row["run_id"]) / "done.json"
        if not done.exists():
            continue
        payload = json.loads(done.read_text())
        rec = {"run_id": str(row["run_id"]), "dataset_id": str(row["dataset_id"]),
               "scenario": str(row["scenario"]), "method": str(row["method"]),
               "seed": int(row["seed"])}
        rec.update({f"metrics_{k}": v for k, v in (payload.get("metrics") or {}).items()})
        recs.append(rec)
    return pd.DataFrame(recs)


def _paired(df: pd.DataFrame, metric: str) -> dict:
    """Paired baseline-vs-residual contrast over datasets both arms completed."""
    wide = df.pivot_table(index="dataset_id", columns="method", values=metric)
    if BASELINE not in wide.columns or TREATMENT not in wide.columns:
        return {}
    wide = wide.dropna(subset=[BASELINE, TREATMENT])
    if wide.empty:
        return {}
    base = wide[BASELINE].to_numpy(float)
    treat = wide[TREATMENT].to_numpy(float)
    diff = treat - base
    rec = {
        "metric": metric, "n_pairs": int(len(diff)),
        "mean_baseline": float(base.mean()), "mean_residual": float(treat.mean()),
        "mean_delta": float(diff.mean()),
        "median_delta": float(np.median(diff)),
        "n_worse": int((diff < 0).sum()), "n_better": int((diff > 0).sum()),
        "wilcoxon_p": None, "ci_low": None, "ci_high": None,
    }
    if np.any(diff != 0) and len(diff) >= 6:
        rec["wilcoxon_p"] = float(stats.wilcoxon(treat, base).pvalue)
    # Paired bootstrap CI on the mean delta; seed fixed so the interval is reproducible.
    rng = np.random.default_rng(13)
    boots = np.array([diff[rng.integers(0, len(diff), len(diff))].mean()
                      for _ in range(2000)])
    rec["ci_low"], rec["ci_high"] = (float(np.quantile(boots, 0.025)),
                                     float(np.quantile(boots, 0.975)))
    return rec


def cmd_summarize() -> int:
    runs = _load_runs()
    if runs.empty:
        print("no completed runs found; run the SLURM array first")
        return 1
    outdir = ensure_dir(rel(*OUT_ROOT))
    runs.to_parquet(outdir / "residual_cost_runs.parquet", index=False)

    rows = []
    for metric in METRICS:
        if metric not in runs.columns:
            continue
        rec = _paired(runs, metric)
        if rec:
            rec["scenario"] = "ALL"
            rows.append(rec)
        for scenario, sub in runs.groupby("scenario"):
            rec = _paired(sub, metric)
            if rec:
                rec["scenario"] = scenario
                rows.append(rec)
    res = pd.DataFrame(rows)
    # Five metrics are tested on the same 100 pairs, so a single raw p is not a result.
    # BH across the pooled metrics only -- the per-scenario rows are descriptive breakdowns
    # of those same tests, not additional hypotheses.
    res["wilcoxon_q"] = np.nan
    pooled_mask = (res.scenario == "ALL") & res.wilcoxon_p.notna()
    if pooled_mask.any():
        p = res.loc[pooled_mask, "wilcoxon_p"].to_numpy(float)
        order = np.argsort(p)
        ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
        q = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
        out_q = np.empty_like(q)
        out_q[order] = q
        res.loc[pooled_mask, "wilcoxon_q"] = out_q
    res.to_parquet(outdir / "residual_cost_summary.parquet", index=False)
    _write_report(runs, res, outdir)
    print(res[res.scenario == "ALL"][
        ["metric", "n_pairs", "mean_baseline", "mean_residual", "mean_delta",
         "ci_low", "ci_high", "wilcoxon_p"]].to_string(index=False))
    return 0


def _write_report(runs: pd.DataFrame, res: pd.DataFrame, outdir) -> None:
    pin_path = rel(*PIN_OUT)
    pin = json.loads(pin_path.read_text()) if pin_path.exists() else {}
    pooled = res[res.scenario == "ALL"]

    def f(v, spec="{:.4f}"):
        return spec.format(v) if isinstance(v, (int, float)) and np.isfinite(v) else "n/a"

    lines = [
        "# Residual-cost ablation: is residualization free on unconfounded data?", "",
        f"`{BASELINE}` vs `{TREATMENT}` on **{runs.dataset_id.nunique()} fresh datasets** "
        f"across {runs.scenario.nunique()} unconfounded scenarios, paired within dataset "
        "(both arms see byte-identical input; only the residualization flag differs).", "",
        "**Scope.** In unconfounded scenarios `RIN`, `neuron_frac` and `batch` are "
        "constants that `build_design_matrix` drops, so this measures the cost of "
        "residualizing **library depth only** (`library_size`). It is not a test of "
        "residualization in general and must not be reported as one.", "",
        "**Why fresh datasets.** The archived unconfounded datasets predate the covariate "
        "columns, and 1,204 of them no longer regenerate from the current generator "
        "(`01_synthetic_benchmark/01_synthetic/_m/SAMPLE_TABLE_REFRESH.md`), so their sample tables "
        "could not be repaired. A paired within-dataset contrast does not need to be "
        "commensurable with the archive, so these were generated fresh and pinned.", "",
        f"Generator pin: commit `{pin.get('git_commit')}`, "
        f"`synthetic_data.py` sha256 `{str(pin.get('synthetic_data_sha256'))[:16]}`, "
        f"dirty={pin.get('git_dirty')}.", "",
        "## Pooled across scenarios", "",
        "`delta` is residual minus baseline, so a **negative** delta is a cost. The CI is "
        "a paired bootstrap (2,000 resamples, seed 13) on the mean delta.", "",
        "| metric | n pairs | baseline | residual | delta | 95% CI | worse/better | "
        "Wilcoxon p | q (BH) |", "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in pooled.itertuples():
        lines.append(
            f"| `{r.metric}` | {r.n_pairs} | {f(r.mean_baseline)} | {f(r.mean_residual)} | "
            f"{f(r.mean_delta, '{:+.4f}')} | {f(r.ci_low, '{:+.4f}')} to "
            f"{f(r.ci_high, '{:+.4f}')} | {r.n_worse}/{r.n_better} | "
            f"{f(r.wilcoxon_p, '{:.3g}')} | {f(r.wilcoxon_q, '{:.3g}')} |")
    lines += ["", "BH is applied across the five pooled metrics, which are five tests on "
                  "the same 100 pairs. The per-scenario table below breaks those same "
                  "tests down and is descriptive; it is not separately corrected."]

    lines += ["", "## Per scenario (`ari_planted`)", "",
              "| scenario | n pairs | baseline | residual | delta | 95% CI | Wilcoxon p |",
              "|---|---|---|---|---|---|---|"]
    per = res[(res.metric == "metrics_ari_planted") & (res.scenario != "ALL")]
    for r in per.sort_values("scenario").itertuples():
        lines.append(
            f"| {r.scenario} | {r.n_pairs} | {f(r.mean_baseline)} | {f(r.mean_residual)} | "
            f"{f(r.mean_delta, '{:+.4f}')} | {f(r.ci_low, '{:+.4f}')} to "
            f"{f(r.ci_high, '{:+.4f}')} | {f(r.wilcoxon_p, '{:.3g}')} |")

    # Verdict keyed on the primary metric: a cost is a CI lying entirely below zero.
    lines += ["", "## Verdict", ""]
    primary = pooled[pooled.metric == "metrics_ari_planted"]
    if primary.empty:
        lines.append("`ari_planted` unavailable; no verdict.")
    else:
        r = primary.iloc[0]
        q = r["wilcoxon_q"]
        survives = q is not None and np.isfinite(q) and q < 0.05
        if r["ci_high"] is not None and r["ci_high"] < 0 and survives:
            lines.append(
                f"Residualization **costs** {abs(r['mean_delta']):.4f} ARI on unconfounded "
                f"data (95% CI {r['ci_low']:+.4f} to {r['ci_high']:+.4f}, BH q="
                f"{q:.3g}). The at-no-cost clause is not supportable; report "
                "residualization as a trade, quantified by this number.")
        elif r["ci_low"] is not None and r["ci_low"] > 0 and survives:
            lines.append(
                "Residualization *improves* recovery even without a designed confound, and "
                "the effect survives correction. That is not automatically an artifact: "
                "library depth varies in every scenario (`sample_depth` is a lognormal "
                "sigma=0.25 draw that scales `expected_total`), so it is a real nuisance "
                "axis even where no confound was planted, and removing it can denoise. "
                "Report the direction, but do not sell it as a benefit -- the effect is "
                "small and the ablation was designed to bound a cost, not to demonstrate a "
                "gain.")
        else:
            direction = ("slightly positive" if r["mean_delta"] > 0 else
                         "slightly negative" if r["mean_delta"] < 0 else "zero")
            lines.append(
                f"The mean ARI change is {r['mean_delta']:+.4f}, 95% CI "
                f"{r['ci_low']:+.4f} to {r['ci_high']:+.4f}, BH q={f(q, '{:.3g}')}: **no "
                f"cost detectable** after correcting across the five metrics. The point "
                f"estimate is {direction}, which is consistent with library depth being a "
                "real nuisance axis in every scenario, but nothing here survives multiple "
                "testing and it must not be reported as a benefit.", )
            lines += ["", "The supportable claim is an **equivalence bounded by the CI**: "
                          f"residualizing library depth changes ARI by at most "
                          f"{max(abs(r['ci_low']), abs(r['ci_high'])):.4f} on unconfounded "
                          "data. That licenses the at-no-cost clause alongside the "
                          "confound-repair result -- stated as a bound, not as a proof of "
                          "exactly zero effect."]
    (outdir / "RESIDUAL_COST.md").write_text("\n".join(lines) + "\n")


def cmd_materialize() -> int:
    """Build every dataset serially, before any array job starts.

    ``save_dataset_bundle`` is not atomic and ``ensure_dataset`` gates only on
    ``manifest.json``, so two array tasks that need the same dataset -- which is exactly
    what a paired ablation does, since both arms share a dataset_id -- can write the same
    .npz concurrently and a third reader hits "No data left in file". Materialising in one
    process removes the race rather than retrying through it.
    """
    from isograph.io.artifacts import load_dataset_bundle
    from isograph_benchmark.benchmark.synthetic_data import ensure_dataset

    grid = pd.read_parquet(rel(*GRID_OUT))
    root = ensure_dir(rel(*DATASET_ROOT))
    ids = grid.drop_duplicates(subset=["dataset_id"], keep="first")
    n_built = n_ok = 0
    for row in ids.itertuples(index=False):
        row = pd.Series(row._asdict())
        d = root / str(row["dataset_id"])
        if (d / "manifest.json").exists():
            try:                       # a partial dir from an earlier race must not stand
                load_dataset_bundle(d)
                n_ok += 1
                continue
            except Exception:
                shutil.rmtree(d, ignore_errors=True)
        ensure_dataset(row, root)
        load_dataset_bundle(d)         # fail loudly here, not inside a 200-task array
        n_built += 1
    print(f"{len(ids)} datasets: {n_ok} already valid, {n_built} built")
    return 0


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("command", choices=["grid", "materialize", "summarize"])
    args = p.parse_args()
    raise SystemExit({"grid": cmd_grid, "materialize": cmd_materialize,
                      "summarize": cmd_summarize}[args.command]())


if __name__ == "__main__":
    main()
