"""Cross-region rollup of the switch coding-consequence enrichment.

Aggregates each region's `switch_consequence/consequence_enrichment.parquet` into a
per-(analysis class, stratum, consequence) summary.

The headline is a **random-effects pooled enrichment ratio with a confidence interval**,
not a combined p-value.  Fisher's method grows arbitrarily significant as analyses are
added, however small the effects are, so a 1.04-fold CDS enrichment measured in ten
analyses acquires an extremely small combined p and reads as though it were large.  Pooling
the per-analysis log enrichments by inverse variance — with Cochran's Q and I² reported —
keeps the effect size in view and shows how consistent it actually is.  The Fisher p is
retained as a secondary column so the previous summary remains reconstructible.

Aging and disease analyses are pooled separately (`analysis_class`): the rollup previously
mixed the nine aging regions with the SCZD analysis into single rows, which pools different
contrasts.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues

from isograph_benchmark.paths import cohort_dir, stage_out
from isograph_benchmark.stats.meta_analysis import meta

_ROOTS = [cohort_dir("brainseq"), cohort_dir("gtex")]

# Regions whose contrast is case/control rather than age.  Everything else is an aging arm.
_DISEASE_REGIONS = {"brainseq_caudate_sczd", "caudate_sczd"}


def _analysis_class(region: str) -> str:
    return "disease" if str(region) in _DISEASE_REGIONS else "aging"


def _collect() -> pd.DataFrame:
    frames = []
    for root in _ROOTS:
        if not root.exists():
            continue
        for region_dir in sorted(p for p in root.iterdir() if p.is_dir()):
            f = region_dir / "_m" / "isograph_vae" / "switch_consequence" / "consequence_enrichment.parquet"
            if f.exists():
                frames.append(pd.read_parquet(f))
    if not frames:
        raise SystemExit("no per-region consequence_enrichment.parquet found; run "
                         "switch_consequence across regions first.")
    df = pd.concat(frames, ignore_index=True)
    df["analysis_class"] = df["region"].map(_analysis_class)
    return df


def run(pooled_all: bool = True) -> pd.DataFrame:
    df = _collect()
    has_se = "log_enrichment_se" in df.columns

    rows = []
    groups = ["analysis_class", "stratum", "consequence"]
    for keys, sub in df.groupby(groups):
        sub = sub.dropna(subset=["enrichment", "p_emp"])
        if sub.empty:
            continue
        # empirical p floored away from 0 (perm resolution) before Fisher combine
        p = np.clip(sub["p_emp"].to_numpy(), 1e-6, 1.0)
        _, p_comb = combine_pvalues(p, method="fisher")
        row = dict(zip(groups, keys))
        row.update({
            "n_regions": len(sub),
            "n_enriched_p05": int(((sub["enrichment"] > 1) & (sub["p_emp"] < 0.05)).sum()),
            "n_depleted_p05": int(((sub["enrichment"] < 1) & (sub["p_emp"] < 0.05)).sum()),
            "median_enrichment": float(sub["enrichment"].median()),
            "median_obs_rate": float(sub["obs_rate"].median()),
            "fisher_p": float(p_comb),
        })
        if has_se:
            g = sub.rename(columns={"log_enrichment": "beta", "log_enrichment_se": "se"})
            row.update(meta(g, effect_name="ratio", count_col=None).to_dict())
        rows.append(row)

    meta_df = pd.DataFrame(rows).sort_values(groups).reset_index(drop=True)
    out = stage_out("mechanism", "switch_consequence_meta.parquet")
    out.parent.mkdir(parents=True, exist_ok=True)
    meta_df.to_parquet(out, index=False)
    _write_report(meta_df, has_se)
    n_by_class = df.groupby("analysis_class")["region"].nunique().to_dict()
    print(f"switch-consequence meta over {n_by_class} regions -> {out}")
    if not has_se:
        print("  NOTE: no log_enrichment_se column found — re-run switch_consequence to "
              "produce the gene-level block-bootstrap SEs, or only Fisher is reported.")
    return meta_df


def _write_report(meta_df: pd.DataFrame, has_se: bool) -> None:
    lines = [
        "# Switch coding-consequence enrichment — cross-region rollup", "",
        "Per consequence class and GO-stratum, pooled separately for the aging and disease "
        "analyses. The primary estimate is the **random-effects pooled enrichment ratio** "
        "(DerSimonian-Laird) over the per-analysis log enrichments, with a 95% CI and I² for "
        "between-analysis heterogeneity. The Fisher-combined p is secondary and should not be "
        "read as an effect size: it grows more significant with every added analysis "
        "regardless of magnitude. `coding_consequence` = CDS change or coding-status change.",
        "",
    ]
    if has_se:
        lines += [
            "| class | stratum | consequence | k | pooled ratio (RE) | 95% CI | p (RE) | I² | median obs rate | Fisher p |",
            "|-------|---------|-------------|---|-------------------|--------|--------|----|-----------------|----------|",
        ]
        for r in meta_df.itertuples():
            ci = (f"{r.ratio_re_low:.2f}–{r.ratio_re_high:.2f}"
                  if np.isfinite(getattr(r, "ratio_re_low", np.nan)) else "n/a")
            ratio = (f"{r.ratio_re:.3f}" if np.isfinite(getattr(r, "ratio_re", np.nan)) else "n/a")
            p_re = (f"{r.p_re:.2e}" if np.isfinite(getattr(r, "p_re", np.nan)) else "n/a")
            i2 = (f"{r.I2:.2f}" if np.isfinite(getattr(r, "I2", np.nan)) else "n/a")
            lines.append(
                f"| {r.analysis_class} | {r.stratum} | {r.consequence} | {int(r.k)} | {ratio} | "
                f"{ci} | {p_re} | {i2} | {r.median_obs_rate:.3f} | {r.fisher_p:.2e} |")
    else:
        lines += [
            "| class | stratum | consequence | regions | enriched (p<.05) | depleted | median enrich | median obs rate | Fisher p |",
            "|-------|---------|-------------|---------|------------------|----------|---------------|-----------------|----------|",
        ]
        for r in meta_df.itertuples():
            lines.append(
                f"| {r.analysis_class} | {r.stratum} | {r.consequence} | {r.n_regions} | "
                f"{r.n_enriched_p05} | {r.n_depleted_p05} | {r.median_enrichment:.2f} | "
                f"{r.median_obs_rate:.3f} | {r.fisher_p:.2e} |")
    (stage_out("mechanism", "SWITCH_CONSEQUENCE_META.md")).write_text("\n".join(lines) + "\n")


def main() -> None:
    argparse.ArgumentParser(description=__doc__).parse_args()
    run()


if __name__ == "__main__":
    main()
