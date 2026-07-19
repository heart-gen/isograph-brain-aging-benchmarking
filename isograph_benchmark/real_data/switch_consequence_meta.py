"""Cross-region rollup of the switch coding-consequence enrichment.

Aggregates each region's `switch_consequence/consequence_enrichment.parquet` into a
per-(stratum, consequence) summary: how many regions show enrichment > 1 at empirical
p < 0.05, the median enrichment, and a Fisher-combined p across regions. Answers whether
the within-gene coding/UTR-consequence signal is consistent across the switch layer.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues

from isograph_benchmark.paths import rel

_ROOTS = [rel("real_data", "brainseq"), rel("real_data", "gtex")]


def _collect() -> pd.DataFrame:
    frames = []
    for root in _ROOTS:
        for f in root.glob("*/_m/isograph_vae/switch_consequence/consequence_enrichment.parquet"):
            frames.append(pd.read_parquet(f))
    if not frames:
        raise SystemExit("no per-region consequence_enrichment.parquet found; run "
                         "switch_consequence across regions first.")
    return pd.concat(frames, ignore_index=True)


def run() -> pd.DataFrame:
    df = _collect()
    rows = []
    for (stratum, cons), sub in df.groupby(["stratum", "consequence"]):
        sub = sub.dropna(subset=["enrichment", "p_emp"])
        if sub.empty:
            continue
        # empirical p floored away from 0 (perm resolution) before Fisher combine
        p = np.clip(sub["p_emp"].to_numpy(), 1e-6, 1.0)
        _, p_comb = combine_pvalues(p, method="fisher")
        rows.append({
            "stratum": stratum, "consequence": cons,
            "n_regions": len(sub),
            "n_enriched_p05": int(((sub["enrichment"] > 1) & (sub["p_emp"] < 0.05)).sum()),
            "n_depleted_p05": int(((sub["enrichment"] < 1) & (sub["p_emp"] < 0.05)).sum()),
            "median_enrichment": float(sub["enrichment"].median()),
            "median_obs_rate": float(sub["obs_rate"].median()),
            "fisher_p": float(p_comb),
        })
    meta = pd.DataFrame(rows).sort_values(["stratum", "consequence"])
    out = rel("real_data", "_m", "switch_consequence_meta.parquet")
    out.parent.mkdir(parents=True, exist_ok=True)
    meta.to_parquet(out, index=False)
    _write_report(meta)
    print(f"switch-consequence meta over {df['region'].nunique()} regions -> {out}")
    return meta


def _write_report(meta: pd.DataFrame) -> None:
    lines = [
        "# Switch coding-consequence enrichment — cross-region rollup", "",
        "Per consequence class and GO-stratum, aggregated over the switch-layer regions: "
        "how many regions show enrichment > 1 at empirical p < 0.05 (vs the within-gene "
        "random-pair null), how many show depletion, the median enrichment, and a "
        "Fisher-combined p. `coding_consequence` = CDS change or coding-status change.", "",
        "| stratum | consequence | regions | enriched (p<.05) | depleted | median enrich | median obs rate | Fisher p |",
        "|---------|-------------|---------|------------------|----------|---------------|-----------------|----------|",
    ]
    for r in meta.itertuples():
        lines.append(
            f"| {r.stratum} | {r.consequence} | {r.n_regions} | {r.n_enriched_p05} | "
            f"{r.n_depleted_p05} | {r.median_enrichment:.2f} | {r.median_obs_rate:.3f} | "
            f"{r.fisher_p:.2e} |")
    (rel("real_data", "_m", "SWITCH_CONSEQUENCE_META.md")).write_text("\n".join(lines) + "\n")


def main() -> None:
    argparse.ArgumentParser(description="Cross-region switch-consequence rollup.").parse_args()
    run()


if __name__ == "__main__":
    main()
