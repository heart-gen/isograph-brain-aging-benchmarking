"""Cross-tissue meta-analysis of the GTEx xQTL anchoring results.

Each per-analysis run (qtl_anchoring.py) produces a matched module-membership odds
ratio for sQTL and eQTL in one brain region/cohort, on that region's own co-switch
modules. The hypothesis — co-switch module genes are enriched for sQTL (splicing)
and not for eQTL (expression) — is shared across regions, so the per-region effects
are pooled with an inverse-variance meta-analysis (fixed effect + DerSimonian-Laird
random effects) per (xqtl_kind, module_set). This recovers the power a single
underpowered tissue lacks.

The per-analysis log-OR standard error is recovered from the saved 95% CI (matched
logistic rows) or from OR + p-value (unmatched Fisher fallback rows).

Writes under real_data/_m/qtl_anchoring_meta/:
  qtl_anchoring_meta.parquet — pooled OR/CI/p, heterogeneity (Q, I2), k per cell.
  per_analysis.parquet       — the collected per-analysis rows with derived se.
  QTL_ANCHORING_META.md      — the sQTL-vs-eQTL pooled contrast writeup.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.sweep_leiden import _artifact_dir
from isograph_benchmark.stats.meta_analysis import meta, meta_keys

# (analysis, region) for the 17 anchoring runs; mirrors 13.qtl_anchoring.sh.
ANALYSES: list[tuple[str, str | None]] = [
    ("brainseq-sczd", None),
    ("brainseq-aging", "caudate"),
    ("brainseq-aging", "hippocampus"),
    ("brainseq-aging", "dlpfc"),
    *[("gtex-aging", r) for r in [
        "amygdala", "anterior_cingulate_cortex_ba24", "caudate_basal_ganglia",
        "cerebellar_hemisphere", "cerebellum", "cortex", "frontal_cortex_ba9",
        "hippocampus", "hypothalamus", "nucleus_accumbens_basal_ganglia",
        "putamen_basal_ganglia", "spinal_cord_cervical_c_1", "substantia_nigra"]],
]
_Z = norm.ppf(0.975)


def _logor_se(row: pd.Series) -> float:
    """log-OR SE from the 95% CI, else from OR + p-value."""
    lo, hi = row.get("or_ci_low"), row.get("or_ci_high")
    if pd.notna(lo) and pd.notna(hi) and lo > 0 and hi > 0:
        return (np.log(hi) - np.log(lo)) / (2 * _Z)
    orr, p = row.get("odds_ratio"), row.get("pvalue")
    if pd.notna(orr) and orr > 0 and pd.notna(p) and 0 < p < 1:
        z = norm.ppf(1 - p / 2)
        if z > 0:
            return abs(np.log(orr)) / z
    return np.nan


def collect(variant: str, methods: tuple[str, ...]) -> pd.DataFrame:
    parts = []
    for method in methods:
        suffix = "" if method == "isograph" else f"_{method}"
        for analysis, region in ANALYSES:
            path = (_artifact_dir(analysis, region, variant).parent
                    / f"qtl_anchoring{suffix}.parquet")
            if path.exists():
                d = pd.read_parquet(path)
                if "graph_method" not in d.columns:
                    d["graph_method"] = method
                parts.append(d)
    if not parts:
        return pd.DataFrame()
    df = pd.concat(parts, ignore_index=True)
    df["beta"] = np.log(df["odds_ratio"])
    df["se"] = df.apply(_logor_se, axis=1)
    return df


_META_KEYS = meta_keys("or")


def _meta(group: pd.DataFrame) -> pd.Series:
    """Thin wrapper over the shared engine; kept so call sites read unchanged."""
    return meta(group, effect_name="or", count_col="n_foreground")


def _diff_rows(per: pd.DataFrame) -> pd.DataFrame:
    """Per (graph_method, analysis, module_set) sQTL-minus-eQTL log-OR difference.

    Removes the shared baseline (co-switch/network genes are cis-QTL-depleted for both
    QTL types) and isolates splicing specificity: ratio = sQTL OR / eQTL OR. The se
    treats the two estimates as independent, which is conservative because they share
    the foreground genes (positively correlated).
    """
    wide = per.pivot_table(index=["graph_method", "analysis", "region", "module_set"],
                           columns="xqtl_kind", values=["beta", "se"])
    rows = []
    for (gm, analysis, region, mset), r in wide.iterrows():
        bs, be = r[("beta", "sQTL")], r[("beta", "eQTL")]
        ss, se_ = r[("se", "sQTL")], r[("se", "eQTL")]
        if np.all(np.isfinite([bs, be, ss, se_])) and ss > 0 and se_ > 0:
            rows.append({"graph_method": gm, "analysis": analysis, "region": region,
                         "module_set": mset, "beta": bs - be,
                         "se": np.sqrt(ss**2 + se_**2)})
    return pd.DataFrame(rows)


def _contrast(diff: pd.DataFrame) -> pd.DataFrame:
    """Meta-analyse the splicing-specificity difference per (graph_method, module set)."""
    if diff.empty:
        return diff
    meta = (diff.groupby(["graph_method", "module_set"], sort=False)
            .apply(_meta, include_groups=False).reset_index())
    return meta.rename(columns={"or_fe": "ratio_fe", "or_fe_low": "ratio_fe_low",
                                "or_fe_high": "ratio_fe_high", "or_re": "ratio_re"})


_MSET_ORDER = {"all_modules": 0, "pheno_sig_modules": 1, "go_invisible_modules": 2,
               "go_visible_modules": 3}
_METHOD_ORDER = {"isograph": 0, "wgcna_switch_only": 1, "wgcna_multiplex": 2}


def _restrict_common(diff: pd.DataFrame) -> pd.DataFrame:
    """Keep only (analysis, region) covered by every graph_method, for a fair
    cross-method contrast (the matched WGCNA baselines run on fewer tissues)."""
    n_methods = diff["graph_method"].nunique()
    if n_methods < 2:
        return diff.iloc[0:0]
    cnt = diff.groupby(["analysis", "region"])["graph_method"].transform("nunique")
    return diff[cnt == n_methods]


def run_meta(variant: str = "standard",
             methods: tuple[str, ...] = ("isograph", "wgcna_switch_only",
                                         "wgcna_multiplex")) -> pd.DataFrame:
    per = collect(variant, methods)
    out_dir = ensure_dir(rel("real_data", "_m", "qtl_anchoring_meta"))
    if per.empty:
        print("no qtl_anchoring*.parquet outputs found")
        return per
    per.to_parquet(out_dir / "per_analysis.parquet", index=False, compression="zstd")
    by_mset = lambda s: s.map(_MSET_ORDER).fillna(9) if s.name == "module_set" else s
    by_method = lambda s: s.map(_METHOD_ORDER).fillna(9) if s.name == "graph_method" else s
    meta = (per.groupby(["graph_method", "xqtl_kind", "module_set"], sort=False)
            .apply(_meta, include_groups=False).reset_index())
    meta = meta.sort_values(["graph_method", "module_set", "xqtl_kind"],
                            key=lambda s: by_method(by_mset(s)))
    meta.to_parquet(out_dir / "qtl_anchoring_meta.parquet", index=False, compression="zstd")

    diff = _diff_rows(per)
    contrast = _contrast(diff)
    contrast_common = _contrast(_restrict_common(diff))
    if not contrast.empty:
        contrast = contrast.sort_values(["graph_method", "module_set"],
                                        key=lambda s: by_method(by_mset(s)))
        contrast.to_parquet(out_dir / "qtl_anchoring_meta_contrast.parquet",
                            index=False, compression="zstd")
    if not contrast_common.empty:
        contrast_common = contrast_common.sort_values(["module_set", "graph_method"],
                                                      key=lambda s: by_method(by_mset(s)))
        contrast_common.to_parquet(out_dir / "qtl_anchoring_meta_contrast_common.parquet",
                                   index=False, compression="zstd")
    _write_report(out_dir, per, meta, contrast, contrast_common)
    return meta


def _markdown_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    rows = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    rows += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(rows)


def _contrast_table(contrast: pd.DataFrame) -> pd.DataFrame:
    csh = contrast.copy()
    for c in ("ratio_fe", "ratio_fe_low", "ratio_fe_high", "ratio_re", "I2"):
        csh[c] = csh[c].round(3)
    csh["p_fe"] = csh["p_fe"].apply(lambda p: f"{p:.2e}" if pd.notna(p) else "NA")
    csh["k"] = csh["k"].astype(int)
    return csh[["graph_method", "module_set", "k", "ratio_fe", "ratio_fe_low",
                "ratio_fe_high", "p_fe", "ratio_re", "I2"]]


def _write_report(out_dir: Path, per: pd.DataFrame, meta: pd.DataFrame,
                  contrast: pd.DataFrame, contrast_common: pd.DataFrame) -> None:
    show = meta.copy()
    for c in ("or_fe", "or_fe_low", "or_fe_high", "or_re", "or_re_low", "or_re_high", "I2"):
        show[c] = show[c].round(2)
    for c in ("p_fe", "p_re"):
        show[c] = show[c].apply(lambda p: f"{p:.2e}" if pd.notna(p) else "NA")
    show["k"] = show["k"].astype(int)
    table = show[["graph_method", "module_set", "xqtl_kind", "k", "n_fg_total", "or_fe",
                  "or_fe_low", "or_fe_high", "p_fe", "or_re", "p_re", "I2"]]

    methods = list(per["graph_method"].unique())
    n_analyses = per[["analysis", "region"]].drop_duplicates().shape[0]
    n_common = (per[per["graph_method"] != "isograph"][["analysis", "region"]]
                .drop_duplicates().shape[0]) if len(methods) > 1 else 0
    lines = [
        "# Cross-tissue meta-analysis — co-switch module xQTL anchoring",
        "",
        f"Inverse-variance meta-analysis of the matched module-membership log-OR across "
        f"{n_analyses} brain region/cohort analyses, per (graph method, module set, xQTL "
        "kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = "
        f"heterogeneity. Graph methods: {', '.join(methods)}.",
        "",
        "Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the "
        "qtl_anchoring arrays for each method complete).",
        "",
        "## Pooled odds ratios (module genes vs background, per QTL type)",
        "",
        _markdown_table(table),
        "",
        "## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)",
        "",
        _markdown_table(_contrast_table(contrast)),
    ]
    if not contrast_common.empty:
        lines += [
            "",
            f"## Cross-method contrast on the {n_common} tissues all methods share",
            "",
            "The matched WGCNA baselines (`wgcna_switch_only`, `wgcna_multiplex`) consume "
            "the SAME switch / switch+abundance features as IsoGraph, so comparing their "
            "splicing-specificity on the SAME tissues isolates whether the signal lives in "
            "the switch features or in IsoGraph's VAE + Leiden inference.",
            "",
            _markdown_table(_contrast_table(contrast_common)),
        ]
    lines += [
        "",
        "## Reading",
        "",
        "- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) "
        "— expected: coordinated/network genes are more constrained and carry fewer "
        "common-variant cis-QTL. This shared baseline is NOT the result.",
        "- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, "
        "strongest for the GO-invisible modules — splicing-QTL is spared relative to "
        "expression-QTL exactly where the DTU-without-DGE value concentrates.",
        "- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL "
        "signal is a property of the switch features (which both methods share), not of "
        "IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. "
        "IsoGraph's value is the switch-feature *representation* and its finer modules, not "
        "a unique network-inference effect.",
        "- **Primary internal control = the matched WGCNA baselines**, not "
        "`go_visible_modules`. The baselines hold the switch features fixed and vary only "
        "the inference, so a null there localises the effect to IsoGraph's inference. "
        "`go_visible_modules` is a secondary control on module CONTENT and is only ever a "
        "relative contrast — read it as the low end of a gradient, not as an on/off null.",
        "- High I2 flags between-tissue heterogeneity; prefer RE there. A nominally "
        "significant ratio carrying high I2 is driven by a few tissues, not by a "
        "consistent effect, and is weaker evidence than a smaller ratio at I2 near 0 — "
        "compare module sets on consistency as well as magnitude. The contrast se "
        "is conservative (treats sQTL/eQTL estimates as independent though they share "
        "the foreground genes).",
        "- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the "
        "co-switching coordination itself.",
    ]
    (out_dir / "QTL_ANCHORING_META.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="Cross-tissue meta-analysis of xQTL anchoring.")
    p.add_argument("--variant", default="standard")
    p.add_argument("--methods", nargs="+",
                   default=["isograph", "wgcna_switch_only", "wgcna_multiplex"])
    args = p.parse_args()
    meta = run_meta(args.variant, tuple(args.methods))
    if not meta.empty:
        for kind in ("sQTL", "eQTL"):
            sub = meta[meta.xqtl_kind == kind]
            print(f"\n{kind} pooled (FE):")
            for _, r in sub.iterrows():
                print(f"  {r['graph_method']:18} {r['module_set']:22} k={int(r['k']):2}  "
                      f"OR={r['or_fe']:.2f}  p={r['p_fe']:.1e}")


if __name__ == "__main__":
    main()
