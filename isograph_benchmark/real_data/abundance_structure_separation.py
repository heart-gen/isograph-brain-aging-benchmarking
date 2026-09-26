"""Abundance vs isoform-structure separation — inputs for the figSeparation panel.

Assembles the three tidy tables behind the "abundance and isoform structure carry
partially non-redundant information" figure. Nothing new is *modelled* here: it
reuses IsoGraph's own per-sample abundance / switch channels (feature_scores.parquet)
and the already-computed de-confounded incremental test — it only derives the light
per-gene orthogonality summary and pulls out one example gene.

Outputs land in <artifact_dir>/abundance_structure/:

1. axis_orthogonality.parquet — per gene with both channels: gene_id, pearson_r
   (abundance channel vs switch channel across samples), abs_r, n_transcripts.
   pearson_r mass near 0 ⇒ the two inferred axes are computationally separable.

2. incremental_summary.parquet — region/cohort × category counts pooled from each
   analysis's incremental_association/summary.json (composition_unique = switch adds
   phenotype signal beyond abundance ⇒ the separation adds information).

3. example_gene.parquet (+ example_gene_stats.json) — per-sample abundance z and
   switch score by phenotype group for one composition-unique gene (default: the
   smallest p_switch_given_abund with a clearly non-significant p_abund_given_switch,
   i.e. stable total abundance but real switching), for the concrete illustration.

Two further modes serve the *general* form of the separability claim, which must not
rest on the one analysis that happens to carry the example gene:

``--orthogonality-only``
    Writes only (1) for the named analysis/region. Run as an array over all 17 stores
    (``03_module_characterization/_h/03b``); it needs feature_scores.parquet and the
    bundle's sample table, nothing from the incremental test.

``--rollup``
    Reads every store's axis_orthogonality.parquet and writes, under the stage-03 ``_m``:
    ``axis_orthogonality_all.parquet`` (one row per gene per analysis),
    ``axis_orthogonality_summary.{parquet,csv}`` (per-analysis n, median |r|, quartiles,
    fraction |r| < 0.1 and > 0.5) and ``AXIS_ORTHOGONALITY.md``. The faceted panel A of
    figSeparation and the manuscript's range statement read the rollup, not one store.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data.incremental_association import (
    _channel_matrix,
    _covariate_cols,  # noqa: F401  (kept for parity; not required here)
    _load,
    _sample_cols,
)
from isograph_benchmark.real_data.run_models import GTEX_REGIONS
from isograph_benchmark.real_data.sweep_leiden import _artifact_dir

# analyses/regions whose incremental_association/summary.json feed Panel B
_INCREMENTAL_TARGETS: list[tuple[str, str | None]] = (
    [("brainseq-sczd", None)]
    + [("brainseq-aging", r) for r in ("caudate", "hippocampus", "dlpfc")]
    + [("gtex-aging", r) for r in GTEX_REGIONS]
)

_LABELS = {"brainseq-sczd": "BrainSEQ SCZD", "brainseq-aging": "BrainSEQ aging",
           "gtex-aging": "GTEx aging"}


# --------------------------------------------------------------------------- #
# 1. Per-gene abundance-vs-switch orthogonality
# --------------------------------------------------------------------------- #
def axis_orthogonality(fs: pd.DataFrame, bundle) -> pd.DataFrame:
    """Row-wise Pearson r between each gene's abundance and switch channels."""
    st = bundle.sample_table
    sample_ids = set(st["sample_id"].astype(str))
    samp = _sample_cols(fs, sample_ids)

    SW = _channel_matrix(fs, "switch", samp)
    AB = _channel_matrix(fs, "abundance", samp)
    genes = SW.index.intersection(AB.index)
    S = SW.loc[genes].to_numpy(float)
    A = AB.loc[genes].to_numpy(float)

    finite = np.isfinite(S).all(axis=1) & np.isfinite(A).all(axis=1)
    S, A, genes = S[finite], A[finite], genes[finite]

    Sc = S - S.mean(axis=1, keepdims=True)
    Ac = A - A.mean(axis=1, keepdims=True)
    denom = np.sqrt((Sc**2).sum(axis=1) * (Ac**2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (Sc * Ac).sum(axis=1) / denom
    keep = np.isfinite(r)  # drops constant rows (zero variance)

    nt = (fs[fs["feature_type"] == "switch"].set_index("gene_id")["n_transcripts"]
          .reindex(genes).to_numpy())
    out = pd.DataFrame({
        "gene_id": np.asarray(genes)[keep],
        "pearson_r": r[keep],
        "abs_r": np.abs(r[keep]),
        "n_transcripts": nt[keep],
    })
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# 2. Incremental-test category counts pooled across analyses (Panel B)
# --------------------------------------------------------------------------- #
def incremental_summary(variant: str) -> pd.DataFrame:
    rows = []
    for analysis, region in _INCREMENTAL_TARGETS:
        sj = _artifact_dir(analysis, region, variant) / "incremental_association" / "summary.json"
        if not sj.exists():
            continue
        g = json.loads(sj.read_text())["gene_level"]
        rows.append({
            "analysis": analysis,
            "cohort": _LABELS.get(analysis, analysis),
            "region": region or "caudate",
            "label": f"{_LABELS.get(analysis, analysis)}: {region or 'caudate'}",
            "n_tested": int(g["n_tested"]),
            "composition_unique": int(g["composition_unique"]),
            "abundance_unique": int(g["abundance_unique"]),
            "both": int(g["both"]),
            "neither": int(g["neither"]),
        })
    return pd.DataFrame(rows)


def incremental_effect_sizes(variant: str) -> pd.DataFrame:
    """Conditional effect-size distributions pooled across every incremental analysis.

    An FDR count ("34 SCZD genes", "43 caudate genes") is strongly sample-size dependent,
    so it cannot say how much unique information the switch channel carries.  This pools
    the per-analysis ``gene_level_effect_sizes`` blocks into one long table: the overall
    quantiles of each conditional partial R², and the same per FDR category, so the
    significant genes can be read against the ``neither`` noise floor of their own analysis.

    One row per (analysis, region, effect, category); ``category == "__all__"`` is the
    unstratified distribution.
    """
    rows = []
    for analysis, region in _INCREMENTAL_TARGETS:
        sj = _artifact_dir(analysis, region, variant) / "incremental_association" / "summary.json"
        if not sj.exists():
            continue
        payload = json.loads(sj.read_text())
        blocks = payload.get("gene_level_effect_sizes") or {}
        base = {
            "analysis": analysis,
            "cohort": _LABELS.get(analysis, analysis),
            "region": region or "caudate",
            "variant": variant,
            "composition_adjusted": bool(payload.get("composition_adjusted", False)),
        }
        for effect, entry in blocks.items():
            q = entry.get("quantiles") or {}
            rows.append({**base, "effect": effect, "category": "__all__",
                         "n": entry.get("n"),
                         "q10": q.get("0.1"), "q25": q.get("0.25"), "median": q.get("0.5"),
                         "q75": q.get("0.75"), "q90": q.get("0.9")})
            for cat, sub in (entry.get("by_category") or {}).items():
                rows.append({**base, "effect": effect, "category": cat,
                             "n": sub.get("n"), "q10": None, "q25": None,
                             "median": sub.get("median"), "q75": None, "q90": sub.get("q90")})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# 3. Example gene: stable total abundance, real isoform switch (Panel C)
# --------------------------------------------------------------------------- #
def _pick_example_gene(artifact_dir, gene: str | None) -> tuple[str, pd.Series]:
    cu = pd.read_parquet(artifact_dir / "composition_unique" / "genes.parquet")
    if gene is not None:
        hit = cu[cu["gene_id"].astype(str).str.startswith(gene)]
        if hit.empty:
            raise SystemExit(f"gene {gene!r} not in composition_unique/genes.parquet")
        row = hit.iloc[0]
        return str(row["gene_id"]), row
    # cleanest: strong switch, most clearly non-significant abundance
    ranked = cu.sort_values(
        ["p_abund_given_switch", "p_switch_given_abund"],
        ascending=[False, True]).reset_index(drop=True)
    row = ranked.iloc[0]
    return str(row["gene_id"]), row


def example_gene(analysis: str, fs: pd.DataFrame, bundle, artifact_dir,
                 gene: str | None) -> tuple[pd.DataFrame, dict]:
    gene_id, stats_row = _pick_example_gene(artifact_dir, gene)
    st = bundle.sample_table
    sample_ids = set(st["sample_id"].astype(str))
    samp = _sample_cols(fs, sample_ids)
    sti = st.set_index(st["sample_id"].astype(str)).loc[samp]

    AB = _channel_matrix(fs, "abundance", samp)
    SW = _channel_matrix(fs, "switch", samp)
    if gene_id not in AB.index or gene_id not in SW.index:
        raise SystemExit(f"example gene {gene_id} missing an IsoGraph channel")

    group_col = "Dx" if analysis == "brainseq-sczd" else "Age"
    df = pd.DataFrame({
        "sample_id": samp,
        "group": sti[group_col].astype(str).to_numpy() if group_col == "Dx"
        else sti[group_col].to_numpy(float),
        "abundance_z": AB.loc[gene_id, samp].to_numpy(float),
        "switch_score": SW.loc[gene_id, samp].to_numpy(float),
    })
    stats_out = {
        "gene_id": gene_id,
        "group_col": group_col,
        "p_switch_given_abund": float(stats_row["p_switch_given_abund"]),
        "p_abund_given_switch": float(stats_row["p_abund_given_switch"]),
        "fdr_switch_given_abund": float(stats_row["fdr_switch_given_abund"]),
        "fdr_abund_given_switch": float(stats_row["fdr_abund_given_switch"]),
    }
    return df, stats_out


def _write_effect_report(eff: pd.DataFrame) -> None:
    """EFFECT_SIZES.md — the conditional effect sizes behind the FDR counts."""
    lines = [
        "# Conditional effect sizes for the unique-switch-information test", "",
        "The de-confounded gene-level test asks whether the switch channel carries "
        "phenotype signal the abundance channel does not, and vice versa. Reporting only "
        "the FDR-significant count makes the answer a function of sample size. These are "
        "the effect-size distributions behind those counts: `partial_r2_switch_given_abund` "
        "is the proportion of residual variance the switch block explains after the "
        "covariates and abundance are already in the model (`partial_r2_abund_given_switch` "
        "is the mirror image).", "",
        "Read the FDR-significant categories against `neither`, which is the same "
        "analysis's own noise floor. `composition_unique` = switch-significant only; "
        "`abundance_unique` = abundance-significant only; `both` = both.", "",
    ]
    for effect, block in eff.groupby("effect"):
        lines += [f"## `{effect}`", "",
                  "| analysis | region | adj | n | median (all) | q90 (all) | "
                  "median: composition_unique | median: neither | ratio vs neither |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for keys, g in block.groupby(["analysis", "region", "composition_adjusted"]):
            analysis, region, adjusted = keys
            allrow = g[g["category"] == "__all__"]
            cu = g[g["category"] == "composition_unique"]["median"]
            ne = g[g["category"] == "neither"]["median"]
            cu_v = float(cu.iloc[0]) if len(cu) and pd.notna(cu.iloc[0]) else np.nan
            ne_v = float(ne.iloc[0]) if len(ne) and pd.notna(ne.iloc[0]) else np.nan
            med = float(allrow["median"].iloc[0]) if len(allrow) else np.nan
            q90 = float(allrow["q90"].iloc[0]) if len(allrow) else np.nan
            n = int(allrow["n"].iloc[0]) if len(allrow) and pd.notna(allrow["n"].iloc[0]) else 0
            ratio = cu_v / ne_v if np.isfinite(cu_v) and np.isfinite(ne_v) and ne_v > 0 else np.nan
            fmt = lambda v, s="{:.4f}": s.format(v) if np.isfinite(v) else "n/a"  # noqa: E731
            lines.append(
                f"| {analysis} | {region} | {'yes' if adjusted else 'no'} | {n:,} | "
                f"{fmt(med)} | {fmt(q90)} | {fmt(cu_v)} | {fmt(ne_v)} | "
                f"{fmt(ratio, '{:.1f}x')} |")
        lines.append("")
    out = stage_out("characterize", "EFFECT_SIZES.md")
    out.write_text("\n".join(lines) + "\n")


def _analysis_label(analysis: str, region: str | None) -> tuple[str, str, str]:
    """(cohort, region, label) exactly as incremental_summary() labels them."""
    cohort = _LABELS.get(analysis, analysis)
    reg = region or "caudate"
    return cohort, reg, f"{cohort}: {reg}"


def _summarise_orthogonality(ortho: pd.DataFrame) -> dict:
    a = ortho["abs_r"].to_numpy(float)
    return {
        "n_genes": int(len(a)),
        "median_abs_r": float(np.median(a)),
        "q25_abs_r": float(np.quantile(a, 0.25)),
        "q75_abs_r": float(np.quantile(a, 0.75)),
        "frac_abs_r_lt_0.1": float((a < 0.1).mean()),
        "frac_abs_r_lt_0.3": float((a < 0.3).mean()),
        "frac_abs_r_gt_0.5": float((a > 0.5).mean()),
        "median_r": float(ortho["pearson_r"].median()),
    }


def write_orthogonality(analysis: str, region: str | None, variant: str) -> pd.DataFrame:
    """Derive + write axis_orthogonality.parquet for one store; returns the table."""
    artifact_dir, fs, _modules, bundle = _load(analysis, region, variant)
    out = ensure_dir(artifact_dir / "abundance_structure")
    ortho = axis_orthogonality(fs, bundle)
    ortho.to_parquet(out / "axis_orthogonality.parquet", index=False, compression="zstd")
    s = _summarise_orthogonality(ortho)
    print(f"[{analysis}{'/' + region if region else ''}] axis_orthogonality: "
          f"n={s['n_genes']} genes, median|r|={s['median_abs_r']:.3f}, "
          f"frac|r|<0.1={s['frac_abs_r_lt_0.1']:.2f}")
    return ortho


def rollup_orthogonality(variant: str) -> pd.DataFrame:
    """Pool every store's axis_orthogonality.parquet; write the long table + summary."""
    long_rows, summ_rows, missing = [], [], []
    for analysis, region in _INCREMENTAL_TARGETS:
        p = _artifact_dir(analysis, region, variant) / "abundance_structure" / "axis_orthogonality.parquet"
        cohort, reg, label = _analysis_label(analysis, region)
        if not p.exists():
            missing.append(label)
            continue
        o = pd.read_parquet(p)
        o.insert(0, "label", label)
        o.insert(0, "region", reg)
        o.insert(0, "cohort", cohort)
        o.insert(0, "analysis", analysis)
        long_rows.append(o)
        summ_rows.append({"analysis": analysis, "cohort": cohort, "region": reg,
                          "label": label, **_summarise_orthogonality(o)})
    if not long_rows:
        raise SystemExit("rollup: no axis_orthogonality.parquet found in any store "
                         "(run --orthogonality-only over the stores first)")
    long = pd.concat(long_rows, ignore_index=True)
    summ = pd.DataFrame(summ_rows)

    out = ensure_dir(stage_out("characterize"))
    long.to_parquet(out / "axis_orthogonality_all.parquet", index=False, compression="zstd")
    summ.to_parquet(out / "axis_orthogonality_summary.parquet", index=False, compression="zstd")
    summ.to_csv(out / "axis_orthogonality_summary.csv", index=False, float_format="%.4f")
    _write_orthogonality_report(summ, missing, out)
    print(f"[rollup] {len(summ)} analyses, {len(long):,} gene rows; "
          f"median|r| {summ['median_abs_r'].min():.3f}-{summ['median_abs_r'].max():.3f}, "
          f"frac|r|<0.1 {summ['frac_abs_r_lt_0.1'].min():.2f}-{summ['frac_abs_r_lt_0.1'].max():.2f}"
          + (f"; MISSING: {', '.join(missing)}" if missing else ""))
    return summ


def _write_orthogonality_report(summ: pd.DataFrame, missing: list[str], out) -> None:
    """AXIS_ORTHOGONALITY.md — the per-analysis separability numbers behind figSeparation A."""
    s = summ
    lines = [
        "# Abundance-vs-switch axis orthogonality, all analyses", "",
        "Per gene, the Pearson correlation across samples between IsoGraph's abundance "
        "channel and its switch channel (`feature_scores.parquet`), computed by "
        "`abundance_structure_separation.py --orthogonality-only` for every store and pooled "
        "by `--rollup`. Genes with a constant channel are dropped. The manuscript's "
        "separability statement quotes the range of the per-analysis medians and of the "
        "fraction of genes with |r| < 0.1; it does not quote one analysis.", "",
        f"- analyses: {len(s)} of {len(_INCREMENTAL_TARGETS)}"
        + (f" (missing: {', '.join(missing)})" if missing else ""),
        f"- genes per analysis: {int(s['n_genes'].min()):,}-{int(s['n_genes'].max()):,}",
        f"- median |r|: {s['median_abs_r'].min():.3f}-{s['median_abs_r'].max():.3f} "
        f"(median of medians {s['median_abs_r'].median():.3f})",
        f"- fraction |r| < 0.1: {s['frac_abs_r_lt_0.1'].min():.2f}-{s['frac_abs_r_lt_0.1'].max():.2f}",
        f"- fraction |r| < 0.3: {s['frac_abs_r_lt_0.3'].min():.2f}-{s['frac_abs_r_lt_0.3'].max():.2f}",
        f"- fraction |r| > 0.5: {s['frac_abs_r_gt_0.5'].min():.3f}-{s['frac_abs_r_gt_0.5'].max():.3f}",
        "",
        "| analysis | region | n genes | median \\|r\\| | IQR | \\|r\\| < 0.1 | \\|r\\| < 0.3 | \\|r\\| > 0.5 | median r |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for _, r in s.iterrows():
        lines.append(
            f"| {r['analysis']} | {r['region']} | {int(r['n_genes']):,} | {r['median_abs_r']:.3f} | "
            f"{r['q25_abs_r']:.3f}-{r['q75_abs_r']:.3f} | {r['frac_abs_r_lt_0.1']:.2f} | "
            f"{r['frac_abs_r_lt_0.3']:.2f} | {r['frac_abs_r_gt_0.5']:.3f} | {r['median_r']:+.3f} |")
    lines += ["", "Outputs: `axis_orthogonality_all.parquet` (per gene per analysis), "
              "`axis_orthogonality_summary.{parquet,csv}` (this table). Figure: "
              "`manuscript/_h/abundance_structure_figure.R` panel A (S-real-4).", ""]
    (out / "AXIS_ORTHOGONALITY.md").write_text("\n".join(lines))


def run(analysis: str, region: str | None, variant: str, gene: str | None) -> None:
    artifact_dir, fs, _modules, bundle = _load(analysis, region, variant)
    out = ensure_dir(artifact_dir / "abundance_structure")

    ortho = axis_orthogonality(fs, bundle)
    ortho.to_parquet(out / "axis_orthogonality.parquet", index=False, compression="zstd")

    incr = incremental_summary(variant)
    incr.to_parquet(out / "incremental_summary.parquet", index=False, compression="zstd")

    # Cross-analysis product, so it lands at the repo _m level rather than under one
    # region's artifact dir.
    eff = incremental_effect_sizes(variant)
    if not eff.empty:
        eff_out = ensure_dir(stage_out("characterize")) / "incremental_effect_sizes.parquet"
        eff.to_parquet(eff_out, index=False, compression="zstd")
        _write_effect_report(eff)

    eg, eg_stats = example_gene(analysis, fs, bundle, artifact_dir, gene)
    eg.to_parquet(out / "example_gene.parquet", index=False, compression="zstd")
    (out / "example_gene_stats.json").write_text(json.dumps(eg_stats, indent=2))

    med = float(ortho["abs_r"].median())
    frac_lt = float((ortho["abs_r"] < 0.1).mean())
    print(f"[{analysis}] axis_orthogonality: n={len(ortho)} genes, "
          f"median|r|={med:.3f}, frac|r|<0.1={frac_lt:.2f}")
    print(f"[{analysis}] incremental_summary: {len(incr)} analyses, "
          f"total composition_unique={int(incr['composition_unique'].sum())}")
    print(f"[{analysis}] example gene {eg_stats['gene_id']}: "
          f"p_switch|abund={eg_stats['p_switch_given_abund']:.2e}, "
          f"p_abund|switch={eg_stats['p_abund_given_switch']:.2f}")
    print(f"[{analysis}] written to {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis", nargs="?", default="brainseq-sczd",
                        choices=["brainseq-sczd", "brainseq-aging", "gtex-aging"])
    parser.add_argument("--region", default=None,
                        help="Region for aging analyses (default caudate).")
    parser.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    parser.add_argument("--gene", default=None,
                        help="Example gene id/prefix (default: auto-pick the cleanest "
                             "composition-unique gene).")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--orthogonality-only", action="store_true",
                      help="Write only axis_orthogonality.parquet for this analysis/region "
                           "(the per-store step of the all-analyses array).")
    mode.add_argument("--rollup", action="store_true",
                      help="Pool every store's axis_orthogonality.parquet into the stage _m "
                           "(ignores the analysis/region arguments).")
    args = parser.parse_args()
    if args.rollup:
        rollup_orthogonality(args.variant)
        return
    region = args.region
    if args.analysis != "brainseq-sczd" and region is None:
        region = "caudate"
    if args.orthogonality_only:
        write_orthogonality(args.analysis, region, args.variant)
        return
    run(args.analysis, region, args.variant, args.gene)


if __name__ == "__main__":
    main()
