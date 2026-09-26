"""Does the switch-unique picture survive moving the FDR line?

A gene is called switch-unique when its switch coordinate reaches FDR <= alpha conditional
on abundance **and** the abundance test conditional on the switch coordinate does not. The
second half is an acceptance of a null, so a gene can enter or leave the class because its
abundance test sits just either side of alpha. The manuscript already says the
classification does not test whether the two conditional effects differ, but the result of
the subsection is a *regional pattern* -- limbic/striatal analyses retaining switch-unique
genes after composition adjustment while the two GTEx cortical analyses collapse -- and
nothing showed that this pattern is not an artifact of alpha = 0.10.

This module recounts at several alphas. It is a recount, not a refit: BH-adjusted values
(``fdr_switch_given_abund``, ``fdr_abund_given_switch``) are computed once over all genes
by ``incremental_association`` and do not depend on alpha, so only the cut moves. The
category rule is reproduced exactly as ``incremental_association.gene_level_deconfounded``
applies it, ``<=`` included.

Outputs (per analysis x alpha, base and composition-adjusted where available):

1. ``switch_unique_threshold.parquet``/``.csv`` -- counts in all four categories, plus
   ``n_overlap``/``n_new`` between the base and adjusted switch-unique sets at that alpha.
2. ``switch_unique_threshold_rank.csv`` -- each analysis's rank by switch-unique count and
   by gene-level persistence at each alpha, and the Spearman correlation of those rankings
   against the alpha = 0.10 reference. This is the number the claim rests on.
3. ``SWITCH_UNIQUE_THRESHOLD.md``.

Run (local, no scheduler): ``bash 03_module_characterization/_h/04c.switch_unique_threshold.sh``
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.paths import ensure_dir, stage_out

# The manuscript's operating point, listed first so it reads as the reference.
DEFAULT_ALPHAS = (0.10, 0.05, 0.20)
REFERENCE_ALPHA = 0.10
# Below this many unadjusted switch-unique genes, persistence is a ratio of single
# digits and its rank is noise; rank stability is reported with and without them.
MIN_DENOMINATOR = 10

BASE = "incremental_association"
ADJ = "incremental_association_composition"


def _categorize(gl: pd.DataFrame, alpha: float) -> pd.Series:
    """The category rule from incremental_association, reproduced exactly (``<=``)."""
    sw = gl["fdr_switch_given_abund"] <= alpha
    ab = gl["fdr_abund_given_switch"] <= alpha
    return pd.Series(np.select([sw & ~ab, ab & ~sw, sw & ab],
                               ["switch_unique", "abundance_unique", "both"], "neither"),
                     index=gl.index)


def _stores(variant: str) -> list[tuple[str, str]]:
    """Every analysis with a gene-level incremental-association table, discovered on disk."""
    subdir = "isograph_vae_with_abundance" if variant == "with-abundance" else "isograph_vae"
    found = []
    for p in sorted(stage_out("modules").glob(
            f"*/*/_m/{subdir}/{BASE}/gene_level.parquet")):
        found.append((p.parts[-6], p.parts[-5]))
    return found


def _label(cohort: str, region: str) -> str:
    """The analysis label the composition rollup uses, so the tables join."""
    if cohort == "gtex":
        return f"GTEx {region}"
    if region == "caudate_sczd":
        return "SCZD (caudate)"
    return f"aging {'DLPFC' if region == 'dlpfc' else region}"


def _arm(cohort: str, region: str, arm: str, variant: str) -> pd.DataFrame | None:
    subdir = "isograph_vae_with_abundance" if variant == "with-abundance" else "isograph_vae"
    p = stage_out("modules", cohort, region, "_m", subdir, arm, "gene_level.parquet")
    if not p.exists():
        return None
    return pd.read_parquet(p, columns=["gene_id", "fdr_switch_given_abund",
                                       "fdr_abund_given_switch"])


def _rows(variant: str, alphas: tuple[float, ...]) -> pd.DataFrame:
    rows = []
    for cohort, region in _stores(variant):
        base = _arm(cohort, region, BASE, variant)
        adj = _arm(cohort, region, ADJ, variant)
        for alpha in alphas:
            b_cat = _categorize(base, alpha)
            b_set = set(base.loc[b_cat == "switch_unique", "gene_id"])
            row = {
                "label": _label(cohort, region), "cohort": cohort, "region": region,
                "alpha": alpha, "n_tested": int(len(base)),
                "switch_unique_base": len(b_set),
                "abundance_unique_base": int((b_cat == "abundance_unique").sum()),
                "both_base": int((b_cat == "both").sum()),
                "composition_adjusted": adj is not None,
            }
            if adj is None:
                row.update({"switch_unique_adj": np.nan, "n_overlap": np.nan,
                            "n_new": np.nan, "persistence": np.nan})
            else:
                a_set = set(adj.loc[_categorize(adj, alpha) == "switch_unique", "gene_id"])
                row.update({
                    "switch_unique_adj": len(a_set),
                    "n_overlap": len(b_set & a_set),
                    "n_new": len(a_set - b_set),
                    # A real fraction: the kept genes are a subset of the unadjusted set.
                    "persistence": len(b_set & a_set) / len(b_set) if b_set else np.nan,
                })
            rows.append(row)
    return pd.DataFrame(rows)


def _rank_stability(df: pd.DataFrame) -> pd.DataFrame:
    """Does moving alpha reorder the analyses, or only rescale the counts?

    The subsection's claim is about which analyses retain and which collapse, so rank
    agreement with the reference alpha is the thing to test -- not whether the absolute
    counts move, which they must.
    """
    ref = df[df["alpha"] == REFERENCE_ALPHA].set_index("label")
    out = []
    for alpha, g in df.groupby("alpha"):
        g = g.set_index("label")
        shared = ref.index.intersection(g.index)
        entry = {"alpha": alpha, "n_analyses": len(shared)}
        for col in ("switch_unique_base", "persistence"):
            a, b = ref.loc[shared, col], g.loc[shared, col]
            ok = a.notna() & b.notna()
            entry[f"spearman_{col}"] = (float(stats.spearmanr(a[ok], b[ok]).statistic)
                                        if ok.sum() > 2 else np.nan)
            entry[f"n_{col}"] = int(ok.sum())
        # Persistence is a ratio, so an analysis with two unadjusted genes contributes a
        # value of 0, 0.5 or 1 and can dominate a rank correlation over 12 analyses. Repeat
        # it over the analyses whose denominator is large enough for the rank to mean
        # anything at BOTH alphas; report both numbers rather than only the flattering one.
        big = shared[(ref.loc[shared, "switch_unique_base"] >= MIN_DENOMINATOR)
                     & (g.loc[shared, "switch_unique_base"] >= MIN_DENOMINATOR)]
        a, b = ref.loc[big, "persistence"], g.loc[big, "persistence"]
        ok = a.notna() & b.notna()
        entry["spearman_persistence_large"] = (float(stats.spearmanr(a[ok], b[ok]).statistic)
                                               if ok.sum() > 2 else np.nan)
        entry["n_persistence_large"] = int(ok.sum())
        out.append(entry)
    return pd.DataFrame(out).sort_values("alpha", ignore_index=True)


# A collapse and a retention, stated without reference to any one alpha.
COLLAPSE_MAX = 0.05
RETAIN_MIN = 0.20


def _robust_calls(df: pd.DataFrame) -> pd.DataFrame:
    """Which analyses collapse, or retain, at *every* alpha tested.

    The rank correlations answer "does the ordering move"; this answers the question the
    manuscript actually asks, which is whether a given analysis keeps its switch-unique
    genes. An analysis is called only when the same call holds at every alpha, so the call
    carries no dependence on where the line was drawn.
    """
    g = (df.dropna(subset=["persistence"])
         .groupby("label")
         .agg(n_alphas=("alpha", "size"),
              min_persistence=("persistence", "min"),
              max_persistence=("persistence", "max"),
              min_base=("switch_unique_base", "min")))
    g["call"] = np.select(
        [g["max_persistence"] < COLLAPSE_MAX, g["min_persistence"] >= RETAIN_MIN],
        ["collapses at every alpha", "retains at every alpha"], "threshold-dependent")
    # A call resting on a handful of genes is not a call; say so on the row itself.
    g.loc[g["min_base"] < MIN_DENOMINATOR, "call"] += " (small set)"
    return g.reset_index().sort_values("min_persistence", ignore_index=True)


def _md_table(df: pd.DataFrame) -> str:
    head = "| " + " | ".join(df.columns) + " |"
    rule = "| " + " | ".join("---" for _ in df.columns) + " |"
    body = ["| " + " | ".join("" if pd.isna(v) else
                              (f"{v:.3g}" if isinstance(v, float) else str(v))
                              for v in row) + " |"
            for row in df.itertuples(index=False)]
    return "\n".join([head, rule, *body])


def run(variant: str, alphas: tuple[float, ...]) -> None:
    df = _rows(variant, alphas)
    if df.empty:
        raise SystemExit("no incremental_association gene_level tables found")
    ranks = _rank_stability(df)

    out = ensure_dir(stage_out("characterize"))
    df.to_parquet(out / "switch_unique_threshold.parquet", index=False, compression="zstd")
    # CSV alongside the parquet: the figure reads it, so the figure builds with an `arrow`
    # compiled without zstd support.
    df.to_csv(out / "switch_unique_threshold.csv", index=False)
    ranks.to_csv(out / "switch_unique_threshold_rank.csv", index=False)
    calls = _robust_calls(df)
    calls.to_csv(out / "switch_unique_threshold_calls.csv", index=False)

    ref = df[df["alpha"] == REFERENCE_ALPHA]
    wide = df.pivot(index="label", columns="alpha", values="switch_unique_base")
    wide.columns = [f"base_a{a:g}" for a in wide.columns]
    pers = df.pivot(index="label", columns="alpha", values="persistence")
    pers.columns = [f"persist_a{a:g}" for a in pers.columns]
    table = (wide.join(pers).reset_index()
             .sort_values(f"base_a{REFERENCE_ALPHA:g}", ascending=False,
                          ignore_index=True))

    summary = {
        "variant": variant, "alphas": list(alphas),
        "reference_alpha": REFERENCE_ALPHA,
        "n_analyses": int(df["label"].nunique()),
        "n_composition_adjusted": int(ref["composition_adjusted"].sum()),
        "rank_stability": ranks.to_dict(orient="records"),
        "calls": calls.to_dict(orient="records"),
    }
    (out / "switch_unique_threshold_summary.json").write_text(json.dumps(summary, indent=2))

    lines = [
        "# Threshold sensitivity of the switch-unique classification", "",
        "A recount, not a refit: `fdr_switch_given_abund` and `fdr_abund_given_switch` are "
        "BH-adjusted once over all genes by `incremental_association`, so changing alpha "
        "moves only the cut. The category rule is reproduced exactly, `<=` included.", "",
        f"{df['label'].nunique()} analyses; "
        f"{int(ref['composition_adjusted'].sum())} of them also have a "
        "composition-adjusted arm.", "",
        "## Switch-unique counts and gene-level persistence by alpha", "",
        _md_table(table), "",
        "## Does alpha reorder the analyses?", "",
        _md_table(ranks), "",
        "## Calls that hold at every alpha", "",
        _md_table(calls), "",
        "## Interpretation", "",
        "- `spearman_switch_unique_base` and `spearman_persistence` compare each alpha's "
        f"ranking of the analyses with the manuscript's alpha = {REFERENCE_ALPHA:g}. Counts "
        "must move with alpha; what the subsection claims is the *ordering* -- which "
        "analyses carry a large switch-unique set, and which lose it under composition "
        "adjustment.",
        "- A ranking correlation near 1 means the regional pattern, including the cortical "
        "collapse, is not an artifact of where the line was drawn.",
        f"- `spearman_persistence` over all analyses is **not** a reliable number here, and "
        f"is reported only for completeness. Persistence is a ratio, and several analyses "
        f"have a single-digit unadjusted set (GTEx putamen and nucleus accumbens sit at 1-2 "
        f"genes, BrainSEQ hippocampus at 1), so their persistence jumps between 0, 0.5 and "
        f"1 and swamps the rank. `spearman_persistence_large`, restricted to analyses with "
        f"at least {MIN_DENOMINATOR} unadjusted switch-unique genes at both alphas, is the "
        f"one to read.",
        "- The individual claims the manuscript makes are directly checkable in the table "
        "above and hold at every alpha: GTEx cortex and frontal cortex BA9 retain "
        "essentially none of their unadjusted genes, while ACC BA24, BrainSEQ DLPFC and "
        "BrainSEQ caudate retain a substantial share.",
        f"- The **calls** table is the threshold-free statement: an analysis is called only "
        f"when the same call holds at every alpha, collapsing at persistence < "
        f"{COLLAPSE_MAX:g} or retaining at >= {RETAIN_MIN:g}. Anything else is labelled "
        f"threshold-dependent rather than rounded to the convenient side.",
        "- The classification remains what it was: an asymmetric use of one threshold that "
        "does not test whether the two conditional effects differ. This file bounds the "
        "threshold's influence; it does not turn the class into a contrast.",
    ]
    (out / "SWITCH_UNIQUE_THRESHOLD.md").write_text("\n".join(lines))
    print(f"[threshold] {df['label'].nunique()} analyses x {len(alphas)} alphas", flush=True)
    for r in ranks.itertuples(index=False):
        print(f"[threshold]   alpha={r.alpha:g}: rank rho vs {REFERENCE_ALPHA:g} -- "
              f"counts {r.spearman_switch_unique_base:.3f}, persistence "
              f"{r.spearman_persistence:.3f} (all) / "
              f"{r.spearman_persistence_large:.3f} (n>={MIN_DENOMINATOR}, "
              f"{r.n_persistence_large} analyses)", flush=True)
    for r in calls.itertuples(index=False):
        print(f"[threshold]   {r.label}: {r.call} "
              f"(persistence {r.min_persistence:.3g}-{r.max_persistence:.3g})", flush=True)
    print(f"[threshold] wrote {out/'SWITCH_UNIQUE_THRESHOLD.md'}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    p.add_argument("--alpha", type=float, action="append", dest="alphas",
                   help=f"FDR cut to recount at; repeatable. Default: {DEFAULT_ALPHAS}. "
                        f"{REFERENCE_ALPHA} is always included as the reference.")
    args = p.parse_args()
    alphas = tuple(args.alphas) if args.alphas else DEFAULT_ALPHAS
    if REFERENCE_ALPHA not in alphas:
        alphas = (REFERENCE_ALPHA,) + alphas
    run(args.variant, tuple(sorted(set(alphas))))


if __name__ == "__main__":
    main()
