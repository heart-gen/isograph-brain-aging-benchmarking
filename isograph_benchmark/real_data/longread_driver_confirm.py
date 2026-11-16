"""Long-read confirmation of the module *driver* switches.

``longread_switch_confirm`` (06/01e) scores every coding switch pair of the age-associated
co-switch modules on the ONT DLPFC matrix (Aguzzoli-Heberle 2024 NBT, Bambu, Zenodo
8180677, n = 12), and ``switch_orthogonal_confirm`` (06/02b) supplies the abundance-matched
nulls for the whole set and for the genetically anchored subset. Neither says anything about
the transcripts the manuscript actually names: the per-module **driver transcripts**, the
top-k isoforms by |r| in each module's ``transcript_polarity_table`` (the same rule
``module_trust._module_drivers`` uses for driver reproducibility). This CLI closes that gap
so a "x% of driver switches replicate in long read" statement has a defined numerator,
denominator and null behind it.

Unit and question
-----------------
The unit is the **distinct driver pair**: a coding switch pair already scored in long read
(``pair_confirmation.parquet``) of which at least one member is a driver transcript of an
age-associated module (eigengene spline F-test FDR <= ``--alpha``) in the same GTEx cortical
region. The two discovery regions (frontal_cortex_ba9, cortex) are scored on the same twelve
long-read samples, so a pair nominated in both regions carries one outcome and is counted
once; its regions are kept as a list.

Three quantities are reported, each with the comparison it needs:

  1. **Existence** -- driver transcripts detected in long read (count >= ``min_count`` in
     >= ``min_samples`` samples) and driver pairs with both members detected, against the
     detection rate of the non-driver switch pairs.
  2. **Switch behaviour** -- the switch-like rate (both members detected, within-gene usage
     negatively rank-correlated across samples; the 02b definition) of detected driver pairs
     against an abundance-decile-matched null drawn from non-driver switch pairs, with
     replacement, ``--n-draws`` times. This is the same ``matched_null`` the anchored arm
     uses, so the two subsets are comparable.
  3. **Usable-abundance qualification** -- the same test restricted to pairs whose driver
     member reaches ``--min-if`` mean isoform fraction, the analogue of ``anchored_usable``.

Two backgrounds are drawn: all non-driver switch pairs (primary; the anchored arm's
universe) and non-driver pairs from the driver genes only (secondary; small, but it controls
for gene). Per-region focal counts are descriptive only.

Outputs land in ``06_switch_mechanism/_m/longread_driver_confirm/`` (``sensitivity/k<K>/``
for a non-default ``--k``):
  * ``driver_transcripts.parquet`` -- per region x module x driver: rank, r, in the coding
    switch-pair universe, detected in long read.
  * ``driver_pair_confirmation.parquet`` -- one row per distinct driver pair with its
    long-read outcome, driver side(s), regions and modules.
  * ``summary.json`` -- counts, rates, matched nulls, thresholds.
  * ``LONGREAD_DRIVER_CONFIRM.md`` -- the tables, written from ``summary.json``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, region_store, stage_out
from isograph_benchmark.real_data.longread_null import _md_table, matched_null

# Deliberately not imported from longread_switch_confirm / validate_switch_splicing: that
# chain needs the IsoGraph package, which this CLI (local, post hoc) should not require.
_DEFAULT_REGIONS = ("frontal_cortex_ba9", "cortex")   # == longread_switch_confirm._DEFAULT_REGIONS
_PAIR_CONFIRMATION = stage_out("mechanism", "longread_switch_confirm", "pair_confirmation.parquet")
DEFAULT_K = 5


def _strip_ver(s: pd.Series) -> pd.Series:
    """Ensembl id without version suffix (== validate_switch_splicing._strip_ver)."""
    return s.astype(str).str.replace(r"\.\d+$", "", regex=True)


def _artifact_dir(region: str, variant: str) -> Path:
    """GTEx per-region artifact store (== validate_switch_splicing._artifact_dir for gtex)."""
    subdir = "isograph_vae_with_abundance" if variant == "with-abundance" else "isograph_vae"
    return region_store("gtex", region, subdir)
_LONGREAD_DATASET = "Aguzzoli-Heberle 2024 NBT; Bambu; Zenodo 8180677; n=12 DLPFC BA9/46"


def out_path(k: int) -> Path:
    base = stage_out("mechanism", "longread_driver_confirm")
    return base if k == DEFAULT_K else base / "sensitivity" / f"k{k}"


# --------------------------------------------------------------------------- #
# Driver transcripts
# --------------------------------------------------------------------------- #
def age_modules(art: Path, alpha: float) -> list[str]:
    """Modules age-associated at eigengene spline F-test FDR <= alpha."""
    sp = pd.read_parquet(art / "age_spline.parquet", columns=["module_id", "fdr_ftest"])
    return sorted(sp.loc[sp["fdr_ftest"] <= alpha, "module_id"].unique())


def top_k_drivers(polarity: pd.DataFrame, k: int) -> pd.DataFrame:
    """Top-k transcripts by |r| of one module's transcript_polarity_table (unversioned ids)."""
    df = polarity.dropna(subset=["r"]).copy()
    df["abs_r"] = pd.to_numeric(df["r"], errors="coerce").abs()
    df = df.dropna(subset=["abs_r"]).sort_values("abs_r", ascending=False).head(k)
    return pd.DataFrame({
        "transcript_id": _strip_ver(df["transcript_id"].astype(str)).to_numpy(),
        "gene": _strip_ver(df["gene_id"].astype(str)).to_numpy(),
        "r": df["r"].astype(float).to_numpy(),
        "rank": np.arange(1, len(df) + 1),
    })


def driver_transcripts(region: str, k: int, alpha: float, variant: str) -> pd.DataFrame:
    art = _artifact_dir(region, variant)
    frames = []
    for m in age_modules(art, alpha):
        tp = art / "module_interpret" / m / "transcript_polarity_table.parquet"
        if not tp.exists():
            continue
        d = top_k_drivers(pd.read_parquet(tp, columns=["transcript_id", "gene_id", "r"]), k)
        d.insert(0, "module_id", m)
        d.insert(0, "region", region)
        frames.append(d)
    if not frames:
        return pd.DataFrame(columns=["region", "module_id", "transcript_id", "gene", "r", "rank"])
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------------------- #
# Driver pairs
# --------------------------------------------------------------------------- #
def load_switch_pairs(path: Path = _PAIR_CONFIRMATION) -> pd.DataFrame:
    """Coding switch pairs with their long-read outcome, one row per region x pair."""
    if not path.exists():
        raise SystemExit(f"missing {path}; run longread_switch_confirm first")
    pc = pd.read_parquet(path).copy()
    pc["max_mean_if"] = pc[["t1_mean_if", "t2_mean_if"]].max(axis=1)
    # 02b definition; the 01e column also demands both IF >= 0.05 and is kept as *_strict.
    pc = pc.rename(columns={"switch_like": "switch_like_strict"})
    pc["switch_like"] = pc["pair_detected"] & (pc["usage_spearman"] < 0)
    return pc


def flag_driver_pairs(pairs: pd.DataFrame, drivers: pd.DataFrame, min_if: float) -> pd.DataFrame:
    """Per region, mark pairs with a driver member; record the driver side's abundance."""
    out = pairs.copy()
    key = drivers.groupby(["region", "transcript_id"])["module_id"].agg(
        lambda s: ";".join(sorted(set(s)))
    )
    mod1 = key.reindex(pd.MultiIndex.from_arrays([out["region"], out["t1"]])).to_numpy()
    mod2 = key.reindex(pd.MultiIndex.from_arrays([out["region"], out["t2"]])).to_numpy()
    out["driver_t1"] = pd.notna(mod1)
    out["driver_t2"] = pd.notna(mod2)
    out["driver_pair"] = out["driver_t1"] | out["driver_t2"]
    out["driver_modules"] = [
        ";".join(sorted({*(str(a).split(";") if isinstance(a, str) else []),
                         *(str(b).split(";") if isinstance(b, str) else [])}))
        for a, b in zip(mod1, mod2)
    ]
    out["driver_mean_if"] = np.where(
        out["driver_t1"] & out["driver_t2"], out["max_mean_if"],
        np.where(out["driver_t1"], out["t1_mean_if"],
                 np.where(out["driver_t2"], out["t2_mean_if"], np.nan)))
    out["driver_usable"] = out["driver_pair"] & (out["driver_mean_if"] >= min_if)
    return out


def distinct_pairs(flagged: pd.DataFrame) -> pd.DataFrame:
    """Collapse region rows to one row per unordered transcript pair.

    The long-read samples are the same for every region, so the outcome is identical across
    rows; the regions and driver modules nominating the pair are kept as lists, and a pair is
    a driver pair if any region names it so.
    """
    d = flagged.copy()
    d["tx_a"] = d[["t1", "t2"]].min(axis=1)
    d["tx_b"] = d[["t1", "t2"]].max(axis=1)
    keys = ["gene", "tx_a", "tx_b"]
    sym = d.groupby(keys)[["pair_detected", "switch_like", "max_mean_if"]].nunique(dropna=False)
    if (sym > 1).any().any():
        raise ValueError("region-dependent long-read outcome for a distinct pair")
    join = lambda s: ";".join(sorted({x for v in s for x in str(v).split(";") if x}))  # noqa: E731
    return d.groupby(keys, as_index=False).agg(
        pair_detected=("pair_detected", "first"),
        switch_like=("switch_like", "first"),
        switch_like_strict=("switch_like_strict", "first"),
        usage_spearman=("usage_spearman", "first"),
        max_mean_if=("max_mean_if", "first"),
        driver_pair=("driver_pair", "any"),
        driver_usable=("driver_usable", "any"),
        driver_mean_if=("driver_mean_if", "max"),
        regions=("region", join),
        driver_regions=("region", lambda s: join(s[d.loc[s.index, "driver_pair"]])),
        driver_modules=("driver_modules", join),
        n_regions=("region", "nunique"),
    )


# --------------------------------------------------------------------------- #
# Summary
# --------------------------------------------------------------------------- #
def _rate(num: int, den: int) -> float:
    return float(num / den) if den else float("nan")


def _null_block(focal: pd.DataFrame, bg: pd.DataFrame, n_draws: int, seed: int) -> dict:
    tested = focal[focal["pair_detected"]]
    usable = tested[tested["driver_usable"]]
    return {
        "n_background_pairs": int(len(bg)),
        "n_background_detected": int(bg["pair_detected"].sum()),
        "background_detection_rate": _rate(int(bg["pair_detected"].sum()), len(bg)),
        "background_switch_like_rate_unmatched": _rate(
            int(bg.loc[bg["pair_detected"], "switch_like"].sum()), int(bg["pair_detected"].sum())),
        "switch_like_rate": matched_null(tested, bg[bg["pair_detected"]], "switch_like_rate", n_draws, seed),
        "mean_usage_spearman": matched_null(tested, bg[bg["pair_detected"]], "mean_usage_spearman", n_draws, seed),
        "switch_like_rate_usable_only": matched_null(usable, bg[bg["pair_detected"]], "switch_like_rate", n_draws, seed),
    }


def summarize(drivers: pd.DataFrame, flagged: pd.DataFrame, dp: pd.DataFrame, args) -> dict:
    focal = dp[dp["driver_pair"]]
    tested = focal[focal["pair_detected"]]
    usable = tested[tested["driver_usable"]]
    bg_all = dp[~dp["driver_pair"]]
    bg_gene = bg_all[bg_all["gene"].isin(set(focal["gene"]))]

    per_region = []
    for region, g in flagged.groupby("region"):
        f = g[g["driver_pair"]]
        t = f[f["pair_detected"]]
        dr = drivers[drivers["region"] == region]
        per_region.append({
            "region": region,
            "n_age_modules": int(dr["module_id"].nunique()),
            "n_driver_transcripts": int(dr["transcript_id"].nunique()),
            "n_drivers_in_coding_pairs": int(dr["in_coding_switch_pairs"].sum()),
            "n_drivers_detected_longread": int(dr["detected_longread"].sum()),
            "n_driver_pairs": int(len(f)),
            "n_driver_pairs_detected": int(len(t)),
            "driver_pair_detection_rate": _rate(len(t), len(f)),
            "n_driver_pairs_switch_like": int(t["switch_like"].sum()),
            "driver_switch_like_rate_of_detected": _rate(int(t["switch_like"].sum()), len(t)),
            "nondriver_pair_detection_rate": _rate(int(g.loc[~g["driver_pair"], "pair_detected"].sum()),
                                                   int((~g["driver_pair"]).sum())),
        })

    n_drv_univ = int(drivers.drop_duplicates("transcript_id")["in_coding_switch_pairs"].sum())
    n_drv_det = int(drivers.drop_duplicates("transcript_id")
                    .loc[lambda d: d["in_coding_switch_pairs"], "detected_longread"].sum())
    return {
        "analysis": "Long-read confirmation of IsoGraph module driver switches",
        "question": ("Do the top-k driver transcripts of the age-associated GTEx cortical co-switch "
                     "modules exist, and behave as switches, in orthogonal ONT long-read DLPFC data?"),
        "longread_dataset": _LONGREAD_DATASET,
        "regions": list(args.regions),
        "driver_rule": f"top-{args.k} transcripts by |r| per module (transcript_polarity_table); "
                       f"modules at age spline F-test FDR <= {args.alpha}",
        "switch_like_definition": "both members detected and within-gene usage Spearman rho < 0 "
                                  "(switch_orthogonal_confirm definition); *_strict additionally "
                                  "requires both mean IF >= min_if (longread_switch_confirm column)",
        "thresholds": {"k": args.k, "alpha": args.alpha, "min_if": args.min_if,
                       "n_draws": args.n_draws, "seed": args.seed, "variant": args.variant},
        "driver_transcripts": {
            "n_distinct": int(drivers["transcript_id"].nunique()),
            "n_in_coding_switch_pairs": n_drv_univ,
            "n_detected_longread": n_drv_det,
            "detection_rate_of_in_pairs": _rate(n_drv_det, n_drv_univ),
        },
        "driver_pairs": {
            "n_distinct": int(len(focal)),
            "n_in_both_regions": int((focal["n_regions"] > 1).sum()),
            "n_genes": int(focal["gene"].nunique()),
            "n_detected": int(len(tested)),
            "detection_rate": _rate(len(tested), len(focal)),
            "n_switch_like": int(tested["switch_like"].sum()),
            "switch_like_rate_of_detected": _rate(int(tested["switch_like"].sum()), len(tested)),
            "n_switch_like_strict": int(tested["switch_like_strict"].sum()),
            "switch_like_strict_rate_of_detected": _rate(int(tested["switch_like_strict"].sum()), len(tested)),
            "n_detected_usable": int(len(usable)),
            "n_switch_like_usable": int(usable["switch_like"].sum()),
            "switch_like_rate_usable_only": _rate(int(usable["switch_like"].sum()), len(usable)),
            "n_genes_with_detected_pair": int(tested["gene"].nunique()),
            "n_genes_switch_like": int(tested.loc[tested["switch_like"], "gene"].nunique()),
            "n_genes_switch_like_usable": int(usable.loc[usable["switch_like"], "gene"].nunique()),
        },
        "matched_null": {
            "nondriver_switch_pairs": _null_block(focal, bg_all, args.n_draws, args.seed),
            "nondriver_pairs_same_genes": _null_block(focal, bg_gene, args.n_draws, args.seed),
        },
        "per_region": per_region,
    }


def write_markdown(summary: dict, out_dir: Path) -> None:
    th, dt, dpr = summary["thresholds"], summary["driver_transcripts"], summary["driver_pairs"]
    lines = [
        "# Long-read confirmation of module driver switches",
        "",
        f"Generated by `longread_driver_confirm.py` (k = {th['k']}, alpha = {th['alpha']}, "
        f"min IF = {th['min_if']}, {th['n_draws']} draws, seed {th['seed']}).",
        "",
        "## Question",
        "",
        summary["question"],
        "",
        f"Data: {summary['longread_dataset']}. Regions: {', '.join(summary['regions'])}. "
        f"Driver rule: {summary['driver_rule']}.",
        "",
        "## Existence",
        "",
        f"- Driver transcripts: {dt['n_distinct']} distinct; {dt['n_in_coding_switch_pairs']} are members "
        f"of a coding switch pair, of which **{dt['n_detected_longread']} "
        f"({100 * dt['detection_rate_of_in_pairs']:.1f}%)** are detected in long read.",
        f"- Distinct driver pairs: {dpr['n_distinct']} in {dpr['n_genes']} genes "
        f"({dpr['n_in_both_regions']} nominated in both regions); **{dpr['n_detected']} "
        f"({100 * dpr['detection_rate']:.1f}%)** have both members detected, against "
        f"{100 * summary['matched_null']['nondriver_switch_pairs']['background_detection_rate']:.1f}% "
        "of non-driver switch pairs.",
        "",
        "## Switch behaviour (detected-pair denominator, abundance-matched null)",
        "",
    ]
    rows = []
    for bg_name, label in (("nondriver_switch_pairs", "all non-driver switch pairs"),
                           ("nondriver_pairs_same_genes", "non-driver pairs, driver genes")):
        blk = summary["matched_null"][bg_name]
        for stat, slabel in (("switch_like_rate", "switch-like rate"),
                             ("switch_like_rate_usable_only", "switch-like rate, driver IF usable"),
                             ("mean_usage_spearman", "mean usage rho")):
            m = blk[stat]
            rows.append({"background": label, "statistic": slabel, "n focal": m.get("n_focal"),
                         "observed": m.get("observed"), "null mean": m.get("null_mean"),
                         "null 95%": (f"{m['null_q025']:.3f}-{m['null_q975']:.3f}"
                                      if "null_q025" in m else ""),
                         "p": m.get("p_empirical_two_sided")})
    lines += [_md_table(pd.DataFrame(rows), floats=3), ""]
    lines += [
        f"- Strict (both IF >= {th['min_if']}) switch-like rate of detected driver pairs: "
        f"{dpr['n_switch_like_strict']}/{dpr['n_detected']} = "
        f"{100 * dpr['switch_like_strict_rate_of_detected']:.1f}%.",
        f"- Genes: {dpr['n_genes_with_detected_pair']} with a detected driver pair; "
        f"{dpr['n_genes_switch_like']} with a switch-like one; "
        f"{dpr['n_genes_switch_like_usable']} at usable driver abundance.",
        "",
        "## Per region (descriptive; pairs shared by both regions appear in both rows)",
        "",
        _md_table(pd.DataFrame(summary["per_region"]), floats=3),
        "",
        "## Reading guide",
        "",
        "- A bare negative usage correlation is near-vacuous evidence of a switch (compositional "
        "closure; see `../switch_orthogonal_confirm/GLOBAL_NULL.md`), so the rate must be read "
        "against the matched null, never against zero.",
        "- The two backgrounds answer different questions: all non-driver switch pairs is the "
        "anchored arm's universe and the primary comparison; the same-gene background controls for "
        "gene but is small.",
        "- n = 12 donors, half AD, DLPFC against GTEx cortical regions; no per-sample age, so "
        "nothing here is a signed usage-vs-age replication.",
        "",
    ]
    (out_dir / "LONGREAD_DRIVER_CONFIRM.md").write_text("\n".join(lines))


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def run(args) -> None:
    out_dir = ensure_dir(out_path(args.k))
    pairs = load_switch_pairs()
    pairs = pairs[pairs["region"].isin(args.regions)].copy()
    if pairs.empty:
        raise SystemExit(f"no scored switch pairs for regions {args.regions}")

    drivers = pd.concat(
        [driver_transcripts(r, args.k, args.alpha, args.variant) for r in args.regions],
        ignore_index=True,
    )
    if drivers.empty:
        raise SystemExit("no driver transcripts: no age-associated modules or polarity tables")

    # transcript-level existence, from either pair slot
    tx = pd.concat([pairs[["region", "t1", "t1_detected"]].set_axis(["region", "tx", "det"], axis=1),
                    pairs[["region", "t2", "t2_detected"]].set_axis(["region", "tx", "det"], axis=1)])
    det = tx.groupby(["region", "tx"])["det"].any()
    key = pd.MultiIndex.from_arrays([drivers["region"], drivers["transcript_id"]])
    drivers["in_coding_switch_pairs"] = key.isin(det.index)
    drivers["detected_longread"] = det.reindex(key).fillna(False).to_numpy(dtype=bool)

    flagged = flag_driver_pairs(pairs, drivers, args.min_if)
    dp = distinct_pairs(flagged)
    summary = summarize(drivers, flagged, dp, args)

    drivers.to_parquet(out_dir / "driver_transcripts.parquet", index=False)
    dp[dp["driver_pair"]].drop(columns=["driver_pair"]).to_parquet(
        out_dir / "driver_pair_confirmation.parquet", index=False)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_markdown(summary, out_dir)

    dpr, mn = summary["driver_pairs"], summary["matched_null"]["nondriver_switch_pairs"]["switch_like_rate"]
    print(f"driver pairs: {dpr['n_distinct']} distinct, {dpr['n_detected']} detected "
          f"({100 * dpr['detection_rate']:.1f}%), switch-like {dpr['n_switch_like']}/{dpr['n_detected']} "
          f"= {dpr['switch_like_rate_of_detected']:.3f} vs matched null {mn.get('null_mean', float('nan')):.3f} "
          f"(p = {mn.get('p_empirical_two_sided', float('nan')):.3g})", flush=True)
    print(f"wrote {out_dir}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--k", type=int, default=DEFAULT_K, help="top-k driver transcripts per module")
    ap.add_argument("--alpha", type=float, default=0.05, help="module age FDR (spline F-test)")
    ap.add_argument("--min-if", type=float, default=0.05, help="usable driver mean isoform fraction")
    ap.add_argument("--regions", nargs="+", default=list(_DEFAULT_REGIONS))
    ap.add_argument("--variant", default="standard", choices=["standard", "with-abundance"])
    ap.add_argument("--n-draws", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=13)
    run(ap.parse_args())


if __name__ == "__main__":
    main()
