"""Display-ready tables for the module-reproducibility results subsection.

Presentation layer only: this reads the committed ledgers written by `module_trust.py`,
`eigengene_projection.py`, `replication_permutation.py` and `replication_functional.py`
and re-emits them as flat CSVs. It fits no model and computes no new statistic beyond
counting and taking medians of columns that are already there -- every number it writes
can be traced to one column of one parquet named in `SOURCES` below.

Two reasons it exists rather than having the figure script read the parquets directly:

1. The figure must build on a laptop. Several ledgers are zstd-compressed and a local
   `arrow` build without zstd cannot open them, so the canonical parquet stays the source
   of record and the figure reads the CSV mirror (the same split used by
   `composition_age_coupling.py` and `switch_unique_threshold.py`).
2. The funnel claims in the figure's summary panel are the SAME four fractions the
   Results text quotes. Computing them once, here, keeps the figure, the supplementary
   tables and the prose from drifting apart -- which is exactly what happened to the
   driver-loading rho range when it was re-typed from a summary written at a superseded
   Leiden resolution.

Writes 04_module_trust/_m/stability/module_trust_tables/:
  funnel_claims.csv           the four reproducibility claims, both methods (Fig 3a)
  split_half_modules.csv      per-module split-half trust ledger, both methods (Fig 3b)
  split_half_pairs.csv        per split-half module pair, both methods (Fig 3c,d)
  projection_modules.csv      per-module frozen-eigengene projection ledger (Fig 3e)
  projection_summary.csv      the four method-by-direction projection rows
  functional_preservation.csv matched-vs-null GO and cell-type similarity (Fig 3f)
  crosscohort_permutation.csv the full permutation grid behind the matched-pair count
  resolution_sensitivity.csv  split-half agreement across the Leiden resolution sweep
  projection_sign_scale.csv   raw vs null-standardised projected-age sign agreement
  region_funnel.csv           per-region funnel, both methods (Table S7)
  driver_structure.csv        structural classes of driver switches (moved out of Fig 3)
  MODULE_TRUST_TABLES.md      provenance and the one-line read of each file

Run: bash 04_module_trust/_h/04e.module_trust_tables.sh
     python -m isograph_benchmark.real_data.module_trust_tables [--out DIR]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from isograph_benchmark.paths import ensure_dir, stage_out

TRUST = stage_out("trust.stability", "module_trust")
PROJ = stage_out("trust.stability", "eigengene_projection")
STAB_SUMMARY = stage_out("trust.stability") / "stability_summary.parquet"
FUNC = stage_out("integration", "functional_preservation")
OUT_DEFAULT = stage_out("trust.stability", "module_trust_tables")

# The published age model. `spline` exists for both methods and is carried in the
# permutation grid, but the projection and split-half arms are reported on the linear
# arm, so the figure shows that one and says so.
PUBLISHED_MODEL = "linear"
# The statistic the Results text quotes for the matched-pair count. The stage's
# pre-registered primary is `spline_f`; both are written, and `MODULE_TRUST_TABLES.md`
# names which is which so neither can be quoted alone by accident.
PUBLISHED_STATISTIC = "pearson"
# Production Leiden resolution (PI decision, 2026-09-16). Artifacts written at 5.0 are
# stale; the sweep arms below exist so the stability claim is not read as a property of
# this one choice.
PRODUCTION_RESOLUTION = 2.0

SOURCES = {
    "split_half_modules": "module_stability__<cohort>__<region>__<method>.parquet",
    "split_half_pairs": "within_cohort__<cohort>__<region>__<method>.parquet",
    "projection_modules": "eigengene_projection/eigengene_projection_all.parquet",
    "projection_summary": "eigengene_projection/eigengene_projection_summary.parquet",
    "functional_preservation": "functional_preservation__<method>__<model>__stats.json",
    "crosscohort_permutation": "replication_permutation__*__stats.json",
    "driver_structure": "module_complementarity__<cohort>__<region>__<method>.parquet",
    "resolution_sensitivity": "stability_summary.parquet",
    "projection_sign_scale": "eigengene_projection/eigengene_projection_all.parquet + "
                             "eigengene_projection/switch_axis_alignment__<region>.parquet",
}


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def _method_of(path: Path) -> str:
    """Method is encoded in the filename suffix, not in a column, for these ledgers."""
    return "wgcna" if path.stem.endswith("__wgcna") else "isograph"


def load_stage(prefix: str, trust_dir: Path = TRUST) -> pd.DataFrame:
    """Concatenate every per-region ledger for one funnel stage.

    Pooled tables are excluded: they are a different unit (one row per method over all
    region pairs) and would double-count modules if concatenated with the per-region ones.
    """
    files = sorted(
        f for f in trust_dir.glob(f"{prefix}*.parquet") if "_pooled__" not in f.name
    )
    if not files:
        raise SystemExit(f"no {prefix}*.parquet under {trust_dir}")
    frames = []
    for f in files:
        df = pd.read_parquet(f)
        if "method" not in df.columns:
            df = df.assign(method=_method_of(f))
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------------------- #
# Derived tables (pure functions of the loaded frames, so they are testable)
# --------------------------------------------------------------------------- #
def funnel_claims(stab: pd.DataFrame, pairs: pd.DataFrame,
                  proj_summary: pd.DataFrame) -> pd.DataFrame:
    """The four reproducibility claims on one fraction scale, per method.

    The denominators differ by claim and that is deliberate, so each row carries its own:
    modules for the trust claim, both-significant split-half pairs for the sign claim, and
    age-testable modules for each projection direction. A fraction without its denominator
    is not quotable here -- the split-half sign claim rests on 23 pairs out of 641.
    """
    rows = []
    for method, d in stab.groupby("method", sort=True):
        rows.append(dict(claim="modules_chance_trusted", method=method,
                         num=int(d["trusted"].sum()), den=int(len(d))))
    for method, d in pairs.groupby("method", sort=True):
        both = d[d["both_age_sig"].astype(bool)]
        rows.append(dict(claim="split_half_age_sign_concordance", method=method,
                         num=int(both["sign_concordant"].astype(bool).sum()),
                         den=int(len(both))))
    for _, r in proj_summary.iterrows():
        rows.append(dict(claim=f"aging_axis_transfers__{r['direction']}",
                         method=str(r["method"]),
                         num=int(r["sign_match"]), den=int(r["n_age_testable"])))
    df = pd.DataFrame(rows)
    df["value"] = df["num"] / df["den"].where(df["den"] > 0)
    # Claim order is the order the Results subsection makes the claims in.
    order = ["modules_chance_trusted", "split_half_age_sign_concordance",
             "aging_axis_transfers__brainseq_to_gtex",
             "aging_axis_transfers__gtex_to_brainseq"]
    df["claim"] = pd.Categorical(df["claim"], order, ordered=True)
    return df.sort_values(["claim", "method"]).reset_index(drop=True)


def region_funnel(stab: pd.DataFrame, pairs: pd.DataFrame,
                  comp: pd.DataFrame | None) -> pd.DataFrame:
    """Per cohort x region x method funnel row (Table S7).

    Both methods, unlike the IsoGraph-only table this replaces: the whole subsection is a
    matched contrast, so a reader cannot check it from IsoGraph rows alone. The
    complementarity columns stay IsoGraph-only because WGCNA has no transcript drivers.
    """
    keys = ["cohort", "region", "method"]
    rows = []
    for (cohort, region, method), s in stab.groupby(keys, sort=True):
        w = pairs[(pairs["cohort"] == cohort) & (pairs["region"] == region)
                  & (pairs["method"] == method)]
        both = w[w["both_age_sig"].astype(bool)]
        rho = w["driver_load_rho"].dropna() if "driver_load_rho" in w else pd.Series(dtype=float)
        row = dict(
            cohort=cohort, region=region, method=method,
            n_modules=int(len(s)), n_trusted=int(s["trusted"].sum()),
            frac_trusted=float(s["trusted"].mean()),
            median_coassign_density=float(s["coassign_density"].median()),
            median_null_density=float(s["null_mean"].median()),
            median_best_match_jaccard=float(s["best_match_jaccard"].median()),
            n_split_half_pairs=int(len(w)),
            n_both_age_sig=int(len(both)),
            n_sign_concordant=int(both["sign_concordant"].astype(bool).sum()),
            median_driver_rho=float(rho.median()) if len(rho) else None,
            frac_positive_rho=float((rho > 0).mean()) if len(rho) else None,
        )
        if comp is not None and method == "isograph":
            # Filter on method as well as region: the complementarity ledgers carry no
            # method column of their own, so an unfiltered slice would take the median
            # over the IsoGraph and WGCNA rows together.
            c = comp[(comp["cohort"] == cohort) & (comp["region"] == region)
                     & (comp["method"] == method)]
            if len(c):
                row["median_frac_dtu_without_dge"] = float(
                    c["frac_dtu_without_dge"].median())
                row["median_frac_in_wgcna_age"] = float(
                    c["frac_in_wgcna_age_modules"].median())
        rows.append(row)
    return pd.DataFrame(rows)


def functional_preservation_rows(func_dir: Path = FUNC) -> pd.DataFrame:
    """Matched-pair similarity against the size-matched null, per measure.

    A measure with `n_finite == 0` is written with null statistics and kept in the table:
    it is an absent upstream input, not a null result, and dropping the row silently would
    turn "untested" into "tested and negative".
    """
    rows = []
    for f in sorted(func_dir.glob("functional_preservation__*__stats.json")):
        st = json.loads(f.read_text())
        for measure, m in st["measures"].items():
            rows.append(dict(
                method=st["method"], model=st["model"], measure=measure,
                n_pairs=st["n_pairs"], n_concordant=st["n_concordant"],
                median_gene_jaccard=st["median_gene_jaccard"],
                n_finite=m["n_finite"], mean_matched=m["mean_matched"],
                null_mean=m["null_mean"], null_sd=m["null_sd"], p_emp=m["p_emp"],
                mean_concordant=m["mean_concordant"],
                mean_discordant=m["mean_discordant"],
                n_permutations=st["n_permutations"], seed=st["seed"],
            ))
    return pd.DataFrame(rows).sort_values(
        ["model", "method", "measure"]).reset_index(drop=True)


def permutation_grid(trust_dir: Path = TRUST) -> pd.DataFrame:
    """The full covariate-mode x statistic x null grid behind the matched-pair count.

    All three statistics are written, flagged, because the covariate-free Pearson arm the
    text quotes (1/38) and the stage's pre-registered covariate-adjusted spline arm (0/38)
    do not agree, and the decision rule in REPLICATION_PERMUTATION.md requires both to be
    reported rather than the more favourable one alone.
    """
    rows = []
    for f in sorted(trust_dir.glob("replication_permutation__*__stats.json")):
        st = json.loads(f.read_text())
        rows.append(dict(
            method=st["method"], statistic=st["statistic"], null=st["null"],
            covariates=st["covariates"], n_pairs=st["n_matched_rows"],
            t_obs=st["T_obs"], null_mean=st["null_mean"], null_sd=st["null_sd"],
            null_q95=st["null_q95"], z=st["z"], p_emp=st["p_emp"],
            n_permutations=st["B"], seed=st["seed"],
            is_published_statistic=st["statistic"] == PUBLISHED_STATISTIC,
        ))
    df = pd.DataFrame(rows)
    return df.sort_values(["method", "statistic", "null", "covariates"]).reset_index(drop=True)


def driver_structure(comp: pd.DataFrame) -> pd.DataFrame:
    """Structural classes of the top driver switches, age-significant IsoGraph modules.

    This was a panel of the trust figure until the reproducibility subsection took the
    figure over; the structural classes belong with the transcript-evidence results, so
    the numbers move to their own display item rather than being dropped.
    """
    cols = {"drv_cds_changed": "CDS change", "drv_utr_changed": "UTR change",
            "drv_biotype_switch": "Biotype switch",
            "drv_coding_status_change": "Coding-status change"}
    d = comp[(comp["method"] == "isograph") & comp["age_sig"].astype(bool)]
    rows = []
    for col, label in cols.items():
        v = d[col].dropna() if col in d else pd.Series(dtype=float)
        rows.append(dict(switch_class=label, column=col, n_modules=int(len(v)),
                         mean_frac=float(v.mean()) if len(v) else None,
                         se_frac=float(v.std(ddof=1) / len(v) ** 0.5) if len(v) > 1 else None,
                         median_frac=float(v.median()) if len(v) else None))
    return pd.DataFrame(rows)


def _resolution_of(method: str) -> float | None:
    """Leiden resolution encoded in the method label; None for the fixed WGCNA baseline.

    `isograph` is the production run at 2.0 (PI decision, 2026-09-16) and the sweep arms
    are `isograph_res0p5` ... `isograph_res20`, where `p` stands in for the decimal point.
    """
    if not method.startswith("isograph"):
        return None
    if "_res" not in method:
        return PRODUCTION_RESOLUTION
    return float(method.split("_res", 1)[1].replace("p", "."))


def resolution_sensitivity(summary: pd.DataFrame) -> pd.DataFrame:
    """Split-half agreement across the Leiden resolution sweep, per region.

    The subsection claims stability across donor subsets at one resolution, so the sweep
    is what shows that claim is not an artefact of that choice. WGCNA rows are kept with a
    null resolution rather than dropped: it is the fixed baseline the whole subsection is
    read against, and a sweep table without it invites the reader to compare IsoGraph
    resolutions with each other only.
    """
    df = summary.copy()
    df["resolution"] = df["method"].map(_resolution_of)
    df["method_family"] = df["method"].where(
        ~df["method"].astype(str).str.startswith("isograph"), "isograph")
    df["is_production"] = df["resolution"].eq(PRODUCTION_RESOLUTION) | df["resolution"].isna()
    cols = ["method_family", "method", "resolution", "is_production", "cohort", "region",
            "comparison", "n_seeds", "mean_ari", "sd_ari", "mean_nmi", "sd_nmi",
            "mean_n_common"]
    cols = [c for c in cols if c in df.columns]
    return df[cols].sort_values(
        ["method_family", "resolution", "cohort", "region"],
        na_position="last").reset_index(drop=True)


def switch_axis_orientation(proj_dir: Path = PROJ) -> pd.DataFrame:
    """Per-region sign convention of the shared-gene switch axis.

    The projection carries a source module's weights into the target cohort unchanged, so
    a gene whose switch axis is fitted with the opposite sign in the two cohorts enters
    the projection flipped. About half of the orientable genes are, which is the reason
    the age statistic is null-standardised rather than compared as a raw coefficient.

    Two medians are written because they answer different questions: `median_abs_cosine`
    over every shared gene (the value quoted in `EIGENGENE_PROJECTION.md`) mixes in genes
    with too few shared transcripts to orient, and `median_abs_cosine_orientable` is the
    alignment among genes that can actually be oriented.
    """
    rows = []
    for f in sorted(proj_dir.glob("switch_axis_alignment__*.parquet")):
        a = pd.read_parquet(f)
        usable = a[a["usable"].astype(bool)]
        rows.append(dict(
            pair=f.stem.split("__", 1)[1],
            n_shared_genes=int(len(a)),
            n_orientable=int(len(usable)),
            n_flipped=int((usable["sign"] < 0).sum()),
            frac_flipped=float((usable["sign"] < 0).mean()) if len(usable) else None,
            median_abs_cosine=float(a["cosine"].abs().median()) if len(a) else None,
            median_abs_cosine_orientable=(
                float(usable["cosine"].abs().median()) if len(usable) else None),
        ))
    return pd.DataFrame(rows)


def projection_sign_scale(proj_all: pd.DataFrame,
                          orientation: pd.DataFrame | None = None) -> pd.DataFrame:
    """Raw against null-standardised projected-age sign agreement, per region.

    The text quotes the standardised counts, and this is the table that says why. A raw
    projected age correlation inherits whatever age-correlated structure the target cohort
    carries -- RNA quality, ischemic time, composition -- so raw sign agreement is close
    to unanimous within a region and can be unanimous in either direction: among the
    modules significant in both cohorts on raw correlations, two of the three
    BrainSEQ-to-GTEx regions agree completely and the third disagrees completely. That is
    a property of the target cohort, not of the modules.

    `raw_sign_match_all` counts every module and `raw_sign_match_both_sig` only those
    detectable in both cohorts; both are written, because the two denominators are what
    make the raw and standardised arms look different.
    """
    d = proj_all.copy()
    for col in ("raw_sign_match", "raw_both_sig", "sign_match", "both_sig"):
        d[col] = d[col].astype(bool)
    rows = []
    for (method, direction, pair), g in d.groupby(["method", "direction", "pair"],
                                                  sort=True):
        raw_both = g[g["raw_both_sig"]]
        std_both = g[g["both_sig"]]
        rows.append(dict(
            method=method, direction=direction, pair=pair,
            n_modules=int(len(g)),
            raw_n_both_sig=int(len(raw_both)),
            raw_sign_match_both_sig=int(raw_both["raw_sign_match"].sum()),
            raw_sign_match_all=int(g["raw_sign_match"].sum()),
            std_n_both_sig=int(len(std_both)),
            std_sign_match_all=int(g["sign_match"].sum()),
            median_age_r_target=float(g["age_r_target"].median()),
            median_age_z_target=float(g["age_z_target"].median()),
        ))
    out = pd.DataFrame(rows)
    if orientation is not None and len(orientation):
        # The switch axis is an IsoGraph object: WGCNA projects abundance features, which
        # have no orientation to lose, so those rows stay empty instead of borrowing
        # IsoGraph's numbers.
        out = out.merge(orientation.assign(method="isograph"),
                        on=["method", "pair"], how="left")
        assert out.loc[out["method"] != "isograph", "n_orientable"].isna().all()
    return out.sort_values(["method", "direction", "pair"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def _fmt_frac(num: int, den: int) -> str:
    return f"{num}/{den}" + (f" ({num / den:.2f})" if den else "")


def write_report(out: Path, claims: pd.DataFrame, func: pd.DataFrame,
                 perm: pd.DataFrame, sign_scale: pd.DataFrame | None = None) -> None:
    lines = [
        "# Module trust: display tables",
        "",
        "Written by `isograph_benchmark/real_data/module_trust_tables.py`. Presentation "
        "only -- every column is copied from the ledger named beside it; regenerate "
        "rather than editing.",
        "",
        "| file | source ledger |",
        "|---|---|",
    ]
    for name, src in SOURCES.items():
        lines.append(f"| `{name}.csv` | `{src}` |")

    lines += ["", "## The four claims", "",
              "| claim | IsoGraph | WGCNA |", "|---|---|---|"]
    for claim, d in claims.groupby("claim", observed=True, sort=True):
        cell = {}
        for _, r in d.iterrows():
            cell[r["method"]] = _fmt_frac(int(r["num"]), int(r["den"]))
        lines.append(f"| {claim} | {cell.get('isograph', 'n/a')} | "
                     f"{cell.get('wgcna', 'n/a')} |")

    pub = func[func["model"] == PUBLISHED_MODEL]
    lines += ["", f"## Functional preservation of matched pairs (`{PUBLISHED_MODEL}` arm)",
              "",
              "| method | measure | n | matched mean | size-matched null | p_emp |",
              "|---|---|---|---|---|---|"]
    for _, r in pub.iterrows():
        if not r["n_finite"]:
            lines.append(f"| {r['method']} | `{r['measure']}` | 0 | untested "
                         "(no computable pair) | n/a | n/a |")
            continue
        lines.append(f"| {r['method']} | `{r['measure']}` | {int(r['n_finite'])} | "
                     f"{r['mean_matched']:.4f} | {r['null_mean']:.4f} | "
                     f"{r['p_emp']:.3g} |")

    lines += ["", "## Matched-pair count: the two arms disagree", ""]
    strict = perm[perm["null"] == "matching"]
    for method, d in strict.groupby("method", sort=True):
        for stat in ("pearson", "spline_f"):
            s = d[(d["statistic"] == stat) & (d["covariates"] == "none")]
            if not len(s):
                continue
            r = s.iloc[0]
            tag = " *(quoted in the text)*" if stat == PUBLISHED_STATISTIC else \
                  " *(stage's pre-registered primary)*"
            lines.append(f"- **{method}**, `{stat}` vs the matching null: "
                         f"{int(r['t_obs'])}/{int(r['n_pairs'])}, "
                         f"p_emp = {r['p_emp']:.3g}{tag}")
    if sign_scale is not None and len(sign_scale):
        lines += ["", "## Why the projected-age statistic is null-standardised", "",
                  "Raw sign agreement is counted among the modules significant in both "
                  "cohorts on raw correlations; the standardised column counts every "
                  "age-testable module. A region that agrees unanimously in one direction "
                  "and a region that disagrees unanimously are the same phenomenon: the "
                  "target cohort's own age-correlated structure enters every projection.",
                  "",
                  "| method | direction | region | raw (both-sig) | standardised (all) |",
                  "|---|---|---|---|---|"]
        for _, r in sign_scale.iterrows():
            lines.append(
                f"| {r['method']} | {r['direction']} | {r['pair']} | "
                f"{_fmt_frac(int(r['raw_sign_match_both_sig']), int(r['raw_n_both_sig']))} | "
                f"{_fmt_frac(int(r['std_sign_match_all']), int(r['n_modules']))} |")

    lines += ["",
              "Report both. The pre-registered rule in `REPLICATION_PERMUTATION.md` "
              "forbids the word *replication* for this arm; the honest wording is "
              "\"matched modules with concordant age effects\".",
              ""]
    (out / "MODULE_TRUST_TABLES.md").write_text("\n".join(lines))


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT,
                    help="output directory (default: the stage's module_trust_tables/)")
    args = ap.parse_args(argv)
    out = ensure_dir(args.out)

    stab = load_stage("module_stability__")
    pairs = load_stage("within_cohort__")
    try:
        comp = load_stage("module_complementarity__")
    except SystemExit:
        comp = None

    proj_all = pd.read_parquet(PROJ / "eigengene_projection_all.parquet")
    proj_sum = pd.read_parquet(PROJ / "eigengene_projection_summary.parquet")

    claims = funnel_claims(stab, pairs, proj_sum)
    sign_scale = projection_sign_scale(proj_all, switch_axis_orientation())
    func = functional_preservation_rows()
    perm = permutation_grid()

    written: list[tuple[str, pd.DataFrame]] = [
        ("funnel_claims", claims),
        ("split_half_modules", stab),
        ("split_half_pairs", pairs),
        ("projection_modules", proj_all),
        ("projection_summary", proj_sum),
        ("functional_preservation", func),
        ("crosscohort_permutation", perm),
        ("region_funnel", region_funnel(stab, pairs, comp)),
        ("projection_sign_scale", sign_scale),
    ]
    if STAB_SUMMARY.exists():
        written.append(("resolution_sensitivity",
                        resolution_sensitivity(pd.read_parquet(STAB_SUMMARY))))
    if comp is not None:
        written.append(("driver_structure", driver_structure(comp)))

    print(f"Writing module-trust display tables to {out}")
    for name, df in written:
        df.to_csv(out / f"{name}.csv", index=False)
        print(f"  wrote {name + '.csv':30s} {df.shape[0]:>4d} x {df.shape[1]}")

    write_report(out, claims, func, perm, sign_scale)
    print("  wrote MODULE_TRUST_TABLES.md")

    for claim, d in claims.groupby("claim", observed=True, sort=True):
        parts = " ".join(f"{r['method']} {int(r['num'])}/{int(r['den'])}"
                         for _, r in d.iterrows())
        print(f"    {claim:42s} {parts}")


if __name__ == "__main__":
    main()
