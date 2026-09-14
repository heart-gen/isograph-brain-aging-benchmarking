"""Transcript-filter arms: how BrainSEQ modules and module detection change under a switching filter.

WHY THIS EXISTS
---------------
The production BrainSEQ aging fits keep transcripts with count > 10 in >= 70% of samples
(`run_models._filter_expressed_transcripts`). That is an abundance, differential-expression
style prefilter. The production SCZD caudate fit applies no transcript filter at all. IsoGraph
models isoform *usage*, and the standard prefilter for usage (Soneson et al. 2016, Genome
Biology 17:12; DRIMSeq `dmFilter`, Nowicka & Robinson 2016) asks for an expressed gene and for
each transcript to be expressed and a non-trivial share of its gene in *some* samples, not in
most. A 70% prevalence requirement structurally removes isoforms used in part of the cohort,
which is what an on/off switch looks like.

Each arm re-runs the production fit function (`run_models.run_brainseq_region` /
`run_brainseq_caudate_sczd`: same VAE, covariates, calibration grid, Leiden resolution, seed)
with only the transcript filter changed, writes the full production artifact set OUTSIDE the
production stores, and `compare` reads the committed fits beside the arms. Nothing under
`<store>/isograph_vae/` is written.

ARMS
----
abundance   count > 10 in >= 70% of samples (the production aging filter). For the aging regions
            this refit is the noise floor; for caudate_sczd it is the BrainSEQ method applied to
            the SCZD fit.
switching   the switching filter below.
none        no transcript filter: the SCZD production setting, refit (its noise floor).

SWITCHING FILTER (phenotype-blind; DRIMSeq semantics with sample fractions for group sizes)
    gene        total count >= 10 in >= 70% of samples, so its proportions are estimable
    transcript  count >= 10 in >= 10% of samples, AND
                share of its gene >= 0.10 in >= 10% of samples (share over samples where the
                gene count is > 0)
Deviation from DRIMSeq, on purpose: genes left with a single transcript are kept, as the
production filter keeps them, so the arms differ in which isoforms represent a gene and not
additionally in which genes carry an abundance channel. The 10% prevalence is not derived
from any trait, so one rule serves age and diagnosis.

STAGES
------
    --list               print the fit grid
    fit --index N        one fit per process (a real-data fit peaks near 30 GB)
    compare              per-fit and pairwise summaries + TRANSCRIPT_FILTER_ARMS.md

Outputs: 02_module_discovery/_m/transcript_filter_arms/brainseq/<region>/<arm>/ and, at the
root, fits.parquet, comparisons.parquet, TRANSCRIPT_FILTER_ARMS.md.
"""
from __future__ import annotations

import argparse
import itertools
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from isograph_benchmark.paths import ensure_dir, region_store, stage_out

COHORT = "brainseq"
FDR = 0.05
GIANT_GENES = 900  # the phenotype-blind giant-module criterion behind resolution 5.0

SWITCHING = dict(min_gene_count=10.0, min_gene_fraction=0.70,
                 min_tx_count=10.0, min_tx_prop=0.10, min_tx_fraction=0.10)

# region -> analysis; the analysis fixes the trait table and the production arm
TARGETS = {"caudate": "aging", "hippocampus": "aging", "dlpfc": "aging", "caudate_sczd": "sczd"}
PRODUCTION_ARM = {"aging": "abundance", "sczd": "none"}
TRAIT_FILE = {"aging": "age_linear.parquet", "sczd": "diagnosis_assoc.parquet"}

# Switching arms first, so a partial array still answers the question.
FITS: list[tuple[str, str]] = [
    ("caudate", "switching"), ("hippocampus", "switching"), ("dlpfc", "switching"),
    ("caudate_sczd", "switching"), ("caudate_sczd", "abundance"),
    ("caudate", "abundance"), ("hippocampus", "abundance"), ("dlpfc", "abundance"),
    ("caudate_sczd", "none"),
]

# Seed floor: every arm refit at more seeds. `random_state` drives the VAE initialisation, the
# minibatch order and the Leiden seed, so this is end-to-end seed variation; a same-seed refit is
# deterministic (identical edges) and cannot serve as a floor.
PRODUCTION_SEED = 13
EXTRA_SEEDS = (14, 15)
SEED_FITS: list[tuple[str, str, int]] = [(r, a, s) for s in EXTRA_SEEDS for r, a in FITS]
GRIDS: dict[str, list[tuple[str, str, int]]] = {
    "filter": [(r, a, PRODUCTION_SEED) for r, a in FITS],
    "seed": SEED_FITS,
}


# --------------------------------------------------------------------------- #
# Filters
# --------------------------------------------------------------------------- #
def switching_filter(counts, table: pd.DataFrame, min_gene_count: float = 10.0,
                     min_gene_fraction: float = 0.70, min_tx_count: float = 10.0,
                     min_tx_prop: float = 0.10, min_tx_fraction: float = 0.10):
    """Usage-oriented transcript filter (see module docstring). Returns (counts, table)."""
    c = np.asarray(counts)
    n = c.shape[1]
    gidx = pd.factorize(table["gene_id"].astype(str))[0]
    gene_tot = pd.DataFrame(c, copy=False).groupby(gidx).sum().sort_index().to_numpy(np.float64)
    gene_ok = (gene_tot >= min_gene_count).sum(axis=1) >= np.ceil(min_gene_fraction * n)
    need = np.ceil(min_tx_fraction * n)
    tx_expr = (c >= min_tx_count).sum(axis=1) >= need
    gt = gene_tot[gidx]
    with np.errstate(divide="ignore", invalid="ignore"):
        share_ok = np.where(gt > 0, c / gt, 0.0) >= min_tx_prop
    del gt
    keep = gene_ok[gidx] & tx_expr & (share_ok.sum(axis=1) >= need)
    kept = table.loc[keep].reset_index(drop=True)
    print(f"  switching filter (gene>={min_gene_count:g} in >={min_gene_fraction:.0%}; tx>="
          f"{min_tx_count:g} and share>={min_tx_prop:g} in >={min_tx_fraction:.0%} of {n}): "
          f"{int(keep.sum())}/{len(table)} transcripts, {kept['gene_id'].nunique()}/"
          f"{table['gene_id'].nunique()} genes retained", flush=True)
    return c[keep], kept


def _filters() -> dict:
    from isograph_benchmark.real_data.run_models import _filter_expressed_transcripts
    return {"abundance": _filter_expressed_transcripts, "switching": switching_filter, "none": None}


def arm_dir(region: str, arm: str, seed: int = PRODUCTION_SEED) -> Path:
    leaf = arm if seed == PRODUCTION_SEED else f"{arm}_seed{seed}"
    return stage_out("modules.filter_arms", COHORT, region, leaf)


def _tx_id_col(table: pd.DataFrame) -> str:
    return "transcript_id" if "transcript_id" in table.columns else table.columns[0]


# --------------------------------------------------------------------------- #
# Fit
# --------------------------------------------------------------------------- #
def fit(index: int, grid: str = "filter") -> None:
    from isograph_benchmark.real_data.run_models import (
        run_brainseq_caudate_sczd, run_brainseq_region,
    )

    fits = GRIDS[grid]
    if not 1 <= index <= len(fits):
        raise SystemExit(f"--index must be 1..{len(fits)} for grid {grid!r}")
    region, arm, seed = fits[index - 1]
    out = ensure_dir(arm_dir(region, arm, seed))
    production = region_store(COHORT, region).resolve()
    if production in out.resolve().parents:
        raise SystemExit(f"refusing to write inside a production store: {out}")
    base = _filters()[arm]
    record: dict = {"region": region, "arm": arm, "seed": seed, "analysis": TARGETS[region],
                    "is_production_setting": arm == PRODUCTION_ARM[TARGETS[region]]}

    def recorded_filter(tc, tt):
        record.update(n_transcripts_before=int(len(tt)), n_genes_before=int(tt["gene_id"].nunique()))
        if base is not None:
            tc, tt = base(tc, tt)
        per_gene = tt["gene_id"].value_counts()
        record.update(n_transcripts=int(len(tt)), n_genes=int(len(per_gene)),
                      n_multi_isoform_genes=int((per_gene >= 2).sum()))
        tt[[_tx_id_col(tt), "gene_id"]].to_parquet(out / "kept_transcripts.parquet", index=False)
        return tc, tt

    print(f"[filter-arm {grid} {index}/{len(fits)}] {COHORT}/{region} arm={arm} seed={seed} -> {out}",
          flush=True)
    t0 = time.time()
    if region == "caudate_sczd":
        run_brainseq_caudate_sczd(transcript_filter=recorded_filter, out=out, random_state=seed)
    else:
        run_brainseq_region(region, transcript_filter=recorded_filter, out=out, random_state=seed)
    record["wall_seconds"] = round(time.time() - t0, 1)
    if arm == "switching":
        record["switching_params"] = SWITCHING
    (out / "fit.json").write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2), flush=True)


# --------------------------------------------------------------------------- #
# Compare
# --------------------------------------------------------------------------- #
def load_fit(region: str, arm: str, seed: int = PRODUCTION_SEED) -> dict | None:
    """`arm == "production"` reads the committed store; otherwise the arm directory."""
    d = (region_store(COHORT, region, "isograph_vae") if arm == "production"
         else arm_dir(region, arm, seed))
    if not (d / "modules.parquet").exists():
        return None
    trait_path = d / TRAIT_FILE[TARGETS[region]]
    fit = {
        "region": region, "arm": arm, "seed": seed, "dir": d,
        "modules": pd.read_parquet(d / "modules.parquet", columns=["gene_id", "module_id"]).astype(str),
        "trait": (pd.read_parquet(trait_path) if trait_path.exists()
                  else pd.DataFrame(columns=["module_id", "effect", "fdr"])),
        "kept": (pd.read_parquet(d / "kept_transcripts.parquet")
                 if (d / "kept_transcripts.parquet").exists() else None),
        "meta": json.loads((d / "fit.json").read_text()) if (d / "fit.json").exists() else {},
    }
    cal = d / "calibration.parquet"
    fit["alpha"] = (pd.read_parquet(cal).get("selected_alpha_abundance", pd.Series([np.nan])).iloc[0]
                    if cal.exists() else np.nan)
    fit["trait"]["module_id"] = fit["trait"]["module_id"].astype(str)
    return fit


def _sig_modules(trait: pd.DataFrame, fdr: float = FDR) -> pd.DataFrame:
    return trait.loc[trait["fdr"] < fdr].set_index("module_id")


def fit_summary(f: dict, fdr: float = FDR) -> dict:
    sizes = f["modules"]["module_id"].value_counts()
    sig = _sig_modules(f["trait"], fdr)
    return {
        "region": f["region"], "arm": f["arm"], "seed": f.get("seed", PRODUCTION_SEED),
        "n_transcripts": f["meta"].get("n_transcripts"), "n_genes_input": f["meta"].get("n_genes"),
        "n_multi_isoform_genes": f["meta"].get("n_multi_isoform_genes"),
        "n_modules": int(len(sizes)), "n_genes_assigned": int(sizes.sum()),
        "median_module_size": float(sizes.median()) if len(sizes) else np.nan,
        "max_module_size": int(sizes.max()) if len(sizes) else 0,
        "largest_module_frac": float(sizes.max() / sizes.sum()) if len(sizes) else np.nan,
        "n_modules_ge_900": int((sizes >= GIANT_GENES).sum()),
        "n_trait_tested": int(len(f["trait"])), "n_trait_sig": int(len(sig)),
        "n_genes_in_trait_sig": int(f["modules"]["module_id"].isin(sig.index).sum()),
        "selected_alpha_abundance": f["alpha"],
    }


def _retention(src_sets, src_sig, dst_sets, dst_trait, fdr):
    """Per significant source module: best-Jaccard destination counterpart, and whether that
    counterpart is significant with the same sign."""
    dst = dst_trait.set_index("module_id")
    best_j, kept = [], 0
    for mid, row in src_sig.iterrows():
        genes = src_sets.get(mid, set())
        if not genes or dst_sets.empty:
            best_j.append(0.0)
            continue
        jac = dst_sets.apply(lambda g: len(genes & g) / len(genes | g))
        best = jac.idxmax()
        best_j.append(float(jac.max()))
        if best in dst.index and dst.loc[best, "fdr"] < fdr and \
                np.sign(dst.loc[best, "effect"]) == np.sign(row["effect"]):
            kept += 1
    return kept, (float(np.median(best_j)) if best_j else np.nan)


def _gene_concordance(ref: dict, qry: dict, m: pd.DataFrame, fdr: float) -> dict:
    """Module-free trait agreement over genes assigned in both fits: each gene inherits its
    module's trait effect and significance. Best-Jaccard retention is base-rate blind (when most
    `qry` modules are significant any counterpart qualifies); here a `ref`-significant gene
    counts as kept only if it is significant with the same sign in `qry`, read against the
    `qry` significance rate, with a Fisher odds ratio for the sig x sig overlap."""
    from scipy.stats import fisher_exact, spearmanr

    def per_gene(f, side):
        t = f["trait"][["module_id", "effect", "fdr"]].rename(
            columns={"module_id": f"module_id_{side}", "effect": f"effect_{side}", "fdr": f"fdr_{side}"})
        return t
    g = m.merge(per_gene(ref, "ref"), on="module_id_ref", how="left") \
         .merge(per_gene(qry, "qry"), on="module_id_qry", how="left")
    sr, sq = (g["fdr_ref"] < fdr).to_numpy(), (g["fdr_qry"] < fdr).to_numpy()
    same = np.sign(g["effect_ref"]).to_numpy() == np.sign(g["effect_qry"]).to_numpy()
    ok = g["effect_ref"].notna() & g["effect_qry"].notna()
    rho = (float(spearmanr(g.loc[ok, "effect_ref"], g.loc[ok, "effect_qry"]).statistic)
           if ok.sum() > 2 else np.nan)
    tab = [[int((sr & sq).sum()), int((sr & ~sq).sum())], [int((~sr & sq).sum()), int((~sr & ~sq).sum())]]
    odds = float(fisher_exact(tab)[0]) if sr.any() and sq.any() else np.nan
    return {
        "gene_effect_spearman": rho,
        "n_ref_sig_genes": int(sr.sum()),
        "frac_ref_sig_genes_kept_same_sign": float((sr & sq & same).sum() / sr.sum()) if sr.any() else np.nan,
        "qry_sig_gene_rate": float(sq.mean()) if len(sq) else np.nan,
        "sig_gene_odds_ratio": odds,
    }


def compare_fits(ref: dict, qry: dict, fdr: float = FDR) -> dict:
    """Partition agreement and trait-module correspondence between two fits of one target."""
    r, q = ref["modules"], qry["modules"]
    m = r.merge(q, on="gene_id", suffixes=("_ref", "_qry"))
    r_genes, q_genes = set(r["gene_id"]), set(q["gene_id"])
    r_sets = r.groupby("module_id")["gene_id"].apply(set)
    q_sets = q.groupby("module_id")["gene_id"].apply(set)
    r_sig, q_sig = _sig_modules(ref["trait"], fdr), _sig_modules(qry["trait"], fdr)
    r_kept, r_bj = _retention(r_sets, r_sig, q_sets, qry["trait"], fdr)
    q_kept, q_bj = _retention(q_sets, q_sig, r_sets, ref["trait"], fdr)
    r_tg = set(r.loc[r["module_id"].isin(r_sig.index), "gene_id"])
    q_tg = set(q.loc[q["module_id"].isin(q_sig.index), "gene_id"])
    out = {
        "region": ref["region"], "ref": ref["arm"], "qry": qry["arm"],
        "ref_seed": ref.get("seed", PRODUCTION_SEED), "qry_seed": qry.get("seed", PRODUCTION_SEED),
        "n_genes_both": int(len(m)),
        "assigned_gene_jaccard": len(r_genes & q_genes) / max(len(r_genes | q_genes), 1),
        "ari": float(adjusted_rand_score(m["module_id_ref"], m["module_id_qry"])) if len(m) else np.nan,
        "nmi": float(normalized_mutual_info_score(m["module_id_ref"], m["module_id_qry"])) if len(m) else np.nan,
        "n_trait_sig_ref": int(len(r_sig)), "n_trait_sig_qry": int(len(q_sig)),
        "n_ref_sig_retained_in_qry": int(r_kept), "median_best_jaccard_ref_sig": r_bj,
        "n_qry_sig_retained_in_ref": int(q_kept), "median_best_jaccard_qry_sig": q_bj,
        "trait_gene_jaccard": (len(r_tg & q_tg) / len(r_tg | q_tg)) if (r_tg | q_tg) else np.nan,
        **_gene_concordance(ref, qry, m, fdr),
        "transcript_jaccard": np.nan,
    }
    if ref["kept"] is not None and qry["kept"] is not None:
        a = set(ref["kept"].iloc[:, 0].astype(str))
        b = set(qry["kept"].iloc[:, 0].astype(str))
        out["transcript_jaccard"] = len(a & b) / max(len(a | b), 1)
    return out


def _md(df: pd.DataFrame, floats: int = 3) -> str:
    def cell(v):
        if isinstance(v, (float, np.floating)):
            return "" if pd.isna(v) else f"{v:.{floats}g}"
        return "" if v is None else str(v)
    cols = list(df.columns)
    return "\n".join(["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols),
                      *["| " + " | ".join(cell(v) for v in row) + " |"
                        for row in df.itertuples(index=False, name=None)]])


def _seed_floor(fdr: float, missing: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Every arm at every seed: per-fit summaries, all pairwise comparisons labelled
    `within_arm` (seed-to-seed variation, the floor) or `between_arm` (filter + seed), and a
    per (region, ref arm, qry arm) aggregate. Pairs keep the production-setting arm as `ref`."""
    rows, pairs = [], []
    for region, analysis in TARGETS.items():
        arms = [PRODUCTION_ARM[analysis], "switching"] + (["abundance"] if analysis == "sczd" else [])
        fits = {}
        for arm in arms:
            for seed in (PRODUCTION_SEED, *EXTRA_SEEDS):
                f = load_fit(region, arm, seed)
                if f is None:
                    missing.append(f"{region}/{arm}/seed{seed}")
                    continue
                fits[(arm, seed)] = f
                rows.append(fit_summary(f, fdr))
        for (ka, fa), (kb, fb) in itertools.combinations(fits.items(), 2):
            pairs.append({"kind": "within_arm" if ka[0] == kb[0] else "between_arm",
                          **compare_fits(fa, fb, fdr)})
    fits_df, pairs_df = pd.DataFrame(rows), pd.DataFrame(pairs)
    if pairs_df.empty:
        return fits_df, pairs_df, pd.DataFrame()

    def by_seed(col):
        return lambda s: " / ".join(str(int(v)) for v in s)
    per_fit = (fits_df.sort_values("seed").groupby(["region", "arm"], sort=False)
               .agg(seeds=("seed", by_seed("seed")), n_modules=("n_modules", by_seed("n_modules")),
                    n_trait_sig=("n_trait_sig", by_seed("n_trait_sig")),
                    n_genes_in_trait_sig=("n_genes_in_trait_sig", by_seed("n_genes_in_trait_sig")),
                    n_modules_ge_900=("n_modules_ge_900", "max"))
               .reset_index())
    agg = (pairs_df.groupby(["region", "kind", "ref", "qry"], sort=False)
           .agg(n_pairs=("ari", "size"), ari_median=("ari", "median"), ari_min=("ari", "min"),
                ari_max=("ari", "max"), gene_effect_spearman_median=("gene_effect_spearman", "median"),
                frac_ref_sig_genes_kept_same_sign_median=("frac_ref_sig_genes_kept_same_sign", "median"),
                qry_sig_gene_rate_median=("qry_sig_gene_rate", "median"),
                trait_gene_jaccard_median=("trait_gene_jaccard", "median"))
           .reset_index())
    return per_fit, pairs_df, agg


def compare(fdr: float = FDR) -> None:
    root = ensure_dir(stage_out("modules.filter_arms"))
    summaries, pairs, missing = [], [], []
    for region, analysis in TARGETS.items():
        arms = ["production", PRODUCTION_ARM[analysis], "switching"]
        if analysis == "sczd":
            arms.append("abundance")
        fits = {}
        for arm in arms:
            f = load_fit(region, arm)
            if f is None:
                missing.append(f"{region}/{arm}")
            else:
                fits[arm] = f
                summaries.append(fit_summary(f, fdr))
        for a, b in itertools.combinations(arms, 2):
            if a in fits and b in fits:
                pairs.append(compare_fits(fits[a], fits[b], fdr))
    fits_df, pairs_df = pd.DataFrame(summaries), pd.DataFrame(pairs)
    fits_df.to_parquet(root / "fits.parquet", index=False)
    pairs_df.to_parquet(root / "comparisons.parquet", index=False)

    lines = [
        "# Transcript-filter arms (BrainSEQ)",
        "",
        "`transcript_filter_arms.py`: the production IsoGraph fit re-run with only the transcript "
        "filter changed, written outside the production stores. `production` is the committed "
        "fit. Aging production uses the abundance filter (count > 10 in ≥ 70%), so its "
        "`abundance` refit is the noise floor; SCZD production is unfiltered, so its `none` refit "
        "is the floor and `abundance` is the BrainSEQ method applied to SCZD.",
        "",
        "Switching filter: gene count ≥ 10 in ≥ 70% of samples; transcript count ≥ 10 and share "
        "of its gene ≥ 0.10, each in ≥ 10% of samples.",
        "",
        f"Trait tables: aging `age_linear` (linear age), SCZD `diagnosis_assoc`; significant = FDR < {fdr}.",
        "",
        "## Fits",
        "",
        _md(fits_df),
        "",
        "## Pairwise",
        "",
        "Read every production-vs-arm row against the production-vs-floor row for the same "
        "region. A floor at ARI = 1 means the refit is deterministic at the production seed "
        "(identical edges), NOT that seed or input-perturbation variation is zero; the floor "
        "bounds refit reproducibility only. `n_ref_sig_retained_in_qry`: significant `ref` "
        "modules whose best-Jaccard `qry` counterpart is significant with the same sign (and the "
        "reverse column); it is base-rate blind, so read it with `median_best_jaccard_*` and the "
        "gene-level columns. `frac_ref_sig_genes_kept_same_sign` against `qry_sig_gene_rate` "
        "(its chance level) and `sig_gene_odds_ratio` are module-free; `gene_effect_spearman` "
        "correlates each gene's module trait effect across the two fits. "
        "`trait_gene_jaccard`: overlap of the genes in significant modules. ARI/NMI and the "
        "gene-level columns are over genes assigned in both fits.",
        "",
        _md(pairs_df),
        "",
    ]
    seed_fits, seed_pairs, seed_agg = _seed_floor(fdr, missing)
    if not seed_agg.empty:
        seed_pairs.to_parquet(root / "seed_pairs.parquet", index=False)
        seed_agg.to_parquet(root / "seed_floor.parquet", index=False)
        lines += [
            "## Seed floor",
            "",
            f"Every arm refit at seeds {PRODUCTION_SEED}, {', '.join(map(str, EXTRA_SEEDS))} "
            "(`random_state` drives VAE initialisation, minibatch order and Leiden). `within_arm` "
            "pairs are seed-to-seed variation for one filter: the floor. `between_arm` pairs "
            "change the filter (and, off the diagonal, the seed). A filter effect is a "
            "`between_arm` value outside the `within_arm` range of both arms.",
            "",
            "Per arm, values in seed order:",
            "",
            _md(seed_fits),
            "",
            _md(seed_agg),
            "",
        ]
    if missing:
        lines += ["Missing fits: " + ", ".join(missing), ""]
    (root / "TRANSCRIPT_FILTER_ARMS.md").write_text("\n".join(lines))
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(fits_df.to_string(index=False))
        print(pairs_df.to_string(index=False))
    if missing:
        print("missing:", ", ".join(missing))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", nargs="?", choices=["fit", "compare"])
    ap.add_argument("--index", type=int)
    ap.add_argument("--grid", choices=sorted(GRIDS), default="filter",
                    help="filter: the arms at the production seed; seed: the same arms at extra seeds")
    ap.add_argument("--fdr", type=float, default=FDR)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for i, (region, arm, seed) in enumerate(GRIDS[args.grid], 1):
            print(i, region, arm, seed)
    elif args.stage == "fit":
        if args.index is None:
            ap.error("fit needs --index")
        fit(args.index, args.grid)
    elif args.stage == "compare":
        compare(args.fdr)
    else:
        ap.error("pass fit --index N, compare, or --list")


if __name__ == "__main__":
    main()
