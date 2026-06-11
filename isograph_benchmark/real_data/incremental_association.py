"""Switch-vs-abundance incremental association — strengths/limitations of IsoGraph.

Two complementary, reproducible analyses that ask whether isoform COMPOSITION
(IsoGraph's switch channel) carries phenotype signal beyond total ABUNDANCE,
using IsoGraph's own per-sample switch and abundance features (identical
preprocessing, so the contrast is apples-to-apples):

1. gene_level_deconfounded  — the rigorous, de-confounded test. Per gene:
     switch_g    ~ phenotype + abundance_g + covariates   (composition | abundance)
     abundance_g ~ phenotype + switch_g    + covariates   (abundance | composition)
   No clustering / eigengenes, so it removes the module-definition circularity.
   "composition-unique" genes (switch significant, abundance not) are genuine
   DTU-without-DGE: isoform-switching regulation a DGE/WGCNA pipeline cannot see.

2. module_level_incremental — the eigengene-level version (switch eigengene vs
   abundance eigengene of the SAME module genes). Reported for completeness/
   contrast: it INFLATES composition-uniqueness because modules are defined by
   switch co-expression (coherent switch eigengene, incoherent abundance
   eigengene). Keep it next to the gene-level result so the circularity is visible.

phenotype = SCZD diagnosis (brainseq-sczd) or non-linear age (brainseq-aging,
df=3 spline F-test). Outputs land in <artifact_dir>/incremental_association/.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from patsy import dmatrix
from scipy import stats

from isograph.io.artifacts import load_dataset_bundle

from isograph_benchmark.paths import ensure_dir
from isograph_benchmark.real_data.run_models import (
    GTEX_REGIONS,
    diagnosis_association,
    spline_age_association,
)
from isograph_benchmark.real_data.sweep_leiden import (
    AGING_COVARIATE_COLS,
    GTEX_AGE_COL,
    GTEX_COVARIATE_COLS,
    SCZD_COVARIATE_COLS,
    _artifact_dir,
    _bundle_path,
    _pivot_eigengenes,
)

_META_COLS = ("feature_id", "gene_id", "feature_type", "n_transcripts")


def _covariate_cols(analysis: str) -> list[str]:
    if analysis == "brainseq-sczd":
        return SCZD_COVARIATE_COLS
    if analysis == "gtex-aging":
        return GTEX_COVARIATE_COLS
    return AGING_COVARIATE_COLS


def _age_col(analysis: str) -> str:
    return GTEX_AGE_COL if analysis == "gtex-aging" else "Age"


def _load(analysis: str, region: str | None, variant: str):
    artifact_dir = _artifact_dir(analysis, region, variant)
    fs = pd.read_parquet(artifact_dir / "feature_scores.parquet")
    modules = pd.read_parquet(artifact_dir / "modules.parquet")
    bundle = load_dataset_bundle(_bundle_path(analysis, region))
    return artifact_dir, fs, modules, bundle


def _sample_cols(fs: pd.DataFrame, sample_ids: set[str]) -> list[str]:
    return [c for c in fs.columns if c not in _META_COLS and c in sample_ids]


def _channel_matrix(fs: pd.DataFrame, feature_type: str, samp: list[str]) -> pd.DataFrame:
    """gene_id-indexed (genes × samples) matrix for one feature channel."""
    sub = fs[fs["feature_type"] == feature_type].set_index("gene_id")
    return sub[samp]


def _eigengene_table(channel: pd.DataFrame, modules: pd.DataFrame, samp: list[str]) -> pd.DataFrame:
    """module_id × samples eigengene = mean of the channel over module genes
    (mirrors compute_trait_associations' switch eigengene construction)."""
    rows = {}
    for mid in sorted(modules["module_id"].unique()):
        genes = modules.loc[modules["module_id"] == mid, "gene_id"]
        present = channel.index.intersection(genes)
        if len(present) == 0:
            continue
        rows[mid] = channel.loc[present].to_numpy(float).mean(axis=0)
    eg = pd.DataFrame(rows, index=samp).T
    eg.index.name = "module_id"
    return eg.reset_index()


# --------------------------------------------------------------------------- #
# 1. Gene-level de-confounded test
# --------------------------------------------------------------------------- #
def _ols_coef_p(X: np.ndarray, y: np.ndarray, idx: int) -> float:
    if not (np.all(np.isfinite(X)) and np.all(np.isfinite(y))):
        return np.nan
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    dof = len(y) - rank
    if dof <= 0:
        return np.nan
    s2 = float((y - X @ beta) @ (y - X @ beta)) / dof
    se = np.sqrt(max(s2 * np.linalg.pinv(X.T @ X)[idx, idx], 1e-30))
    return float(2 * stats.t.sf(abs(beta[idx] / se), dof))


def _ftest_block_p(Xfull: np.ndarray, Xred: np.ndarray, y: np.ndarray) -> float:
    if not (np.all(np.isfinite(Xfull)) and np.all(np.isfinite(y))):
        return np.nan
    bf, _, rf, _ = np.linalg.lstsq(Xfull, y, rcond=None)
    br, _, rr, _ = np.linalg.lstsq(Xred, y, rcond=None)
    rss_f = float(((y - Xfull @ bf) ** 2).sum())
    rss_r = float(((y - Xred @ br) ** 2).sum())
    dfn, dfd = max(rf - rr, 1), max(len(y) - rf, 1)
    F = ((rss_r - rss_f) / dfn) / (rss_f / dfd)
    return float(stats.f.sf(F, dfn, dfd)) if np.isfinite(F) else np.nan


def gene_level_deconfounded(
    analysis: str,
    fs: pd.DataFrame,
    bundle,
    fdr_alpha: float = 0.10,
    min_samples: int = 30,
) -> pd.DataFrame:
    st = bundle.sample_table
    sample_ids = set(st["sample_id"].astype(str))
    samp = _sample_cols(fs, sample_ids)
    st = st.set_index(st["sample_id"].astype(str)).loc[samp]
    covs = _covariate_cols(analysis)

    SW = _channel_matrix(fs, "switch", samp)
    AB = _channel_matrix(fs, "abundance", samp)
    genes = SW.index.intersection(AB.index)
    SWm = SW.loc[genes].to_numpy(float).T   # samples × genes
    ABm = AB.loc[genes].to_numpy(float).T

    covdf = pd.get_dummies(st[[c for c in covs if c in st.columns]].copy(), drop_first=True)
    C = np.hstack([np.ones((len(st), 1)), covdf.to_numpy(dtype=float)])
    C_ok = np.all(np.isfinite(C), axis=1)

    pA = np.full(len(genes), np.nan)
    pB = np.full(len(genes), np.nan)

    if analysis == "brainseq-sczd":
        pheno = (st["Dx"].astype(str).to_numpy() == "SCZD").astype(float)[:, None]
        dx_idx = C.shape[1]
        for j in range(len(genes)):
            ys, ya = SWm[:, j], ABm[:, j]
            m = np.isfinite(ys) & np.isfinite(ya) & C_ok
            if m.sum() < min_samples:
                continue
            pA[j] = _ols_coef_p(np.hstack([C[m], pheno[m], ya[m, None]]), ys[m], dx_idx)
            pB[j] = _ols_coef_p(np.hstack([C[m], pheno[m], ys[m, None]]), ya[m], dx_idx)
    else:
        age = st[_age_col(analysis)].to_numpy(float)
        age_z = (age - age.mean()) / age.std()
        knots = np.quantile(age_z, [0.5])
        bd = [age_z.min(), age_z.max()]
        B = np.asarray(dmatrix(
            f"cr(age_z, knots={list(knots)}, lower_bound={bd[0]}, upper_bound={bd[1]}) - 1",
            {"age_z": age_z}, return_type="matrix"))
        for j in range(len(genes)):
            ys, ya = SWm[:, j], ABm[:, j]
            m = np.isfinite(ys) & np.isfinite(ya) & C_ok
            if m.sum() < min_samples:
                continue
            pA[j] = _ftest_block_p(np.hstack([C[m], B[m], ya[m, None]]), np.hstack([C[m], ya[m, None]]), ys[m])
            pB[j] = _ftest_block_p(np.hstack([C[m], B[m], ys[m, None]]), np.hstack([C[m], ys[m, None]]), ya[m])

    ok = np.isfinite(pA) & np.isfinite(pB)
    res = pd.DataFrame({
        "gene_id": genes[ok],
        "p_switch_given_abund": pA[ok],
        "p_abund_given_switch": pB[ok],
    })
    res["fdr_switch_given_abund"] = stats.false_discovery_control(res["p_switch_given_abund"], method="bh")
    res["fdr_abund_given_switch"] = stats.false_discovery_control(res["p_abund_given_switch"], method="bh")
    sw = res["fdr_switch_given_abund"] <= fdr_alpha
    ab = res["fdr_abund_given_switch"] <= fdr_alpha
    res["category"] = np.select(
        [sw & ~ab, ab & ~sw, sw & ab], ["composition_unique", "abundance_unique", "both"], "neither"
    )
    return res.sort_values("p_switch_given_abund").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# 2. Module-level incremental (eigengene switch vs abundance) — for contrast
# --------------------------------------------------------------------------- #
def module_level_incremental(
    analysis: str,
    fs: pd.DataFrame,
    modules: pd.DataFrame,
    bundle,
    fdr_alpha: float = 0.10,
) -> pd.DataFrame:
    st = bundle.sample_table
    sample_ids = set(st["sample_id"].astype(str))
    samp = _sample_cols(fs, sample_ids)
    covs = _covariate_cols(analysis)

    eg_switch = _pivot_eigengenes(_eigengene_table(_channel_matrix(fs, "switch", samp), modules, samp))
    eg_abund = _pivot_eigengenes(_eigengene_table(_channel_matrix(fs, "abundance", samp), modules, samp))

    def assoc(egp):
        if analysis == "brainseq-sczd":
            d = diagnosis_association(egp, st, covariate_cols=covs)
            return d.set_index("module_id")["fdr"]
        d = spline_age_association(egp, st, covariate_cols=covs, age_col=_age_col(analysis))
        return d.drop_duplicates("module_id").set_index("module_id")["fdr_ftest"]

    res = pd.DataFrame({"switch_fdr": assoc(eg_switch), "abund_fdr": assoc(eg_abund)}).dropna()
    res["n_genes"] = [int((modules["module_id"] == m).sum()) for m in res.index]
    sw = res["switch_fdr"] <= fdr_alpha
    ab = res["abund_fdr"] <= fdr_alpha
    res["category"] = np.select(
        [sw & ~ab, ab & ~sw, sw & ab], ["composition_unique", "abundance_unique", "both"], "neither"
    )
    return res.reset_index().rename(columns={"index": "module_id"})


def _summary(gene_level: pd.DataFrame, module_level: pd.DataFrame) -> dict:
    gc = gene_level["category"].value_counts().to_dict()
    mc = module_level["category"].value_counts().to_dict()
    return {
        "gene_level": {"n_tested": int(len(gene_level)), **{k: int(gc.get(k, 0)) for k in
                       ["composition_unique", "abundance_unique", "both", "neither"]}},
        "module_level": {"n_tested": int(len(module_level)), **{k: int(mc.get(k, 0)) for k in
                         ["composition_unique", "abundance_unique", "both", "neither"]}},
    }


def run_analysis(analysis: str, region: str | None, variant: str, fdr_alpha: float) -> dict:
    label = f"{analysis}/{region}" if region else analysis
    artifact_dir, fs, modules, bundle = _load(analysis, region, variant)
    print(f"[{label}] gene-level de-confounded test ...", flush=True)
    gene_level = gene_level_deconfounded(analysis, fs, bundle, fdr_alpha=fdr_alpha)
    print(f"[{label}] module-level incremental ...", flush=True)
    module_level = module_level_incremental(analysis, fs, modules, bundle, fdr_alpha=fdr_alpha)

    out = ensure_dir(artifact_dir / "incremental_association")
    gene_level.to_parquet(out / "gene_level.parquet", index=False, compression="zstd")
    module_level.to_parquet(out / "module_level.parquet", index=False, compression="zstd")
    summary = _summary(gene_level, module_level)
    summary.update({"analysis": analysis, "region": region, "variant": variant, "fdr_alpha": fdr_alpha})
    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    g = summary["gene_level"]
    m = summary["module_level"]
    print(f"[{label}] GENE-LEVEL  (n={g['n_tested']}): composition_unique={g['composition_unique']} "
          f"abundance_unique={g['abundance_unique']} both={g['both']}")
    print(f"[{label}] MODULE-LEVEL(n={m['n_tested']}): composition_unique={m['composition_unique']} "
          f"abundance_unique={m['abundance_unique']} both={m['both']}")
    print(f"[{label}] written to {out}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis", choices=["brainseq-sczd", "brainseq-aging", "gtex-aging"])
    parser.add_argument("--region", action="append", dest="regions",
                        help="For brainseq-aging/gtex-aging: region(s). Repeatable. "
                             "Default: all brainseq (3) or all gtex (13) regions.")
    parser.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    parser.add_argument("--fdr", type=float, default=0.10)
    args = parser.parse_args()

    if args.analysis == "brainseq-sczd":
        run_analysis("brainseq-sczd", None, args.variant, args.fdr)
    elif args.analysis == "gtex-aging":
        for region in (args.regions or GTEX_REGIONS):
            run_analysis("gtex-aging", region, args.variant, args.fdr)
    else:
        for region in (args.regions or ["caudate", "hippocampus", "dlpfc"]):
            run_analysis("brainseq-aging", region, args.variant, args.fdr)


if __name__ == "__main__":
    main()
