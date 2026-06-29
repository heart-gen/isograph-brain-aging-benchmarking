from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
from patsy import dmatrix
from scipy import stats

from isograph.io.artifacts import load_dataset_bundle
from isograph.models.vae import VaeNetworkModel
from isograph.workflow.config import VaeModelConfig
from isograph_benchmark.paths import ensure_dir, rel

# Age evaluation probabilities for spline projection (p10, p25, p50, p75, p90)
AGE_PROBS = np.array([0.25, 0.50, 0.75])

# Best Leiden resolution per region, selected data-driven by the Part 1 sweep
# (isograph_benchmark.real_data.sweep_leiden --write-best; GO-enrichment count
# under the giant_fraction <= 0.30 constraint). Used as the default resolution
# for the with-abundance refit so it matches the chosen standard partition.
# Override per run with --leiden-resolution. Keys: brainseq region dir names.
BEST_LEIDEN_RESOLUTION = {
    "caudate": 2.25,       # brainseq aging caudate
    "hippocampus": 2.0,    # brainseq aging hippocampus
    "dlpfc": 3.0,          # brainseq aging dlpfc
    "caudate_sczd": 2.0,   # brainseq SCZD+Control caudate
}

# Promoted single-config VAE defaults for every real-data fit (validation gate B,
# 2026-06-22). Spread into each production VaeModelConfig so the shipped pipeline
# matches the validated config and a stranger needs no per-dataset tuning:
#   - grad_clip_norm=1.0 (B.2 PASS): one fixed lr=1e-3 + global grad-norm clipping
#     trained all 6 trust-funnel regions AND the GTEx nucleus_accumbens that diverged
#     at lr=1e-3 pre-clip (0/7 diverged, no OOM). Removes the hand-tuned GTEx lr=3e-4.
#   - estimability switch reliability (S5): the sole *positive* split-half-stability
#     lever in the A/B (covariate-free minor-isoform-usage downweighting). degradation/
#     differential-TIN were negative and are not used here.
# The synthetic benchmark grid is intentionally NOT changed (its config defaults stay
# off until synthetic + Type-I are re-confirmed unaffected).
_PROMOTED_VAE = dict(
    grad_clip_norm=1.0,
    switch_reliability_weighting=True,
    switch_reliability_source="estimability",
    switch_estimability_min_minor_usage=0.1,
)


# Canonical production Leiden resolution for the standard variant. See the wiki
# (Tuning-and-Stability-Selection) for the rationale and the biology-driven sweep.
CANONICAL_LEIDEN_RESOLUTION = 5.0


def _isograph_out_subdir(leiden_resolution: float | None) -> str:
    """Output subdir name for a standard isograph_vae fit.

    With no override (None) or the canonical resolution (5.0) this is the
    canonical ``isograph_vae`` dir that the GWAS and trust-funnel cascades
    consume. Any *other* explicit ``leiden_resolution`` is written to a
    resolution-suffixed sibling (e.g. ``isograph_vae_res5`` for 5.0 is the
    canonical dir, ``isograph_vae_res2`` for 2.0) so a non-canonical resolution
    is a side-by-side comparison set and never clobbers the canonical modules.
    The suffix encodes the resolution with '.' -> 'p' (e.g. 2.25 ->
    isograph_vae_res2p25).
    """
    if leiden_resolution is None or leiden_resolution == CANONICAL_LEIDEN_RESOLUTION:
        return "isograph_vae"
    token = f"{leiden_resolution:g}".replace(".", "p")
    return f"isograph_vae_res{token}"


def _filter_expressed_transcripts(
    transcript_counts: np.ndarray,
    transcript_table: pd.DataFrame,
    min_count: float = 10.0,
    min_fraction: float = 0.70,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Keep transcripts with count > min_count in >= min_fraction of samples.

    Mirrors the filterByExpr-style filter used in the SCZD bundle creation,
    adapted for a continuous covariate (aging spline) by treating all samples
    as a single group.  Drops lowly-expressed isoforms that add noise to the
    switch-coordinate computation and would otherwise inflate the (n_features²)
    similarity matrix.
    """
    n_samples = transcript_counts.shape[1]
    n_pass = int(min_fraction * n_samples)
    tx_pass = (transcript_counts > min_count).sum(axis=1) >= n_pass
    n_orig_tx = len(transcript_table)
    n_kept_tx = int(tx_pass.sum())
    n_orig_genes = transcript_table["gene_id"].nunique()
    n_kept_genes = transcript_table.loc[tx_pass, "gene_id"].nunique()
    print(
        f"  transcript filter (count>{min_count:.0f}, ≥{min_fraction:.0%} of {n_samples} samples): "
        f"{n_kept_tx}/{n_orig_tx} transcripts, {n_kept_genes}/{n_orig_genes} genes retained",
        flush=True,
    )
    return transcript_counts[tx_pass], transcript_table.loc[tx_pass].reset_index(drop=True)


AGE_LABELS = ["p25", "p50", "p75"]

GTEX_REGIONS = [
    "amygdala",
    "anterior_cingulate_cortex_ba24",
    "caudate_basal_ganglia",
    "cerebellar_hemisphere",
    "cerebellum",
    "cortex",
    "frontal_cortex_ba9",
    "hippocampus",
    "hypothalamus",
    "nucleus_accumbens_basal_ganglia",
    "putamen_basal_ganglia",
    "spinal_cord_cervical_c_1",
    "substantia_nigra",
]


def _standardize(x: np.ndarray) -> np.ndarray:
    return (x - x.mean()) / x.std()


def _natural_spline_basis(age_z: np.ndarray, knots: np.ndarray, boundary: np.ndarray) -> np.ndarray:
    """Natural cubic spline basis via patsy cr() — linear beyond boundary knots."""
    basis = dmatrix(
        f"cr(age_z, knots={list(knots)}, lower_bound={boundary[0]}, upper_bound={boundary[1]}) - 1",
        {"age_z": age_z},
        return_type="matrix",
    )
    return np.asarray(basis)


def linear_age_association(
    eigengenes: pd.DataFrame, sample_table: pd.DataFrame, age_col: str = "Age"
) -> pd.DataFrame:
    """Pearson correlation of each module eigengene with the age column."""
    merged = eigengenes.merge(sample_table[["sample_id", age_col]], on="sample_id", how="inner")
    age = merged[age_col].to_numpy(dtype=float)
    module_cols = [c for c in eigengenes.columns if c != "sample_id"]
    rows = []
    for col in module_cols:
        eg = merged[col].to_numpy(dtype=float)
        mask = np.isfinite(eg) & np.isfinite(age)
        if mask.sum() < 10:
            continue
        r, p = stats.pearsonr(age[mask], eg[mask])
        rows.append({"module_id": col, "trait": "Age_linear", "effect": r, "pvalue": p, "n": mask.sum()})
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = stats.false_discovery_control(result["pvalue"], method="bh")
    return result


def spline_age_association(
    eigengenes: pd.DataFrame, sample_table: pd.DataFrame, covariate_cols: list[str],
    age_col: str = "Age",
) -> pd.DataFrame:
    """
    Non-linear age association via natural cubic spline, projected to 3 age points.

    Fits: eigengene ~ cr(age_z, knots=[median]) + covariates  (df=3, K=3 coeffs)
    Projects spline coefficients to age_eval = qnorm([0.25, 0.50, 0.75]) = early/mid/late.

    df=3 (one interior knot) was selected empirically: it matches df=4 and beats
    df=5 on module-level F-test power across all three aging regions while being the
    most parsimonious / interpretable trajectory (one bend). See AGE_PROBS/AGE_LABELS.

    Returns one row per (module, age_label) with per-point effect/se/z, the F-test
    columns (pvalue_ftest, fdr_ftest) testing the joint spline component against the
    covariate-only reduced model (USE fdr_ftest for module-level significance — the
    per-point z-test is over-conservative because the points are one spline), and
    the raw spline coefficients (coef_1..K) + full coefficient covariance
    (cov_age_ij), repeated per module row, so the trajectory can be re-projected to
    any age grid (or fed to mash) downstream.
    """
    keep_cols = ["sample_id", age_col] + [c for c in covariate_cols if c in sample_table.columns]
    merged = eigengenes.merge(sample_table[keep_cols], on="sample_id", how="inner")
    merged = merged.dropna(subset=[age_col])

    age_z = _standardize(merged[age_col].to_numpy(dtype=float))
    ns_knots = np.quantile(age_z, [0.5])
    ns_boundary = np.array([age_z.min(), age_z.max()])

    # Build spline basis for observed samples
    B_obs = _natural_spline_basis(age_z, ns_knots, ns_boundary)
    n_spline = B_obs.shape[1]

    # Build projection basis at 5 evaluation points
    age_eval = stats.norm.ppf(AGE_PROBS)
    age_eval_clipped = np.clip(age_eval, ns_boundary[0], ns_boundary[1])
    B_proj = _natural_spline_basis(age_eval_clipped, ns_knots, ns_boundary)

    # Build covariate design matrix
    available_covariates = [c for c in covariate_cols if c in merged.columns and c != age_col]
    covariate_df = merged[available_covariates].copy()
    covariate_df = pd.get_dummies(covariate_df, drop_first=True)
    covariate_arr = covariate_df.to_numpy(dtype=float)

    intercept = np.ones((len(merged), 1))
    X = np.hstack([intercept, B_obs, covariate_arr])
    X_reduced = np.hstack([intercept, covariate_arr])
    spline_idx = slice(1, 1 + n_spline)

    module_cols = [c for c in eigengenes.columns if c != "sample_id"]
    rows = []
    # Per-module F-test accumulator for later BH correction
    ftest_pvals: dict[str, float] = {}

    for col in module_cols:
        y = merged[col].to_numpy(dtype=float)
        mask = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
        n_obs = int(mask.sum())
        if n_obs < 10:
            continue
        X_fit, y_fit = X[mask], y[mask]
        coef, _, rank, _ = np.linalg.lstsq(X_fit, y_fit, rcond=None)
        residuals = y_fit - X_fit @ coef
        # Correct residual df: n_obs - rank (not rank - 1)
        df_res = max(n_obs - rank, 1)
        sigma2 = (residuals @ residuals) / df_res
        V = sigma2 * np.linalg.pinv(X_fit.T @ X_fit)

        beta_spline = coef[spline_idx]
        V_spline = V[spline_idx, :][:, spline_idx]

        # Persist raw coefficients + full covariance (repeated per module row) so the
        # trajectory can be re-projected to any age grid / fed to mash downstream.
        coef_cols = {f"coef_{i + 1}": float(beta_spline[i]) for i in range(n_spline)}
        cov_cols = {
            f"cov_age_{i + 1}{j + 1}": float(V_spline[i, j])
            for i in range(n_spline)
            for j in range(n_spline)
        }

        beta_proj = B_proj @ beta_spline
        V_proj = B_proj @ V_spline @ B_proj.T
        se_proj = np.sqrt(np.maximum(np.diag(V_proj), 0))
        z_proj = np.where(se_proj > 0, beta_proj / se_proj, 0.0)
        p_proj = 2 * stats.norm.sf(np.abs(z_proj))

        # F-test: full (spline + covariates) vs reduced (covariates only)
        X_red_fit = X_reduced[mask]
        coef_red, _, rank_red, _ = np.linalg.lstsq(X_red_fit, y_fit, rcond=None)
        rss_full = float(residuals @ residuals)
        rss_red = float(((y_fit - X_red_fit @ coef_red) ** 2).sum())
        df_num = max(rank - rank_red, 1)
        f_stat = ((rss_red - rss_full) / df_num) / (rss_full / df_res)
        p_ftest = float(stats.f.sf(f_stat, df_num, df_res)) if np.isfinite(f_stat) else np.nan
        ftest_pvals[col] = p_ftest

        for label, age_p, beta, se, z, p in zip(
            AGE_LABELS, AGE_PROBS, beta_proj, se_proj, z_proj, p_proj
        ):
            rows.append({
                "module_id": col,
                "trait": "Age_spline",
                "age_label": label,
                "age_prob": age_p,
                "effect": beta,
                "se": se,
                "z": z,
                "pvalue": p,
                "pvalue_ftest": p_ftest,
                "n": n_obs,
                **coef_cols,
                **cov_cols,
            })

    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = result.groupby("age_label")["pvalue"].transform(
            lambda pv: stats.false_discovery_control(pv, method="bh")
        )
        # Module-level FDR using the F-test p-value (one p-value per module)
        mod_ids = list(ftest_pvals.keys())
        fdr_vals = stats.false_discovery_control(
            [ftest_pvals[m] for m in mod_ids], method="bh"
        )
        fdr_ftest_map = dict(zip(mod_ids, fdr_vals))
        result["fdr_ftest"] = result["module_id"].map(fdr_ftest_map)
    return result


def diagnosis_association(
    eigengenes: pd.DataFrame, sample_table: pd.DataFrame, covariate_cols: list[str],
    dx_col: str = "Dx", control_label: str = "Control", case_label: str = "SCZD",
) -> pd.DataFrame:
    """Module eigengene association with case/control diagnosis."""
    if dx_col not in sample_table.columns:
        raise ValueError(f"Missing diagnosis column: {dx_col}")

    keep_cols = ["sample_id", dx_col] + [c for c in covariate_cols if c in sample_table.columns]
    merged = eigengenes.merge(sample_table[keep_cols], on="sample_id", how="inner")
    merged = merged.loc[merged[dx_col].isin([control_label, case_label])].copy()
    merged[dx_col] = pd.Categorical(merged[dx_col], categories=[control_label, case_label])
    if merged[dx_col].nunique(dropna=True) < 2:
        raise ValueError(f"Diagnosis association requires both {control_label} and {case_label} samples.")

    available_covariates = [c for c in covariate_cols if c in merged.columns]
    module_cols = [c for c in eigengenes.columns if c != "sample_id"]
    rows = []
    for col in module_cols:
        analysis = pd.DataFrame(
            {
                "eigengene": pd.to_numeric(merged[col], errors="coerce"),
                "Diagnosis": merged[dx_col],
            }
        )
        if available_covariates:
            analysis = pd.concat([analysis, merged[available_covariates].copy()], axis=1)
        analysis = analysis.replace([np.inf, -np.inf], np.nan).dropna()
        if len(analysis) < 10 or analysis["Diagnosis"].nunique(dropna=True) < 2:
            continue

        y = analysis["eigengene"].to_numpy(dtype=float)
        dx = (analysis["Diagnosis"] == case_label).astype(float).to_numpy()[:, None]
        covariate_df = pd.get_dummies(analysis[available_covariates], drop_first=True)
        covariate_arr = covariate_df.to_numpy(dtype=float) if not covariate_df.empty else np.empty((len(analysis), 0))
        X = np.hstack([np.ones((len(analysis), 1)), dx, covariate_arr])
        finite = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
        y_fit = y[finite]
        X_fit = X[finite]
        if len(y_fit) < 10 or np.linalg.matrix_rank(X_fit) < 2:
            continue

        coef, _, rank, _ = np.linalg.lstsq(X_fit, y_fit, rcond=None)
        residuals = y_fit - X_fit @ coef
        df_res = max(len(y_fit) - rank, 1)
        sigma2 = (residuals @ residuals) / df_res
        V = sigma2 * np.linalg.pinv(X_fit.T @ X_fit)
        se = float(np.sqrt(max(V[1, 1], 0.0)))
        effect = float(coef[1])
        t_stat = effect / se if se > 0 else np.nan
        pvalue = float(2 * stats.t.sf(abs(t_stat), df_res)) if np.isfinite(t_stat) else np.nan
        used_dx = analysis.loc[finite, "Diagnosis"]
        rows.append(
            {
                "module_id": col,
                "trait": f"{dx_col}_{case_label}_vs_{control_label}",
                "effect": effect,
                "se": se,
                "t": t_stat,
                "pvalue": pvalue,
                "n": int(len(y_fit)),
                "n_control": int((used_dx == control_label).sum()),
                "n_case": int((used_dx == case_label).sum()),
            }
        )

    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = stats.false_discovery_control(result["pvalue"], method="bh")
    return result


def _eigengenes_to_sample_table(artifacts) -> pd.DataFrame | None:
    eigengenes = artifacts.eigengene_table
    if eigengenes is None or eigengenes.empty:
        return None
    # eigengene_table: rows=modules, cols=[module_id, sample1, sample2, ...]
    # Pivot to: rows=samples, cols=[sample_id, M000, M001, ...]
    eg_pivot = eigengenes.set_index("module_id").T.reset_index().rename(columns={"index": "sample_id"})
    eg_pivot.columns.name = None
    return eg_pivot


def _save_core_artifacts(artifacts, out: Path) -> None:
    artifacts.module_table.to_parquet(out / "modules.parquet", index=False, compression="zstd")
    artifacts.edge_table.to_parquet(out / "edges.parquet", index=False, compression="zstd")
    artifacts.trait_table.to_parquet(out / "traits.parquet", index=False, compression="zstd")
    artifacts.feature_scores.to_parquet(out / "feature_scores.parquet", index=False, compression="zstd")
    pd.DataFrame([artifacts.calibration or {}]).to_parquet(out / "calibration.parquet", index=False, compression="zstd")
    if artifacts.module_gene_roles is not None and not artifacts.module_gene_roles.empty:
        artifacts.module_gene_roles.to_parquet(out / "module_gene_roles.parquet", index=False, compression="zstd")
    node_diag = getattr(artifacts, "node_diagnostics", None)
    if node_diag is not None and not node_diag.empty:
        node_diag.to_parquet(out / "node_diagnostics.parquet", index=False, compression="zstd")
    resid_qc = getattr(artifacts, "residualization_qc", None)
    if resid_qc is not None and not resid_qc.empty:
        resid_qc.to_parquet(out / "residualization_qc.parquet", index=False, compression="zstd")
    recon = getattr(artifacts, "feature_reconstruction", None)
    if recon is not None and not recon.empty:
        # VAE reconstruction of the multiplex feature matrix; enables post-hoc
        # re-projection of the gene graph under switch-only / switch-primary /
        # full-multiplex channel rules without re-fitting (see project_tiers).
        recon.to_parquet(out / "feature_reconstruction.parquet", index=False, compression="zstd")


# RNA-quality covariates aligned with the BrainSEQ DE/DTU aging model
# (limma/satuRn design: RIN + mito_mapping_rate + percent_assigned + SVA). The IsoGraph
# metrics parquet carries the both-ancestry analogs: exonic_rate (== featureCounts-style
# percent_assigned) and x3_bias_75th_percentile (3' coverage bias / degradation). These
# are the dominant confounds of the caudate aging modules (|corr|~0.73 with exonic_rate),
# stronger and more principled than median TIN, and -- unlike TIN / the AA-only DE
# phenotype.tsv -- available for every sample in all three BrainSEQ regions.
RNASEQC_QC_COVARIATES = ["exonic_rate", "x3_bias_75th_percentile"]


def _rnaseqc_covariate_table(region: str | None) -> pd.DataFrame | None:
    """DE-aligned RNA-quality covariates (exonic_rate, x3_bias_75th_percentile) for this
    BrainSEQ region keyed to sample_id, or None if the metrics parquet is absent."""
    if region is None:
        return None
    p = rel("inputs", "processed", "brainseq", "metadata", f"{region}_rnaseq_metrics.parquet")
    if not p.exists():
        return None
    df = pd.read_parquet(p, columns=["sample_rnum", *RNASEQC_QC_COVARIATES])
    df = df.rename(columns={"sample_rnum": "sample_id"})
    df["sample_id"] = df["sample_id"].astype(str)
    return df


# GTEx has no PEER factors in the bundle, but the GTEx RNAseQC suite (already in the
# sample_table) carries the direct analogs of the brainseq DE-aligned QC covariates:
# SMEXNCRT == exonic rate (percent_assigned analog), SM3PB75P == 3' bias 75th percentile
# (x3_bias_75th analog). Used as RNA-quality proxies in the GTEx aging association,
# matching the brainseq covariate philosophy cross-cohort.
GTEX_QC_COVARIATES = ["SMEXNCRT", "SM3PB75P"]

# Covariate adjustment is split (2026-06-28): the *_COVARIATES sets adjust trait
# inference downstream; the *_DISCOVERY_COVARIATES sets are the upstream residualization
# (discovery knob). Discovery keeps technical/topology confounds only (RNA quality,
# alignment, ancestry structure); biological Sex/MoD and the trait are NOT regressed out
# before clustering -- they are adjusted jointly with the trait in the association test.
BRAINSEQ_COVARIATES = [
    "Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
    "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5",
]
BRAINSEQ_DISCOVERY_COVARIATES = [
    "RIN", "mapping_rate", "mito_rate",
    "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5",
]
GTEX_COVARIATES = ["SEX", "SMRIN", "SMTSISCH", "SMMAPRT"]
GTEX_DISCOVERY_COVARIATES = ["SMRIN", "SMTSISCH", "SMMAPRT"]


def _gtex_qc_covariate_table(sample_table: pd.DataFrame) -> pd.DataFrame | None:
    """DE-aligned GTEx RNA-quality covariates pulled from the bundle sample_table
    (numeric), keyed to sample_id; None if the columns are absent."""
    cols = [c for c in GTEX_QC_COVARIATES if c in sample_table.columns]
    if not cols:
        return None
    out = pd.DataFrame({"sample_id": sample_table["sample_id"].astype(str)})
    for c in cols:
        out[c] = pd.to_numeric(sample_table[c], errors="coerce").to_numpy()
    return out


def _save_age_artifacts(artifacts, out, sample_table, covariate_cols, age_col, label,
                        qc_table: pd.DataFrame | None = None):
    _save_core_artifacts(artifacts, out)
    eigengenes = _eigengenes_to_sample_table(artifacts)
    if eigengenes is not None:
        linear = linear_age_association(eigengenes, sample_table, age_col=age_col)
        linear.to_parquet(out / "age_linear.parquet", index=False, compression="zstd")

        spline = spline_age_association(eigengenes, sample_table, covariate_cols, age_col=age_col)
        spline.to_parquet(out / "age_spline.parquet", index=False, compression="zstd")

        msg = (f"  {label}: {len(artifacts.module_table['module_id'].unique())} modules | "
               f"linear n={len(linear)} | spline n={len(spline)}")

        # DE-aligned QC-adjusted spline (the spline to report): keeps the baseline
        # age_spline.parquet so the degradation attenuation is visible before/after.
        if qc_table is not None:
            qc_cols = [c for c in qc_table.columns if c != "sample_id"]
            st_qc = sample_table.copy()
            st_qc["sample_id"] = st_qc["sample_id"].astype(str)
            # Drop any same-named columns first so native QC cols (GTEx SMxxx already in
            # sample_table) are not duplicated by the merge.
            st_qc = st_qc.drop(columns=[c for c in qc_cols if c in st_qc.columns])
            st_qc = st_qc.merge(qc_table, on="sample_id", how="left")
            spline_qc = spline_age_association(
                eigengenes, st_qc, covariate_cols + qc_cols, age_col=age_col,
            )
            spline_qc.to_parquet(out / "age_spline_qc_adjusted.parquet", index=False, compression="zstd")
            msg += f" | spline+QC({'+'.join(qc_cols)}) n={len(spline_qc)}"
        print(msg)


def _save_diagnosis_artifacts(
    artifacts, out: Path, bundle, covariate_cols: list[str], label: str,
    dx_col: str = "Dx", control_label: str = "Control", case_label: str = "SCZD",
) -> None:
    _save_core_artifacts(artifacts, out)
    for stale in ("age_linear.parquet", "age_spline.parquet"):
        path = out / stale
        if path.exists():
            path.unlink()

    eigengenes = _eigengenes_to_sample_table(artifacts)
    if eigengenes is not None:
        diagnosis = diagnosis_association(
            eigengenes,
            bundle.sample_table,
            covariate_cols=covariate_cols,
            dx_col=dx_col,
            control_label=control_label,
            case_label=case_label,
        )
        diagnosis.to_parquet(out / "diagnosis_assoc.parquet", index=False, compression="zstd")
        print(f"  {label}: {len(artifacts.module_table['module_id'].unique())} modules | "
              f"diagnosis n={len(diagnosis)}")


def run_brainseq_aging(regions: list[str] | None = None,
                       leiden_resolution: float | None = None) -> None:
    for region in regions or ["caudate", "hippocampus", "dlpfc"]:
        run_brainseq_region(region, leiden_resolution=leiden_resolution)


def run_gtex_aging(regions: list[str] | None = None,
                   leiden_resolution: float | None = None) -> None:
    gtex_bundle_root = rel("inputs", "bundles", "gtex_v11_brain")
    if regions is None:
        if gtex_bundle_root.exists():
            regions = sorted(path.name for path in gtex_bundle_root.iterdir() if path.is_dir())
        else:
            regions = GTEX_REGIONS
    for region in regions:
        run_gtex_region(region, leiden_resolution=leiden_resolution)




def run_brainseq_region(region: str, leiden_resolution: float | None = None) -> None:
    bundle = load_dataset_bundle(rel("inputs", "bundles", "brainseq_v1", region))
    sample_table = bundle.sample_table
    tc, tt = _filter_expressed_transcripts(
        bundle.matrices["transcript_counts"],
        bundle.feature_tables["transcript"],
    )
    del bundle  # free ~415 MB transcript_counts before VAE feature computation

    covariate_cols = BRAINSEQ_COVARIATES

    # full-multiplex is the PRIMARY production model (task #27): abundance-abundance
    # edges enabled with grid-calibrated alpha_abundance. The VAE fit is unchanged
    # (these settings only affect post-fit graph construction), so the saved
    # feature_reconstruction is identical and project_tiers still derives all three
    # tiers from it. Abundance acts as a conservative safety net on clean real data
    # (calibrated alpha ~0.90-0.95); see project_tiers / tier_checks.
    cfg = VaeModelConfig(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        residualize_covariates=BRAINSEQ_DISCOVERY_COVARIATES,
        min_module_size=20, trait_columns=["Age"], random_state=13,
        allow_abundance_abundance=True,
        alpha_switch=0.5,
        alpha_abundance_grid=[0.70, 0.75, 0.80, 0.85, 0.90, 0.95],
        leiden_resolution=CANONICAL_LEIDEN_RESOLUTION if leiden_resolution is None else leiden_resolution,
        **_PROMOTED_VAE,
    )
    print(f"[{region}] fitting model (res={cfg.leiden_resolution}) ...", flush=True)
    _t0 = time.time()
    artifacts = VaeNetworkModel(cfg).fit(
        transcript_counts=tc,
        transcript_table=tt,
        sample_table=sample_table,
    )
    del tc, tt  # free filtered transcript data; no longer needed after fit
    print(f"[{region}] fit done in {time.time() - _t0:.0f}s", flush=True)

    out = ensure_dir(rel("real_data", "brainseq", region, "_m",
                         _isograph_out_subdir(leiden_resolution)))
    _save_age_artifacts(artifacts, out, sample_table, covariate_cols, age_col="Age", label=region,
                        qc_table=_rnaseqc_covariate_table(region))


def run_brainseq_region_with_abundance(region: str, leiden_resolution: float | None = None) -> None:
    """Legacy comparison arm: same full-multiplex channels as the standard primary fit,
    but at the region's Part 1 BEST_LEIDEN_RESOLUTION instead of the canonical resolution.

    Writes artifacts to isograph_vae_with_abundance/ (separate from the standard run)
    so the best-resolution arm can be compared without overwriting the canonical modules.
    When leiden_resolution is None, uses the region's BEST_LEIDEN_RESOLUTION (Part 1 sweep).
    """
    if leiden_resolution is None:
        leiden_resolution = BEST_LEIDEN_RESOLUTION.get(region, 2.0)
    bundle = load_dataset_bundle(rel("inputs", "bundles", "brainseq_v1", region))
    sample_table = bundle.sample_table
    tc, tt = _filter_expressed_transcripts(
        bundle.matrices["transcript_counts"],
        bundle.feature_tables["transcript"],
    )
    del bundle

    covariate_cols = BRAINSEQ_COVARIATES

    cfg = VaeModelConfig(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        residualize_covariates=BRAINSEQ_DISCOVERY_COVARIATES,
        min_module_size=20, trait_columns=["Age"], random_state=13,
        allow_abundance_abundance=True,
        alpha_switch=0.5,
        alpha_abundance_grid=[0.70, 0.75, 0.80, 0.85, 0.90, 0.95],
        leiden_resolution=leiden_resolution,
        **_PROMOTED_VAE,
    )
    print(f"[{region}+abundance] fitting model ...", flush=True)
    _t0 = time.time()
    artifacts = VaeNetworkModel(cfg).fit(
        transcript_counts=tc,
        transcript_table=tt,
        sample_table=sample_table,
    )
    del tc, tt
    print(f"[{region}+abundance] fit done in {time.time() - _t0:.0f}s | "
          f"alpha_abundance={artifacts.calibration.get('alpha_abundance') if artifacts.calibration else 'n/a'}", flush=True)

    out = ensure_dir(rel("real_data", "brainseq", region, "_m", "isograph_vae_with_abundance"))
    _save_age_artifacts(artifacts, out, sample_table, covariate_cols, age_col="Age",
                        label=f"{region}+abundance", qc_table=_rnaseqc_covariate_table(region))


def run_gtex_region(region_dir_name: str, leiden_resolution: float | None = None) -> None:
    bundle = load_dataset_bundle(rel("inputs", "bundles", "gtex_v11_brain", region_dir_name))

    # GTEx QC covariates: RIN (SMRIN), ischemic time (SMTSISCH), mapping rate (SMMAPRT), sex (SEX)
    covariate_cols = GTEX_COVARIATES

    cfg = VaeModelConfig(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        # GTEx used to need a hand-tuned lr=3e-4: at the default lr=1e-3 the VAE
        # diverged (val ELBO -> ~1e8) on some regions. Validation gate B.2 (2026-06-22)
        # showed grad_clip_norm=1.0 (in _PROMOTED_VAE) lets the single default lr=1e-3
        # train every GTEx region — including nucleus_accumbens, which diverged pre-clip
        # — with no divergence/OOM. The per-cohort LR babysitting is therefore retired.
        residualize_covariates=GTEX_DISCOVERY_COVARIATES,
        min_module_size=20, trait_columns=["AGE"], random_state=13,
        # full-multiplex primary (task #27): abundance-abundance edges with grid
        # calibration; only affects post-fit graph, VAE reconstruction is unchanged.
        allow_abundance_abundance=True,
        alpha_switch=0.5,
        alpha_abundance_grid=[0.70, 0.75, 0.80, 0.85, 0.90, 0.95],
        leiden_resolution=CANONICAL_LEIDEN_RESOLUTION if leiden_resolution is None else leiden_resolution,
        **_PROMOTED_VAE,
    )
    print(f"[{region_dir_name}] fitting model (res={cfg.leiden_resolution}) ...", flush=True)
    _t0 = time.time()
    artifacts = VaeNetworkModel(cfg).fit(
        transcript_counts=bundle.matrices["transcript_counts"],
        transcript_table=bundle.feature_tables["transcript"],
        sample_table=bundle.sample_table,
    )
    print(f"[{region_dir_name}] fit done in {time.time() - _t0:.0f}s", flush=True)

    out = ensure_dir(rel("real_data", "gtex", region_dir_name, "_m",
                         _isograph_out_subdir(leiden_resolution)))
    _save_age_artifacts(artifacts, out, bundle.sample_table, covariate_cols, age_col="AGE",
                        label=region_dir_name, qc_table=_gtex_qc_covariate_table(bundle.sample_table))


def _drd2_gene_id(transcript_table: pd.DataFrame) -> str | None:
    """Resolve DRD2 ENSEMBL gene_id from transcript names (e.g. 'DRD2-201')."""
    if "transcript_name" not in transcript_table.columns:
        return None
    mask = transcript_table["transcript_name"].str.startswith("DRD2-", na=False)
    matches = transcript_table.loc[mask, "gene_id"]
    return str(matches.iloc[0]) if not matches.empty else None


def check_drd2(artifacts, transcript_table: pd.DataFrame) -> bool:
    """Return True if DRD2 is assigned to any module; print a pass/fail summary."""
    gene_id = _drd2_gene_id(transcript_table)
    if gene_id is None:
        print("WARNING: DRD2 check — could not resolve gene_id from transcript_name column")
        return False
    if artifacts.module_table.empty or "gene_id" not in artifacts.module_table.columns:
        print("WARNING: DRD2 check — module_table is empty or missing gene_id")
        return False
    hits = artifacts.module_table[artifacts.module_table["gene_id"] == gene_id]
    if hits.empty:
        print(f"WARNING: DRD2 check FAILED — {gene_id} not found in any module")
        return False
    module_ids = hits["module_id"].unique().tolist()
    role = None
    if artifacts.module_gene_roles is not None and not artifacts.module_gene_roles.empty:
        role_hits = artifacts.module_gene_roles[artifacts.module_gene_roles["gene_id"] == gene_id]
        if not role_hits.empty:
            role = role_hits["module_role"].iloc[0]
    print(f"DRD2 check PASSED — {gene_id} found in module(s): {module_ids} (role: {role})")
    return True


def run_brainseq_caudate_sczd(leiden_resolution: float | None = None) -> None:
    """Run IsoGraph VAE on the SCZD+Control caudate bundle (Dx as trait)."""
    bundle = load_dataset_bundle(rel("inputs", "bundles", "brainseq_sczd", "caudate"))

    covariate_cols = BRAINSEQ_COVARIATES

    # full-multiplex primary (task #27): unified with the aging canonical fits.
    # Abundance-abundance edges with grid calibration; only affects post-fit graph,
    # VAE reconstruction unchanged. Note: this can recover DRD2 (previously
    # isolated_below_alpha at the fixed 0.95) if the calibrated alpha admits its
    # 0.93 abundance edge — the DRD2 case study reflects the unified model.
    cfg = VaeModelConfig(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        residualize_covariates=BRAINSEQ_DISCOVERY_COVARIATES,
        min_module_size=20, trait_columns=["Dx"],
        random_state=13,
        allow_abundance_abundance=True,
        alpha_switch=0.5,
        alpha_abundance_grid=[0.70, 0.75, 0.80, 0.85, 0.90, 0.95],
        leiden_resolution=CANONICAL_LEIDEN_RESOLUTION if leiden_resolution is None else leiden_resolution,
        **_PROMOTED_VAE,
    )
    print(f"[caudate_sczd] fitting model (res={cfg.leiden_resolution}) ...", flush=True)
    _t0 = time.time()
    artifacts = VaeNetworkModel(cfg).fit(
        transcript_counts=bundle.matrices["transcript_counts"],
        transcript_table=bundle.feature_tables["transcript"],
        sample_table=bundle.sample_table,
    )
    print(f"[caudate_sczd] fit done in {time.time() - _t0:.0f}s", flush=True)

    out = ensure_dir(rel("real_data", "brainseq", "caudate_sczd", "_m",
                         _isograph_out_subdir(leiden_resolution)))
    diagnosis_covariates = ["Age"] + covariate_cols
    _save_diagnosis_artifacts(
        artifacts, out, bundle, covariate_cols=diagnosis_covariates, label="caudate_sczd",
    )
    check_drd2(artifacts, bundle.feature_tables["transcript"])


def run_brainseq_caudate_sczd_with_abundance(leiden_resolution: float | None = None) -> None:
    """Re-enable abundance-abundance edges for the SCZD+Control caudate bundle.

    Writes artifacts to isograph_vae_with_abundance/ (separate from the baseline)
    so the with-abundance partition can be compared without overwriting it. When
    leiden_resolution is None, uses BEST_LEIDEN_RESOLUTION["caudate_sczd"] (Part 1 sweep).
    """
    if leiden_resolution is None:
        leiden_resolution = BEST_LEIDEN_RESOLUTION.get("caudate_sczd", 2.0)
    bundle = load_dataset_bundle(rel("inputs", "bundles", "brainseq_sczd", "caudate"))

    covariate_cols = BRAINSEQ_COVARIATES

    cfg = VaeModelConfig(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        residualize_covariates=BRAINSEQ_DISCOVERY_COVARIATES,
        min_module_size=20, trait_columns=["Dx"],
        random_state=13,
        allow_abundance_abundance=True,
        alpha_switch=0.5,
        alpha_abundance_grid=[0.70, 0.75, 0.80, 0.85, 0.90, 0.95],
        leiden_resolution=leiden_resolution,
        **_PROMOTED_VAE,
    )
    print("[caudate_sczd+abundance] fitting model ...", flush=True)
    _t0 = time.time()
    artifacts = VaeNetworkModel(cfg).fit(
        transcript_counts=bundle.matrices["transcript_counts"],
        transcript_table=bundle.feature_tables["transcript"],
        sample_table=bundle.sample_table,
    )
    print(f"[caudate_sczd+abundance] fit done in {time.time() - _t0:.0f}s | "
          f"alpha_abundance={artifacts.calibration.get('alpha_abundance') if artifacts.calibration else 'n/a'}",
          flush=True)

    out = ensure_dir(rel("real_data", "brainseq", "caudate_sczd", "_m", "isograph_vae_with_abundance"))
    diagnosis_covariates = ["Age"] + covariate_cols
    _save_diagnosis_artifacts(
        artifacts, out, bundle, covariate_cols=diagnosis_covariates, label="caudate_sczd+abundance",
    )
    check_drd2(artifacts, bundle.feature_tables["transcript"])


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )
    parser = argparse.ArgumentParser(description="Run real-data IsoGraph analyses.")
    parser.add_argument(
        "analysis", nargs="?", default="all",
        choices=["all", "brainseq-aging", "brainseq-sczd", "gtex-aging"],
    )
    parser.add_argument(
        "--region", action="append",
        help="Run only the named region. Can be repeated for brainseq-aging or gtex-aging.",
    )
    parser.add_argument(
        "--variant", default="standard", choices=["standard", "with-abundance"],
        help="standard: full-multiplex primary config (allow_abundance_abundance=True, "
             "alpha_abundance_grid calibration) at the canonical resolution, written to "
             "isograph_vae. with-abundance: legacy comparison arm with the same channels "
             "at the Part 1 BEST_LEIDEN_RESOLUTION, written to isograph_vae_with_abundance. "
             "Applies to brainseq-aging and brainseq-sczd.",
    )
    parser.add_argument(
        "--leiden-resolution", type=float, default=None,
        help="Override leiden_resolution. Default (None) and the canonical value "
             f"({CANONICAL_LEIDEN_RESOLUTION}) both write to the canonical "
             "isograph_vae dir that the GWAS + trust-funnel cascades consume. "
             "Any OTHER value writes to a resolution-suffixed sibling dir (e.g. "
             "isograph_vae_res2 for 2.0, isograph_vae_res2p25 for 2.25) so a "
             "biology-driven resolution sweep is a set of side-by-side comparison "
             "dirs that never clobber the canonical modules (with-abundance uses "
             "BEST_LEIDEN_RESOLUTION from the Part 1 sweep).",
    )
    args = parser.parse_args()

    if args.analysis == "brainseq-aging":
        if args.variant == "with-abundance":
            for region in (args.region or ["caudate", "hippocampus", "dlpfc"]):
                run_brainseq_region_with_abundance(region, leiden_resolution=args.leiden_resolution)
        else:
            run_brainseq_aging(args.region, leiden_resolution=args.leiden_resolution)
    elif args.analysis == "brainseq-sczd":
        if args.variant == "with-abundance":
            run_brainseq_caudate_sczd_with_abundance(leiden_resolution=args.leiden_resolution)
        else:
            run_brainseq_caudate_sczd(leiden_resolution=args.leiden_resolution)
    elif args.analysis == "gtex-aging":
        run_gtex_aging(args.region, leiden_resolution=args.leiden_resolution)
    else:
        print("BrainSEQ aging regions:")
        run_brainseq_aging()

        print("\nBrainSEQ caudate SCZD+Control:")
        run_brainseq_caudate_sczd()

        print("\nGTEx v11 brain aging regions:")
        run_gtex_aging()


if __name__ == "__main__":
    main()
