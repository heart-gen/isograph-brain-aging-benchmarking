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
AGE_PROBS = np.array([0.10, 0.25, 0.50, 0.75, 0.90])
AGE_LABELS = ["p10", "p25", "p50", "p75", "p90"]

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
    Non-linear age association via natural cubic spline, projected to 5 age points.

    Fits: eigengene ~ ns(age_z, knots=[q1/3, q2/3]) + covariates
    Projects spline coefficients to age_eval = qnorm([0.10, 0.25, 0.50, 0.75, 0.90]).
    """
    keep_cols = ["sample_id", age_col] + [c for c in covariate_cols if c in sample_table.columns]
    merged = eigengenes.merge(sample_table[keep_cols], on="sample_id", how="inner")
    merged = merged.dropna(subset=[age_col])

    age_z = _standardize(merged[age_col].to_numpy(dtype=float))
    ns_knots = np.quantile(age_z, [1 / 3, 2 / 3])
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
    spline_idx = slice(1, 1 + n_spline)

    module_cols = [c for c in eigengenes.columns if c != "sample_id"]
    rows = []

    for col in module_cols:
        y = merged[col].to_numpy(dtype=float)
        mask = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
        if mask.sum() < 10:
            continue
        X_fit, y_fit = X[mask], y[mask]
        coef, _, rank, _ = np.linalg.lstsq(X_fit, y_fit, rcond=None)
        residuals = y_fit - X_fit @ coef
        df_res = max(rank - 1, 1)
        sigma2 = (residuals @ residuals) / df_res
        V = sigma2 * np.linalg.pinv(X_fit.T @ X_fit)

        beta_spline = coef[spline_idx]
        V_spline = V[spline_idx, :][:, spline_idx]

        beta_proj = B_proj @ beta_spline
        V_proj = B_proj @ V_spline @ B_proj.T
        se_proj = np.sqrt(np.maximum(np.diag(V_proj), 0))
        z_proj = np.where(se_proj > 0, beta_proj / se_proj, 0.0)
        p_proj = 2 * stats.norm.sf(np.abs(z_proj))

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
                "n": int(mask.sum()),
            })

    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = result.groupby("age_label")["pvalue"].transform(
            lambda pv: stats.false_discovery_control(pv, method="bh")
        )
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
    eg_long = eigengenes.T.reset_index().rename(columns={"index": "sample_id"})
    eg_long.columns.name = None
    return eg_long


def _save_core_artifacts(artifacts, out: Path) -> None:
    artifacts.module_table.to_parquet(out / "modules.parquet", index=False, compression="zstd")
    artifacts.edge_table.to_parquet(out / "edges.parquet", index=False, compression="zstd")
    artifacts.trait_table.to_parquet(out / "traits.parquet", index=False, compression="zstd")
    artifacts.feature_scores.to_parquet(out / "feature_scores.parquet", index=False, compression="zstd")
    pd.DataFrame([artifacts.calibration or {}]).to_parquet(out / "calibration.parquet", index=False, compression="zstd")
    if artifacts.module_gene_roles is not None and not artifacts.module_gene_roles.empty:
        artifacts.module_gene_roles.to_parquet(out / "module_gene_roles.parquet", index=False, compression="zstd")


def _save_age_artifacts(artifacts, out, bundle, covariate_cols, age_col, label):
    _save_core_artifacts(artifacts, out)
    eigengenes = _eigengenes_to_sample_table(artifacts)
    if eigengenes is not None:
        linear = linear_age_association(eigengenes, bundle.sample_table, age_col=age_col)
        linear.to_parquet(out / "age_linear.parquet", index=False, compression="zstd")

        spline = spline_age_association(eigengenes, bundle.sample_table, covariate_cols, age_col=age_col)
        spline.to_parquet(out / "age_spline.parquet", index=False, compression="zstd")

        print(f"  {label}: {len(artifacts.module_table['module_id'].unique())} modules | "
              f"linear n={len(linear)} | spline n={len(spline)}")


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


def run_brainseq_aging(regions: list[str] | None = None) -> None:
    for region in regions or ["caudate", "hippocampus", "dlpfc"]:
        run_brainseq_region(region)


def run_gtex_aging(regions: list[str] | None = None) -> None:
    gtex_bundle_root = rel("inputs", "bundles", "gtex_v11_brain")
    if regions is None:
        if gtex_bundle_root.exists():
            regions = sorted(path.name for path in gtex_bundle_root.iterdir() if path.is_dir())
        else:
            regions = GTEX_REGIONS
    for region in regions:
        run_gtex_region(region)


_MULTIPLEX_ABUNDANCE_GRID = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]
_MULTIPLEX_SWITCH_GRID = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]


def run_brainseq_region(region: str) -> None:
    bundle = load_dataset_bundle(rel("inputs", "bundles", "brainseq_v1", region))

    covariate_cols = [
        "Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
        "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5",
    ]

    cfg = VaeModelConfig(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        residualize_covariates=covariate_cols,
        min_module_size=30, trait_columns=["Age"], random_state=13,
        allow_abundance_abundance=True,
        alpha_switch_grid=_MULTIPLEX_SWITCH_GRID,
        alpha_abundance_grid=_MULTIPLEX_ABUNDANCE_GRID,
    )
    print(f"[{region}] fitting model ...", flush=True)
    _t0 = time.time()
    artifacts = VaeNetworkModel(cfg).fit(
        transcript_counts=bundle.matrices["transcript_counts"],
        transcript_table=bundle.feature_tables["transcript"],
        sample_table=bundle.sample_table,
    )
    print(f"[{region}] fit done in {time.time() - _t0:.0f}s", flush=True)

    out = ensure_dir(rel("real_data", "brainseq", region, "_m", "isograph_vae"))
    _save_age_artifacts(artifacts, out, bundle, covariate_cols, age_col="Age", label=region)


def run_gtex_region(region_dir_name: str) -> None:
    bundle = load_dataset_bundle(rel("inputs", "bundles", "gtex_v11_brain", region_dir_name))

    # GTEx QC covariates: RIN (SMRIN), ischemic time (SMTSISCH), mapping rate (SMMAPRT), sex (SEX)
    covariate_cols = ["SEX", "SMRIN", "SMTSISCH", "SMMAPRT"]

    cfg = VaeModelConfig(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        residualize_covariates=covariate_cols,
        min_module_size=30, trait_columns=["AGE"], random_state=13,
        allow_abundance_abundance=True,
        alpha_switch_grid=_MULTIPLEX_SWITCH_GRID,
        alpha_abundance_grid=_MULTIPLEX_ABUNDANCE_GRID,
    )
    print(f"[{region_dir_name}] fitting model ...", flush=True)
    _t0 = time.time()
    artifacts = VaeNetworkModel(cfg).fit(
        transcript_counts=bundle.matrices["transcript_counts"],
        transcript_table=bundle.feature_tables["transcript"],
        sample_table=bundle.sample_table,
    )
    print(f"[{region_dir_name}] fit done in {time.time() - _t0:.0f}s", flush=True)

    out = ensure_dir(rel("real_data", "gtex", region_dir_name, "_m", "isograph_vae"))
    _save_age_artifacts(artifacts, out, bundle, covariate_cols, age_col="AGE", label=region_dir_name)


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


def run_brainseq_caudate_sczd() -> None:
    """Run IsoGraph VAE on the SCZD+Control caudate bundle (Dx as trait)."""
    bundle = load_dataset_bundle(rel("inputs", "bundles", "brainseq_sczd", "caudate"))

    covariate_cols = [
        "Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
        "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5",
    ]

    cfg = VaeModelConfig(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        residualize_covariates=covariate_cols,
        min_module_size=30, trait_columns=["Dx"],
        random_state=13,
        allow_abundance_abundance=True,
        alpha_switch_grid=_MULTIPLEX_SWITCH_GRID,
        alpha_abundance_grid=_MULTIPLEX_ABUNDANCE_GRID,
    )
    print("[caudate_sczd] fitting model ...", flush=True)
    _t0 = time.time()
    artifacts = VaeNetworkModel(cfg).fit(
        transcript_counts=bundle.matrices["transcript_counts"],
        transcript_table=bundle.feature_tables["transcript"],
        sample_table=bundle.sample_table,
    )
    print(f"[caudate_sczd] fit done in {time.time() - _t0:.0f}s", flush=True)

    out = ensure_dir(rel("real_data", "brainseq", "caudate_sczd", "_m", "isograph_vae"))
    diagnosis_covariates = ["Age"] + covariate_cols
    _save_diagnosis_artifacts(
        artifacts, out, bundle, covariate_cols=diagnosis_covariates, label="caudate_sczd",
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
    args = parser.parse_args()

    if args.analysis == "brainseq-aging":
        run_brainseq_aging(args.region)
    elif args.analysis == "brainseq-sczd":
        run_brainseq_caudate_sczd()
    elif args.analysis == "gtex-aging":
        run_gtex_aging(args.region)
    else:
        print("BrainSEQ aging regions:")
        run_brainseq_aging()

        print("\nBrainSEQ caudate SCZD+Control:")
        run_brainseq_caudate_sczd()

        print("\nGTEx v11 brain aging regions:")
        run_gtex_aging()


if __name__ == "__main__":
    main()
