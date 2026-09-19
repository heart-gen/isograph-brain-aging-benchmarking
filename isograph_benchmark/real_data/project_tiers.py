"""Project the multiplex gene graph into the 4-tier ablation hierarchy from ONE VAE fit.

The VAE denoises the same switch+abundance feature matrix regardless of which edge
channels we later admit, so the three IsoGraph tiers differ *only* in edge-projection
rules and can be derived from a single fit's saved reconstruction -- a rigorous
ablation (the latent representation is held constant) that is also ~3x cheaper than
re-fitting per tier. Reuses ``feature_reconstruction.parquet`` (and feature_scores)
without retraining, mirroring ``sweep_leiden``'s reuse-saved-artifacts pattern.

Tiers (IsoGraph):
    switch_only      -- switch-switch edges only           (pure isoform-switch ablation)
    switch_primary   -- + cross to abundance-only genes     (switch-focused ablation)
    full_multiplex   -- + calibrated abundance-abundance     (PRIMARY production model)
WGCNA is the separate gene-abundance baseline (run elsewhere).

Promotion rationale: under synthetic 3' RNA degradation the switch-primary model's
module recovery collapses (~0.97 clean -> ~0.20 coupled-degradation) while the
abundance-abundance channel rescues it to near-WGCNA levels (~0.68 vs WGCNA ~0.76).
The full multiplex is therefore the robust-across-conditions production model; the
switch-restricted tiers are retained as sensitivity analyses.

Usage (from project root, isograph env):
    python -m isograph_benchmark.real_data.project_tiers brainseq-aging --region caudate
    python -m isograph_benchmark.real_data.project_tiers brainseq-sczd
    python -m isograph_benchmark.real_data.project_tiers gtex-aging --region caudate_basal_ganglia
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from isograph.features.reliability import gene_switch_estimability
from isograph.io.artifacts import load_dataset_bundle
from isograph.models.base import compute_module_gene_roles, compute_trait_associations
from isograph.models.multiplex import (
    project_feature_similarity_to_gene_graph,
    reconstruction_to_similarity,
    select_alpha_abundance,
)
from isograph.models.vae import VaeNetworkModel, _build_node_diagnostics
from isograph.workflow.config import VaeModelConfig
from isograph_benchmark.paths import analysis_store, ensure_dir, rel
from isograph_benchmark.real_data.run_models import (
    _PROMOTED_VAE,
    filter_production_transcripts,
    _gtex_qc_covariate_table,
    _rnaseqc_covariate_table,
    diagnosis_association,
    linear_age_association,
    spline_age_association,
)

META_COLS = ["feature_id", "gene_id", "feature_type", "n_transcripts"]

# Tier -> projection-rule overrides for project_feature_similarity_to_gene_graph.
TIERS: dict[str, dict] = {
    "switch_only": dict(switch_only=True, allow_abundance_abundance=False),
    "switch_primary": dict(switch_only=False, allow_abundance_abundance=False),
    "full_multiplex": dict(switch_only=False, allow_abundance_abundance=True),
}
PRIMARY_TIER = "full_multiplex"
# Output subdir per tier. The canonical 'isograph_vae' dir is the full_multiplex primary
# (run_models); these per-tier dirs are the apples-to-apples ablation replicas projected
# from the same fit, so writing them never clobbers the canonical cascade inputs.
TIER_DIRS: dict[str, str] = {
    "switch_only": "isograph_vae_switch_only",
    "switch_primary": "isograph_vae_switch_primary",
    "full_multiplex": "isograph_vae_full_multiplex",
}
# alpha_abundance calibration grid for the full-multiplex abundance-abundance channel
# (mirrors run_brainseq_region_with_abundance).
ALPHA_ABUNDANCE_GRID = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]


def _bundle_path(analysis: str, region: str | None) -> Path:
    if analysis == "brainseq-sczd":
        return rel("inputs", "bundles", "brainseq_sczd", "caudate")
    if analysis == "brainseq-aging":
        assert region is not None
        return rel("inputs", "bundles", "brainseq_v1", region)
    if analysis == "gtex-aging":
        assert region is not None
        return rel("inputs", "bundles", "gtex_v11_brain", region)
    raise ValueError(f"Unknown analysis: {analysis!r}")


def _region_dir(analysis: str, region: str | None) -> Path:
    return analysis_store(analysis, region)


def _trait_spec(analysis: str) -> tuple[str, list[str]]:
    """Return (trait_col, covariate_cols) matching the production run for this analysis."""
    if analysis == "brainseq-sczd":
        return "Dx", [
            "Age", "Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
            "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5",
        ]
    if analysis == "brainseq-aging":
        return "Age", [
            "Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
            "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5",
        ]
    if analysis == "gtex-aging":
        return "AGE", ["SEX", "SMRIN", "SMTSISCH", "SMMAPRT"]
    raise ValueError(f"Unknown analysis: {analysis!r}")


def _discovery_covariates(analysis: str) -> list[str]:
    """Upstream residualization (discovery knob), mirroring run_models.*_DISCOVERY_COVARIATES:
    technical/topology confounds only. Biological covariates and the trait are adjusted
    downstream in _associate, not regressed out before clustering."""
    if analysis in ("brainseq-sczd", "brainseq-aging"):
        return ["RIN", "mapping_rate", "mito_rate",
                "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5"]
    if analysis == "gtex-aging":
        return ["SMRIN", "SMTSISCH", "SMMAPRT"]
    raise ValueError(f"Unknown analysis: {analysis!r}")


def _pivot_eigengenes(eigengene_table: pd.DataFrame) -> pd.DataFrame:
    return (
        eigengene_table.set_index("module_id").T.reset_index()
        .rename(columns={"index": "sample_id"}).rename_axis(None, axis=1)
    )


def _base_cfg(analysis: str, discovery_covariates: list[str], trait_col: str) -> VaeModelConfig:
    """Config carrying the production module-detection knobs (leiden_resolution,
    min_module_size, random_state, estimability) so tier modules match production."""
    return VaeModelConfig(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        residualize_covariates=discovery_covariates,
        min_module_size=20, trait_columns=[trait_col], random_state=13,
        allow_abundance_abundance=False, alpha_switch=0.5, leiden_resolution=2.0,
        **_PROMOTED_VAE,
    )


def _associate(analysis, module_table, feature_scores, sample_table, covariate_cols, trait_col,
               qc_table=None):
    """Trait association for the tier's modules; returns (assoc_df, label, n_sig, extras)
    where extras is a list of (label, df) for additional associations (e.g. the DE-aligned
    QC-adjusted aging spline, mirroring run_models._save_age_artifacts)."""
    trait_table, eigengene_table = compute_trait_associations(
        module_table, feature_scores, sample_table, trait_columns=[trait_col],
    )
    if eigengene_table.empty:
        return pd.DataFrame(), "none", 0, []
    eg = _pivot_eigengenes(eigengene_table)
    if analysis == "brainseq-sczd":
        assoc = diagnosis_association(eg, sample_table, covariate_cols=covariate_cols)
        n_sig = int((assoc["fdr"] <= 0.10).sum()) if not assoc.empty and "fdr" in assoc else 0
        return assoc, "diagnosis_assoc", n_sig, []
    age_col = "AGE" if analysis == "gtex-aging" else "Age"
    spline = spline_age_association(eg, sample_table, covariate_cols=covariate_cols, age_col=age_col)
    n_sig = (
        int(spline[spline["fdr_ftest"] <= 0.10]["module_id"].nunique())
        if not spline.empty and "fdr_ftest" in spline else 0
    )
    extras = []
    if qc_table is not None:
        qc_cols = [c for c in qc_table.columns if c != "sample_id"]
        st_qc = sample_table.copy()
        st_qc["sample_id"] = st_qc["sample_id"].astype(str)
        # Collision-safe: GTEx QC cols (SMEXNCRT/SM3PB75P) are native to the bundle
        # sample_table, so drop them before merging or the join suffixes to _x/_y
        # and the qc covariates silently vanish from the spline design.
        st_qc = st_qc.drop(columns=[c for c in qc_cols if c in st_qc.columns])
        st_qc = st_qc.merge(qc_table, on="sample_id", how="left")
        spline_qc = spline_age_association(
            eg, st_qc, covariate_cols=covariate_cols + qc_cols, age_col=age_col,
        )
        extras.append(("age_spline_qc_adjusted", spline_qc))
    return spline, "age_spline", n_sig, extras


def project_region(analysis: str, region: str | None, source_subdir: str = "isograph_vae",
                   dry_run: bool = False) -> pd.DataFrame:
    label = f"{analysis}/{region}" if region else analysis
    rdir = _region_dir(analysis, region)
    src = rdir / source_subdir
    recon_path = src / "feature_reconstruction.parquet"
    fs_path = src / "feature_scores.parquet"
    if not recon_path.exists():
        print(f"[{label}] ERROR: missing {recon_path}. Re-run the fit to emit the reconstruction.")
        return pd.DataFrame()

    print(f"[{label}] loading reconstruction + feature_scores ...", flush=True)
    recon = pd.read_parquet(recon_path)
    feature_info = recon[META_COLS].reset_index(drop=True)
    sample_cols = [c for c in recon.columns if c not in META_COLS]
    x_recon = recon[sample_cols].to_numpy(dtype=np.float32).T  # (n_samples, n_features)
    feature_scores = pd.read_parquet(fs_path)

    trait_col, covariate_cols = _trait_spec(analysis)
    bundle = load_dataset_bundle(_bundle_path(analysis, region))
    sample_table = bundle.sample_table

    # Per-gene estimability reliability (promoted production lever), recomputed from
    # the bundle so switch-edge downweighting matches the production fit.
    tc, tt = filter_production_transcripts(
        bundle.matrices["transcript_counts"], bundle.feature_tables["transcript"],
    )
    gene_reliability = gene_switch_estimability(
        tc, tt, min_minor_usage=_PROMOTED_VAE["switch_estimability_min_minor_usage"],
    )
    del tc, tt, bundle

    print(f"[{label}] recomputing feature similarity ({x_recon.shape[1]} features) ...", flush=True)
    t0 = time.time()
    sim = reconstruction_to_similarity(x_recon)
    print(f"  similarity in {time.time()-t0:.0f}s", flush=True)

    cfg = _base_cfg(analysis, _discovery_covariates(analysis), trait_col)
    model = VaeNetworkModel(cfg)
    # DE-aligned QC covariates for the QC-adjusted aging spline: brainseq pulls from the
    # metrics parquet; GTEx pulls the native RNAseQC analogs (SMEXNCRT, SM3PB75P) from the
    # bundle sample_table. None for sczd (diagnosis association, not aging).
    if analysis == "brainseq-aging":
        qc_table = _rnaseqc_covariate_table(region)
    elif analysis == "gtex-aging":
        qc_table = _gtex_qc_covariate_table(sample_table)
    else:
        qc_table = None
    summary_rows = []

    for tier, rules in TIERS.items():
        t1 = time.time()
        alpha_abundance = None
        if tier == "full_multiplex":
            # Calibrate the abundance-abundance threshold to avoid a giant component.
            alpha_abundance = select_alpha_abundance(
                sim, feature_info, cfg.alpha, ALPHA_ABUNDANCE_GRID, alpha_switch=cfg.alpha_switch,
            )
        node_stats: dict = {}
        net_graph, edge_rows = project_feature_similarity_to_gene_graph(
            sim, feature_info, cfg.alpha,
            alpha_switch=cfg.alpha_switch,
            alpha_abundance=alpha_abundance,
            gene_reliability=gene_reliability,
            node_stats=node_stats,
            **rules,
        )
        module_table = model._module_table(net_graph)
        roles = compute_module_gene_roles(module_table, feature_scores, sample_table)
        node_diag = _build_node_diagnostics(
            feature_info=feature_info, edge_rows=edge_rows, module_table=module_table,
            gene_reliability=gene_reliability, node_stats=node_stats,
            alpha_switch=cfg.alpha_switch, min_module_size=cfg.min_module_size,
            reliability_on=cfg.switch_reliability_weighting,
        )
        assoc, assoc_label, n_sig, extra_assocs = _associate(
            analysis, module_table, feature_scores, sample_table, covariate_cols, trait_col,
            qc_table=qc_table,
        )

        n_modules = int(module_table["module_id"].nunique()) if not module_table.empty else 0
        sizes = (module_table.groupby("module_id").size().sort_values(ascending=False)
                 if not module_table.empty else pd.Series(dtype=int))
        giant = int(sizes.iloc[0]) if len(sizes) else 0
        n_assigned = int(len(module_table))
        summary_rows.append({
            "analysis": analysis, "region": region or "caudate_sczd", "tier": tier,
            "is_primary": tier == PRIMARY_TIER, "n_modules": n_modules,
            "n_assigned": n_assigned, "giant_size": giant,
            "giant_fraction": round(giant / n_assigned, 4) if n_assigned else 0.0,
            "n_edges": len(edge_rows), "alpha_abundance": alpha_abundance,
            "n_sig_trait_fdr10": n_sig, "elapsed_s": round(time.time() - t1, 1),
        })
        print(f"  [{tier}] modules={n_modules} giant={giant} "
              f"({summary_rows[-1]['giant_fraction']:.1%}) edges={len(edge_rows)} "
              f"alpha_ab={alpha_abundance} n_sig={n_sig} | {summary_rows[-1]['elapsed_s']:.0f}s",
              flush=True)

        if not dry_run:
            out = ensure_dir(rdir / TIER_DIRS[tier])
            module_table.to_parquet(out / "modules.parquet", index=False, compression="zstd")
            pd.DataFrame(edge_rows).to_parquet(out / "edges.parquet", index=False, compression="zstd")
            if not roles.empty:
                roles.to_parquet(out / "module_gene_roles.parquet", index=False, compression="zstd")
            if not node_diag.empty:
                node_diag.to_parquet(out / "node_diagnostics.parquet", index=False, compression="zstd")
            if not assoc.empty:
                assoc.to_parquet(out / f"{assoc_label}.parquet", index=False, compression="zstd")
            for xlabel, xdf in extra_assocs:
                if not xdf.empty:
                    xdf.to_parquet(out / f"{xlabel}.parquet", index=False, compression="zstd")

    summary = pd.DataFrame(summary_rows)
    print(f"\n[{label}] tier summary:")
    print(summary[["tier", "is_primary", "n_modules", "giant_fraction",
                   "n_edges", "alpha_abundance", "n_sig_trait_fdr10"]].to_string(index=False))
    if not dry_run:
        summary.to_parquet(rdir / "tier_projection_summary.parquet", index=False, compression="zstd")
    return summary


def main() -> None:
    p = argparse.ArgumentParser(description="Project one VAE fit into the 3 IsoGraph tiers.")
    p.add_argument("analysis", choices=["brainseq-sczd", "brainseq-aging", "gtex-aging"])
    p.add_argument("--region", default=None, help="Region (required for aging analyses).")
    p.add_argument("--source-subdir", default="isograph_vae",
                   help="Fit dir holding feature_reconstruction.parquet (default: isograph_vae).")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    project_region(args.analysis, args.region, source_subdir=args.source_subdir, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
