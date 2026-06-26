"""A/B test: which RNA-quality covariate set should the IsoGraph aging association use?

For each module we run the production spline age-association under a panel of covariate
sets and compare FDR(F-test) at alpha. The BrainSEQ DE/DTU aging model (limma/satuRn)
uses RIN + mito_mapping_rate + percent_assigned + SVA -- NOT median TIN -- so we test
the both-ancestry analogs from the IsoGraph metrics parquet:
    exonic_rate              == featureCounts-style percent_assigned (assigned fraction)
    x3_bias_75th_percentile  == 3' coverage bias / degradation
and compare against median_TIN (AA-only, indirect; over-aggressive).

Finding (caudate): the read-distribution / 3'-bias metrics correlate ~0.7 with the
significant aging modules and attenuate them to borderline (FDR ~0.06-0.08); median_TIN
annihilates them (FDR ~0.69). The adopted set is the DE-aligned ``exonic_rate +
x3_bias_75th_percentile`` (see run_models.RNASEQC_QC_COVARIATES).

Usage (project root, isograph env):
    python -m isograph_benchmark.real_data.qc_covariate_test brainseq-aging --region caudate
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from isograph.io.artifacts import load_dataset_bundle
from isograph.models.base import compute_trait_associations
from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.project_tiers import (
    TIER_DIRS,
    _bundle_path,
    _pivot_eigengenes,
    _region_dir,
    _trait_spec,
)
from isograph_benchmark.real_data.run_models import (
    RNASEQC_QC_COVARIATES,
    _rnaseqc_covariate_table,
    spline_age_association,
)

FDR_ALPHA = 0.10

# Covariate-set panel to compare (label -> extra covariate columns beyond baseline).
# "adopted" is the production choice wired into run_models.
PANEL: dict[str, list[str]] = {
    "baseline": [],
    "+exonic_rate": ["exonic_rate"],
    "+x3_bias_75th": ["x3_bias_75th_percentile"],
    "adopted (exonic+x3bias)": list(RNASEQC_QC_COVARIATES),
    "+median_TIN (compare)": ["median_TIN"],
}


def _tin_table(region: str | None) -> pd.DataFrame | None:
    if region is None:
        return None
    p = rel("inputs", "tin", f"brainseq__{region}__sample_median_tin.csv")
    if not p.exists():
        return None
    df = pd.read_csv(p)
    scol = "sample_id" if "sample_id" in df.columns else df.columns[0]
    vcol = [c for c in df.columns if c != scol][0]
    return pd.DataFrame({"sample_id": df[scol].astype(str), "median_TIN": df[vcol].astype(float)})


def _augmented_sample_table(region: str | None, bundle) -> pd.DataFrame:
    """Bundle sample_table + DE-aligned QC metrics + median_TIN (for comparison)."""
    st = bundle.sample_table.copy()
    st["sample_id"] = st["sample_id"].astype(str)
    qc = _rnaseqc_covariate_table(region)
    if qc is not None:
        st = st.merge(qc, on="sample_id", how="left")
    tin = _tin_table(region)
    if tin is not None:
        st = st.merge(tin, on="sample_id", how="left")
    return st


def _eligible(extra: list[str], st: pd.DataFrame) -> bool:
    return all(c in st.columns and st[c].notna().any() for c in extra)


def test_region(analysis: str, region: str | None) -> dict | None:
    label = f"{analysis}/{region}" if region else analysis
    if analysis != "brainseq-aging":
        print(f"[{label}] QC-covariate aging test currently targets brainseq-aging only.")
        return None

    rdir = _region_dir(analysis, region)
    _, base_covs = _trait_spec(analysis)
    bundle = load_dataset_bundle(_bundle_path(analysis, region))
    st = _augmented_sample_table(region, bundle)

    sources = {"isograph_vae (production)": "isograph_vae",
               **{f"tier:{t}": d for t, d in TIER_DIRS.items()}}
    rows, summary = [], {}
    for sl, sub in sources.items():
        sdir = rdir / sub
        if not (sdir / "modules.parquet").exists():
            continue
        fs_dir = sdir if (sdir / "feature_scores.parquet").exists() else rdir / "isograph_vae"
        if not (fs_dir / "feature_scores.parquet").exists():
            continue
        fs = pd.read_parquet(fs_dir / "feature_scores.parquet")
        modules = pd.read_parquet(sdir / "modules.parquet")

        _, eg = compute_trait_associations(modules, fs, st, trait_columns=["Age"])
        if eg.empty:
            continue
        egp = _pivot_eigengenes(eg)

        per_source = {}
        for plabel, extra in PANEL.items():
            if extra and not _eligible(extra, st):
                continue
            d = spline_age_association(egp, st, covariate_cols=base_covs + extra, age_col="Age")
            if d.empty or "fdr_ftest" not in d:
                continue
            fdr = d.drop_duplicates("module_id").set_index("module_id")["fdr_ftest"]
            nsig = int((fdr <= FDR_ALPHA).sum())
            per_source[plabel] = nsig
            for mid, v in fdr.items():
                rows.append({"source": sl, "panel": plabel, "module_id": mid,
                             "fdr_ftest": round(float(v), 4), "sig": bool(v <= FDR_ALPHA)})
        summary[sl] = per_source
        cells = "  ".join(f"{k}={v}" for k, v in per_source.items())
        print(f"[{label}] {sl}: n_sig  {cells}", flush=True)

    if not rows:
        print(f"[{label}] no eligible module tables / covariates.")
        return None

    out = ensure_dir(rdir / "qc_covariate_test")
    pd.DataFrame(rows).to_parquet(out / "module_age_qc_panel.parquet", index=False, compression="zstd")
    result = {"analysis": analysis, "region": region, "fdr_alpha": FDR_ALPHA,
              "adopted_covariates": list(RNASEQC_QC_COVARIATES), "n_sig_by_source_panel": summary}
    (out / "summary.json").write_text(json.dumps(result, indent=2))
    print(f"\n[{label}] adopted covariates: {RNASEQC_QC_COVARIATES} -> {out}")
    return result


def main() -> None:
    p = argparse.ArgumentParser(description="RNA-quality covariate-set A/B for module-aging associations.")
    p.add_argument("analysis", choices=["brainseq-aging"])
    p.add_argument("--region", default=None)
    args = p.parse_args()
    test_region(args.analysis, args.region)


if __name__ == "__main__":
    main()
