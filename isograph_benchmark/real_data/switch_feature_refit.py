"""Full refit per preprocessing setting: does the *network* survive the feature choices?

``switch_feature_sensitivity`` rebuilds the switch channel under each pseudocount, expression
filter and minor-isoform setting but holds the published module partition fixed, so it
measures whether the representation and its trait signal are stable -- not whether an
independently refit network would find the same modules. This is that follow-on: for each
setting, the production IsoGraph fit is re-run from the perturbed transcript counts with the
production configuration (same VAE, covariates, calibration grid, Leiden resolution and seed),
and the resulting partition and module--age associations are compared with the published fit.

**The noise floor is measured, not assumed.** The published setting is refit too. Its agreement
with the committed partition is the ceiling any perturbed setting can reach, so a perturbed ARI
is read against it rather than against 1.

One fit per process (``--index``), because a real-data fit peaks near 30 GB; ``--aggregate``
combines the per-setting comparisons. Outputs land in
``06_switch_mechanism/_m/switch_feature_sensitivity/refit/<cohort>/<region>/``.
"""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from isograph.io.artifacts import load_dataset_bundle
from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.run_models import (
    BRAINSEQ_DISCOVERY_COVARIATES,
    CANONICAL_LEIDEN_RESOLUTION,
    GTEX_DISCOVERY_COVARIATES,
    _PROMOTED_VAE,
    _eigengenes_to_sample_table,
    linear_age_association,
)
from isograph_benchmark.real_data.switch_feature_sensitivity import (
    COHORTS,
    PUBLISHED,
    PUBLISHED_BY_COHORT,
    _artifact_dir,
    _drop_minor_isoforms,
    apply_transcript_filter,
    _label,
    _root_out,
    _settings,
)

_DISCOVERY_COVARIATES = {"brainseq": BRAINSEQ_DISCOVERY_COVARIATES,
                         "gtex": GTEX_DISCOVERY_COVARIATES}


def setting_grid(cohort: str) -> list[dict]:
    """Every distinct setting across the three perturbation axes, published first, once."""
    published = PUBLISHED_BY_COHORT[cohort]
    grid = [{"axis": "published", "label": "published", "setting": dict(published)}]
    for axis in ("pseudocount", "expression", "minor_isoform"):
        for s in _settings(axis, published):
            if all(s[k] == published[k] for k in published):
                continue
            grid.append({"axis": axis, "label": _label(s, axis), "setting": s})
    return grid


def slug(entry: dict) -> str:
    return re.sub(r"[^A-Za-z0-9.=_-]+", "_", f"{entry['axis']}__{entry['label']}")


def _out_dir(cohort: str, region: str) -> Path:
    return ensure_dir(_root_out() / "refit" / cohort / region)


def counts_for_setting(counts: np.ndarray, table: pd.DataFrame, setting: dict):
    """Transcript counts under one setting, exactly as ``build_features`` prepares them."""
    tc, tt = apply_transcript_filter(np.asarray(counts), table, setting)
    tc, tt = _drop_minor_isoforms(tc, tt, setting["min_usage"])
    if setting["pseudocount"] != PUBLISHED["pseudocount"]:
        tc = tc.astype(np.float64) + (setting["pseudocount"] - PUBLISHED["pseudocount"])
    return tc, tt


def production_config(cohort: str):
    """The production VaeModelConfig (run_models.run_{brainseq,gtex}_region), canonical resolution."""
    from isograph.workflow.config import VaeModelConfig

    return VaeModelConfig(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        residualize_covariates=_DISCOVERY_COVARIATES[cohort],
        min_module_size=20, trait_columns=[COHORTS[cohort]["age_col"]], random_state=13,
        allow_abundance_abundance=True,
        alpha_switch=0.5,
        alpha_abundance_grid=[0.70, 0.75, 0.80, 0.85, 0.90, 0.95],
        leiden_resolution=CANONICAL_LEIDEN_RESOLUTION,
        **_PROMOTED_VAE,
    )


def compare_partitions(published: pd.DataFrame, refit: pd.DataFrame,
                       published_age: pd.DataFrame, refit_age: pd.DataFrame,
                       fdr: float = 0.05) -> dict:
    """Partition agreement and survival of the published module--age associations.

    Agreement is over genes assigned in both partitions. For each published FDR-significant
    module, the best-Jaccard refit module is its counterpart; the association is *retained* when
    that counterpart is itself FDR-significant with the same sign.
    """
    p = published[["gene_id", "module_id"]].astype(str)
    r = refit[["gene_id", "module_id"]].astype(str)
    m = p.merge(r, on="gene_id", suffixes=("_pub", "_ref"))
    out = {
        "n_modules_published": int(p["module_id"].nunique()),
        "n_modules_refit": int(r["module_id"].nunique()),
        "n_genes_published": int(len(p)),
        "n_genes_refit": int(len(r)),
        "n_genes_both": int(len(m)),
        "ari": float(adjusted_rand_score(m["module_id_pub"], m["module_id_ref"])) if len(m) else np.nan,
        "nmi": float(normalized_mutual_info_score(m["module_id_pub"], m["module_id_ref"])) if len(m) else np.nan,
    }
    pub_sets = p.groupby("module_id")["gene_id"].apply(set)
    ref_sets = r.groupby("module_id")["gene_id"].apply(set)
    pa = published_age.set_index(published_age["module_id"].astype(str))
    ra = refit_age.set_index(refit_age["module_id"].astype(str)) if not refit_age.empty else pd.DataFrame()
    sig = pa.index[pa["fdr"] < fdr]
    best_j, retained = [], 0
    for mid in sig:
        genes = pub_sets.get(mid, set())
        if not genes or ref_sets.empty:
            best_j.append(0.0)
            continue
        jac = ref_sets.apply(lambda g: len(genes & g) / len(genes | g))
        best = jac.idxmax()
        best_j.append(float(jac.max()))
        if best in ra.index and ra.loc[best, "fdr"] < fdr and \
                np.sign(ra.loc[best, "effect"]) == np.sign(pa.loc[mid, "effect"]):
            retained += 1
    out.update({
        "n_age_sig_published": int(len(sig)),
        "n_age_sig_refit": int((ra["fdr"] < fdr).sum()) if not ra.empty else 0,
        "n_age_sig_retained": int(retained),
        "median_best_jaccard_age_sig": float(np.median(best_j)) if best_j else np.nan,
    })
    return out


def fit_one(args) -> None:
    grid = setting_grid(args.cohort)
    if not 1 <= args.index <= len(grid):
        raise SystemExit(f"--index must be 1..{len(grid)} for {args.cohort}")
    entry = grid[args.index - 1]
    out = ensure_dir(_out_dir(args.cohort, args.region) / slug(entry))
    print(f"[refit {args.index}/{len(grid)}] {args.cohort}/{args.region} {entry['axis']}: "
          f"{entry['label']} -> {out}", flush=True)

    from isograph.models.vae import VaeNetworkModel

    bundle = load_dataset_bundle(rel(*COHORTS[args.cohort]["bundle"], args.region))
    sample_table = bundle.sample_table
    tc, tt = counts_for_setting(
        bundle.matrices["transcript_counts"], bundle.feature_tables["transcript"], entry["setting"]
    )
    del bundle
    t0 = time.time()
    artifacts = VaeNetworkModel(production_config(args.cohort)).fit(
        transcript_counts=tc, transcript_table=tt, sample_table=sample_table,
    )
    del tc, tt
    elapsed = time.time() - t0

    modules = artifacts.module_table
    eig = _eigengenes_to_sample_table(artifacts)
    age = (linear_age_association(eig, sample_table, COHORTS[args.cohort]["age_col"])
           if eig is not None else pd.DataFrame(columns=["module_id", "effect", "fdr"]))
    modules.to_parquet(out / "modules.parquet", index=False)
    age.to_parquet(out / "age_linear.parquet", index=False)

    art = _artifact_dir(args.cohort, args.region)
    comparison = compare_partitions(
        pd.read_parquet(art / "modules.parquet"), modules,
        pd.read_parquet(art / "age_linear.parquet"), age, fdr=args.fdr,
    )
    record = {
        "cohort": args.cohort, "region": args.region, "index": args.index,
        "axis": entry["axis"], "setting": entry["label"], **entry["setting"],
        "fit_seconds": round(elapsed, 1),
        "selected_alpha_abundance": (artifacts.calibration or {}).get("selected_alpha_abundance"),
        **comparison,
    }
    (out / "comparison.json").write_text(json.dumps(record, indent=2, default=float))
    print(json.dumps(record, indent=2, default=float), flush=True)


def aggregate(args) -> None:
    root = _out_dir(args.cohort, args.region)
    rows = [json.loads(p.read_text()) for p in sorted(root.glob("*/comparison.json"))]
    if not rows:
        raise SystemExit(f"no refit comparisons under {root}")
    df = pd.DataFrame(rows).sort_values("index")
    df.to_parquet(root / "refit_summary.parquet", index=False)
    floor = df[df["axis"] == "published"]
    cols = ["axis", "setting", "n_modules_refit", "n_genes_both", "ari", "nmi",
            "n_age_sig_published", "n_age_sig_refit", "n_age_sig_retained",
            "median_best_jaccard_age_sig", "selected_alpha_abundance"]

    def md(frame):
        def cell(v):
            if isinstance(v, float):
                return "" if pd.isna(v) else f"{v:.3g}"
            return "" if v is None else str(v)
        return "\n".join(["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols),
                          *["| " + " | ".join(cell(v) for v in r) + " |"
                            for r in frame[cols].itertuples(index=False, name=None)]])

    lines = [
        "# Full refit per preprocessing setting",
        "",
        f"`switch_feature_refit.py` on **{args.cohort}/{args.region}**: {len(df)} fits of the "
        "production IsoGraph configuration, one per distinct setting of the pseudocount, "
        "expression-filter and minor-isoform axes, each compared with the committed fit.",
        "",
        "## The noise floor",
        "",
        "The `published` row refits the published setting itself. Its agreement with the "
        "committed partition is the ceiling for every other row; a perturbed setting that "
        "matches it has cost nothing beyond refit noise.",
        "",
    ]
    if not floor.empty:
        f = floor.iloc[0]
        lines += [f"Published-setting refit: ARI **{f['ari']:.3f}**, NMI **{f['nmi']:.3f}**, "
                  f"{int(f['n_age_sig_retained'])}/{int(f['n_age_sig_published'])} published "
                  "age-significant modules retained.", ""]
    lines += ["## Every setting", "", md(df), "",
              "`n_age_sig_retained`: published FDR-significant module–age associations whose "
              "best-Jaccard refit counterpart is FDR-significant with the same sign. Agreement "
              "metrics are over genes assigned in both partitions.", ""]
    (root / "REFIT_SENSITIVITY.md").write_text("\n".join(lines))
    print(df[cols].to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cohort", choices=list(COHORTS), default="brainseq")
    ap.add_argument("--region", default="caudate")
    ap.add_argument("--index", type=int, help="1-based setting index (one fit per process)")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--list", action="store_true", help="print the setting grid and exit")
    args = ap.parse_args()
    if args.list:
        for i, e in enumerate(setting_grid(args.cohort), 1):
            print(i, slug(e))
    elif args.aggregate:
        aggregate(args)
    elif args.index is not None:
        fit_one(args)
    else:
        ap.error("pass --index N, --aggregate or --list")


if __name__ == "__main__":
    main()
