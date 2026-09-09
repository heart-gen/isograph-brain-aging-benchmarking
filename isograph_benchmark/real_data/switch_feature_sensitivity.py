"""Do the conclusions survive the arbitrary choices in building switch features?

A reviewer asked for sensitivity analyses over five preprocessing choices that are set by
convention rather than derived: the compositional pseudocount, the transcript-expression
filter, the minor-isoform threshold, transcript number / identifiability, and the
quantification pipeline. None of them had a harness. This is that harness.

**What is varied and what is held fixed.** The switch channel is PC1 of each gene's
within-gene CLR composition, and it is fully determined by the transcript count matrix plus
those preprocessing choices -- so it can be rebuilt from persisted inputs without refitting
the VAE. This CLI rebuilds it under each setting, then holds the **published module
partition fixed** and recomputes the module eigengenes and their age association. Varying
preprocessing while re-clustering would confound "the features moved" with "the clustering
moved"; holding the partition fixed isolates the first, which is what the question asks.

That scoping is a real limitation and is reported as one: a full refit under each setting
would additionally let the network change, and is the expensive follow-on. What is measured
here is whether the *representation and its trait signal* are stable, not whether an
independently refit network would recover the same modules.

**The baseline is verified, not assumed.** Before any perturbation, the harness rebuilds the
switch channel at the published settings and requires it to reproduce
``feature_scores.parquet`` exactly (max |diff| = 0). If that gate fails the run aborts,
because every downstream comparison would otherwise be against a lookalike.

Axes
----
``pseudocount``      the constant added to transcript counts before the composition is taken
``expression``       the (min_count, min_fraction) transcript filter
``minor_isoform``    drop transcripts below a mean within-gene usage before the CLR
``identifiability``  no perturbation -- stratifies genes by transcript number and asks
                     whether the age signal concentrates where isoforms are hardest to
                     resolve, which is what an identifiability artefact would look like
``quantification``   Salmon (BrainSEQ) versus RSEM (GTEx) on shared genes and a shared
                     region; the pipelines cannot be swapped within a cohort, so this is a
                     cross-cohort concordance, and it is reported as such

Outputs land in ``06_switch_mechanism/_m/switch_feature_sensitivity/``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from isograph.features.channels import gene_feature_channels
from isograph.io.artifacts import load_dataset_bundle
from isograph_benchmark.paths import (
    OUTPUT_DIRS,
    ensure_dir,
    region_store,
    rel,
    stage_out,
)
from isograph_benchmark.real_data.run_models import (
    _filter_expressed_transcripts,
    linear_age_association,
)

# Published preprocessing (run_models.run_brainseq_region / run_gtex_region).
PUBLISHED = {"pseudocount": 0.5, "min_count": 10.0, "min_fraction": 0.70, "min_usage": 0.0}

COHORTS = {
    "brainseq": {
        "bundle": ("inputs", "bundles", "brainseq_v1"),
        "artifacts": (*OUTPUT_DIRS["modules"], "brainseq"),
        "age_col": "Age",
        "quantifier": "Salmon",
    },
    "gtex": {
        "bundle": ("inputs", "bundles", "gtex_v11_brain"),
        "artifacts": (*OUTPUT_DIRS["modules"], "gtex"),
        "age_col": "AGE",
        "quantifier": "RSEM",
    },
}

_META = {"feature_id", "gene_id", "feature_type", "n_transcripts"}


def _out_dir() -> Path:
    return ensure_dir(stage_out("mechanism", "switch_feature_sensitivity"))


def _artifact_dir(cohort: str, region: str) -> Path:
    return region_store(cohort, region, "isograph_vae")


# --------------------------------------------------------------------------- #
# Feature construction under a given setting
# --------------------------------------------------------------------------- #
def _drop_minor_isoforms(
    counts: np.ndarray, table: pd.DataFrame, min_usage: float
) -> tuple[np.ndarray, pd.DataFrame]:
    """Drop transcripts whose mean within-gene usage falls below ``min_usage``.

    A gene reduced to a single transcript loses its switch channel entirely, which is the
    honest consequence of a strict threshold and is counted rather than hidden.
    """
    if min_usage <= 0:
        return counts, table
    keep = np.ones(len(table), dtype=bool)
    totals = np.zeros(counts.shape[1])
    for _, idx in table.groupby("gene_id", sort=False).indices.items():
        sub = counts[idx]
        totals = sub.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            usage = np.where(totals > 0, sub / np.where(totals > 0, totals, np.nan), np.nan)
        keep[idx] = np.nanmean(usage, axis=1) >= min_usage
    return counts[keep], table.loc[keep].reset_index(drop=True)


def build_features(
    counts: np.ndarray,
    table: pd.DataFrame,
    pseudocount: float,
    min_count: float,
    min_fraction: float,
    min_usage: float,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Rebuild the multiplex feature matrix under one preprocessing setting.

    ``gene_feature_channels`` is called with no ``switch_design`` because the production fit
    leaves the composition unresidualized (``residualize_composition`` is False) and
    persists the raw channels; matching that is what makes the baseline reproduce exactly.
    """
    tc, tt = _filter_expressed_transcripts(
        counts, table, min_count=min_count, min_fraction=min_fraction
    )
    tc, tt = _drop_minor_isoforms(tc, tt, min_usage)
    if pseudocount != PUBLISHED["pseudocount"]:
        # transcript_usage adds PUBLISHED["pseudocount"] internally, so shifting by the
        # difference yields counts + pseudocount exactly. No clipping: the shifted value is
        # >= pseudocount - 0.5 > -0.5, so the internal +0.5 keeps every entry positive, and
        # clipping here would silently restore 0.5 for zero counts whenever pseudocount<0.5.
        tc = tc.astype(np.float64) + (pseudocount - PUBLISHED["pseudocount"])
    return gene_feature_channels(tc, tt)


def _to_frame(matrix: np.ndarray, info: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    out = pd.DataFrame(matrix, columns=samples)
    out.insert(0, "feature_type", info["feature_type"].to_numpy())
    out.insert(0, "gene_id", info["gene_id"].astype(str).to_numpy())
    return out


# --------------------------------------------------------------------------- #
# Downstream: fixed modules -> eigengenes -> age association
# --------------------------------------------------------------------------- #
def eigengenes_from_features(
    features: pd.DataFrame, modules: pd.DataFrame, samples: list[str]
) -> pd.DataFrame:
    gene_to_module = dict(
        zip(modules["gene_id"].astype(str), modules["module_id"].astype(str))
    )
    f = features[features["gene_id"].isin(gene_to_module)].copy()
    f["module_id"] = f["gene_id"].map(gene_to_module)
    eig = f.groupby("module_id")[samples].mean().sort_index()
    out = eig.T.reset_index().rename(columns={"index": "sample_id"})
    out["sample_id"] = out["sample_id"].astype(str)
    return out


def _compare_to_published(
    assoc: pd.DataFrame, published: pd.DataFrame, fdr: float
) -> dict:
    m = published.merge(assoc, on="module_id", suffixes=("_pub", "_new"))
    if m.empty:
        return {"n_modules": 0}
    sig_pub = m["fdr_pub"] < fdr
    return {
        "n_modules": int(len(m)),
        "effect_pearson_vs_published": float(
            stats.pearsonr(m["effect_pub"], m["effect_new"]).statistic
        ),
        "effect_spearman_vs_published": float(
            stats.spearmanr(m["effect_pub"], m["effect_new"]).statistic
        ),
        "median_abs_effect_shift": float((m["effect_pub"] - m["effect_new"]).abs().median()),
        "n_fdr_sig_published": int(sig_pub.sum()),
        "n_fdr_sig": int((m["fdr_new"] < fdr).sum()),
        "n_published_sig_retained": int((sig_pub & (m["fdr_new"] < fdr)).sum()),
        "n_sign_flips_among_published_sig": int(
            (sig_pub & (np.sign(m["effect_pub"]) != np.sign(m["effect_new"]))).sum()
        ),
    }


# --------------------------------------------------------------------------- #
# Axes
# --------------------------------------------------------------------------- #
def _settings(axis: str) -> list[dict]:
    base = dict(PUBLISHED)
    if axis == "pseudocount":
        return [dict(base, pseudocount=v) for v in (0.1, 0.25, 0.5, 1.0, 2.0)]
    if axis == "expression":
        return [
            dict(base, min_count=c, min_fraction=f)
            for c, f in ((5.0, 0.50), (10.0, 0.50), (10.0, 0.70), (20.0, 0.70), (10.0, 0.90))
        ]
    if axis == "minor_isoform":
        return [dict(base, min_usage=v) for v in (0.0, 0.01, 0.05, 0.10)]
    raise ValueError(axis)


def _label(setting: dict, axis: str) -> str:
    if axis == "pseudocount":
        return f"pseudocount={setting['pseudocount']:g}"
    if axis == "expression":
        return f"count>{setting['min_count']:g},frac>={setting['min_fraction']:g}"
    return f"min_usage={setting['min_usage']:g}"


def run_perturbation_axes(args) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    art = _artifact_dir(args.cohort, args.region)
    bundle = load_dataset_bundle(rel(*COHORTS[args.cohort]["bundle"], args.region))
    sample_table = bundle.sample_table
    counts = np.asarray(bundle.matrices["transcript_counts"])
    table = bundle.feature_tables["transcript"]

    published_fs = pd.read_parquet(art / "feature_scores.parquet")
    samples = [c for c in published_fs.columns if c not in _META]
    modules = pd.read_parquet(art / "modules.parquet")
    published_age = pd.read_parquet(art / "age_linear.parquet")

    # ---- gate: the baseline must reproduce the published switch channel exactly ----
    m0, i0 = build_features(counts, table, **PUBLISHED)
    base = _to_frame(m0, i0, samples)
    pub_sw = published_fs[published_fs["feature_type"] == "switch"].set_index("gene_id")
    new_sw = base[base["feature_type"] == "switch"].set_index("gene_id")
    shared = pub_sw.index.intersection(new_sw.index)
    worst = float(
        np.nanmax(
            np.abs(
                pub_sw.loc[shared, samples].to_numpy(float)
                - new_sw.loc[shared, samples].to_numpy(float)
            )
        )
    )
    if worst > args.tol or len(shared) != len(pub_sw):
        raise SystemExit(
            f"baseline switch channel does not reproduce feature_scores.parquet "
            f"(max |diff| = {worst:.3g} over {len(shared)}/{len(pub_sw)} genes); every "
            "sensitivity below would be measured against the wrong baseline."
        )
    print(f"[verify] baseline reproduces the published switch channel exactly "
          f"(max |diff| = {worst:g}, {len(shared)} genes)")

    base_sw = new_sw[samples]
    rows, per_gene_rows = [], []
    for axis in ("pseudocount", "expression", "minor_isoform"):
        for setting in _settings(axis):
            label = _label(setting, axis)
            is_base = all(setting[k] == PUBLISHED[k] for k in PUBLISHED)
            m, i = (m0, i0) if is_base else build_features(counts, table, **setting)
            feats = _to_frame(m, i, samples)
            sw = feats[feats["feature_type"] == "switch"].set_index("gene_id")

            common = base_sw.index.intersection(sw.index)
            A = base_sw.loc[common].to_numpy(float)
            B = sw.loc[common, samples].to_numpy(float)
            Az = A - A.mean(1, keepdims=True)
            Bz = B - B.mean(1, keepdims=True)
            denom = np.sqrt((Az**2).sum(1) * (Bz**2).sum(1))
            with np.errstate(invalid="ignore", divide="ignore"):
                r = np.where(denom > 0, (Az * Bz).sum(1) / np.where(denom > 0, denom, np.nan), np.nan)

            eig = eigengenes_from_features(feats, modules, samples)
            assoc = linear_age_association(eig, sample_table, COHORTS[args.cohort]["age_col"])
            rec = {
                "axis": axis,
                "setting": label,
                "is_published": is_base,
                **{k: setting[k] for k in PUBLISHED},
                "n_switch_genes": int(len(sw)),
                "n_switch_genes_published": int(len(base_sw)),
                "n_switch_genes_lost": int(len(base_sw) - len(common)),
                "median_abs_feature_r_vs_published": float(np.nanmedian(np.abs(r))),
                "frac_features_r_above_0.95": float(np.nanmean(np.abs(r) > 0.95)),
                "frac_features_sign_flipped": float(np.nanmean(r < 0)),
                **_compare_to_published(assoc, published_age, args.fdr),
            }
            rows.append(rec)
            per_gene_rows.append(
                pd.DataFrame({"axis": axis, "setting": label, "gene_id": common, "feature_r": r})
            )
            print(f"  [{axis}] {label}: {rec['n_switch_genes']} switch genes, "
                  f"median |r| {rec['median_abs_feature_r_vs_published']:.3f}, "
                  f"age r vs published {rec.get('effect_pearson_vs_published', float('nan')):.3f}, "
                  f"sig {rec.get('n_fdr_sig')}/{rec.get('n_fdr_sig_published')}", flush=True)

    meta = {
        "cohort": args.cohort,
        "region": args.region,
        "quantifier": COHORTS[args.cohort]["quantifier"],
        "artifacts": str(art),
        "published_settings": PUBLISHED,
        "baseline_max_abs_diff_vs_feature_scores": worst,
        "n_samples": int(len(sample_table)),
        "n_modules": int(modules["module_id"].nunique()),
        "fdr": args.fdr,
        "scope": (
            "module partition held fixed at the published one; preprocessing is varied and "
            "the features, eigengenes and age association are recomputed. A full refit per "
            "setting would additionally let the network change and is not done here."
        ),
    }
    return pd.DataFrame(rows), pd.concat(per_gene_rows, ignore_index=True), meta


# --------------------------------------------------------------------------- #
# Identifiability (stratification, not perturbation)
# --------------------------------------------------------------------------- #
def run_identifiability(args) -> pd.DataFrame:
    """Does the age signal concentrate in genes with many transcripts?

    A gene with many annotated isoforms is harder to quantify, so if the switch signal were
    a quantification artefact it should grow with transcript number. Per-gene switch--age
    correlation is computed on the published features and summarised by transcript-count
    stratum, with the module-membership rate alongside.
    """
    art = _artifact_dir(args.cohort, args.region)
    fs = pd.read_parquet(art / "feature_scores.parquet")
    samples = [c for c in fs.columns if c not in _META]
    modules = pd.read_parquet(art / "modules.parquet")
    bundle = load_dataset_bundle(rel(*COHORTS[args.cohort]["bundle"], args.region))
    st = bundle.sample_table.set_index(bundle.sample_table["sample_id"].astype(str))
    age = pd.to_numeric(st.loc[samples, COHORTS[args.cohort]["age_col"]], errors="coerce").to_numpy(float)

    sw = fs[fs["feature_type"] == "switch"].copy()
    Y = sw[samples].to_numpy(float)
    ok = np.isfinite(age)
    a = age[ok] - age[ok].mean()
    Yz = Y[:, ok] - Y[:, ok].mean(axis=1, keepdims=True)
    denom = np.sqrt((Yz**2).sum(1) * (a**2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(denom > 0, (Yz * a).sum(1) / np.where(denom > 0, denom, np.nan), np.nan)

    sw = sw.assign(age_r=r, in_module=sw["gene_id"].astype(str).isin(
        set(modules["gene_id"].astype(str))))
    bins = [1, 2, 3, 5, 10, 20, 10_000]
    sw["n_tx_stratum"] = pd.cut(sw["n_transcripts"], bins=bins, right=False)
    out = (
        sw.groupby("n_tx_stratum", observed=True)
        .agg(
            n_genes=("gene_id", "size"),
            median_abs_age_r=("age_r", lambda s: float(np.nanmedian(np.abs(s)))),
            frac_abs_age_r_above_0_2=("age_r", lambda s: float(np.nanmean(np.abs(s) > 0.2))),
            module_membership_rate=("in_module", "mean"),
        )
        .reset_index()
    )
    out["n_tx_stratum"] = out["n_tx_stratum"].astype(str)
    out["cohort"] = args.cohort
    out["region"] = args.region
    return out


# --------------------------------------------------------------------------- #
# Quantification pipeline (cross-cohort concordance)
# --------------------------------------------------------------------------- #
def run_quantification(args) -> pd.DataFrame:
    """Salmon (BrainSEQ) versus RSEM (GTEx) on a shared region and shared genes.

    The pipelines cannot be swapped within a cohort -- neither cohort releases a second
    quantification -- so this is a concordance between cohorts that differ in quantifier,
    and it therefore confounds quantifier with cohort. It is reported as an upper bound on
    how much of the cross-cohort attenuation the quantifier could explain, not as an
    isolated quantifier effect.
    """
    pairs = [("brainseq", "caudate", "gtex", "caudate_basal_ganglia"),
             ("brainseq", "hippocampus", "gtex", "hippocampus")]
    rows = []
    for c1, r1, c2, r2 in pairs:
        try:
            f1 = pd.read_parquet(_artifact_dir(c1, r1) / "feature_scores.parquet")
            f2 = pd.read_parquet(_artifact_dir(c2, r2) / "feature_scores.parquet")
        except FileNotFoundError:
            continue
        s1 = [c for c in f1.columns if c not in _META]
        s2 = [c for c in f2.columns if c not in _META]
        b1 = load_dataset_bundle(rel(*COHORTS[c1]["bundle"], r1)).sample_table
        b2 = load_dataset_bundle(rel(*COHORTS[c2]["bundle"], r2)).sample_table

        def gene_age_r(f, samples, bundle_st, age_col):
            st = bundle_st.set_index(bundle_st["sample_id"].astype(str))
            age = pd.to_numeric(st.loc[samples, age_col], errors="coerce").to_numpy(float)
            sw = f[f["feature_type"] == "switch"]
            Y = sw[samples].to_numpy(float)
            ok = np.isfinite(age)
            a = age[ok] - age[ok].mean()
            Yz = Y[:, ok] - Y[:, ok].mean(axis=1, keepdims=True)
            den = np.sqrt((Yz**2).sum(1) * (a**2).sum())
            with np.errstate(invalid="ignore", divide="ignore"):
                r = np.where(den > 0, (Yz * a).sum(1) / np.where(den > 0, den, np.nan), np.nan)
            return pd.Series(r, index=sw["gene_id"].astype(str).str.split(".").str[0])

        g1 = gene_age_r(f1, s1, b1, COHORTS[c1]["age_col"])
        g2 = gene_age_r(f2, s2, b2, COHORTS[c2]["age_col"])
        common = g1.index.intersection(g2.index)
        v1, v2 = g1.loc[common].to_numpy(float), g2.loc[common].to_numpy(float)
        ok = np.isfinite(v1) & np.isfinite(v2)
        rows.append(
            {
                "cohort_1": c1, "region_1": r1, "quantifier_1": COHORTS[c1]["quantifier"],
                "cohort_2": c2, "region_2": r2, "quantifier_2": COHORTS[c2]["quantifier"],
                "n_shared_genes": int(ok.sum()),
                "pearson_gene_age_effect": float(stats.pearsonr(v1[ok], v2[ok]).statistic),
                "spearman_gene_age_effect": float(stats.spearmanr(v1[ok], v2[ok]).statistic),
                "sign_concordance": float(np.mean(np.sign(v1[ok]) == np.sign(v2[ok]))),
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
def run(args) -> None:
    out_dir = _out_dir()
    summary, per_gene, meta = run_perturbation_axes(args)
    ident = run_identifiability(args)
    quant = run_quantification(args)

    summary.to_parquet(out_dir / "sensitivity_summary.parquet", index=False)
    per_gene.to_parquet(out_dir / "sensitivity_per_gene.parquet", index=False)
    ident.to_parquet(out_dir / "identifiability_strata.parquet", index=False)
    if not quant.empty:
        quant.to_parquet(out_dir / "quantification_concordance.parquet", index=False)
    (out_dir / "sensitivity.json").write_text(json.dumps(meta, indent=2, default=str))
    _write_report(out_dir, summary, ident, quant, meta)
    print(summary.to_string(index=False))


def _md(frame: pd.DataFrame, floats: int = 4) -> str:
    def cell(v):
        if isinstance(v, float):
            return "" if pd.isna(v) else f"{v:.{floats}g}"
        return "" if v is None else str(v)

    cols = list(frame.columns)
    return "\n".join(
        [
            "| " + " | ".join(cols) + " |",
            "|" + "|".join("---" for _ in cols) + "|",
            *["| " + " | ".join(cell(v) for v in row) + " |"
              for row in frame.itertuples(index=False, name=None)],
        ]
    )


def _write_report(out_dir, summary, ident, quant, meta) -> None:
    keep = [
        "axis", "setting", "is_published", "n_switch_genes", "n_switch_genes_lost",
        "median_abs_feature_r_vs_published", "frac_features_sign_flipped",
        "effect_pearson_vs_published", "n_fdr_sig", "n_fdr_sig_published",
        "n_published_sig_retained", "n_sign_flips_among_published_sig",
    ]
    lines = [
        "# Preprocessing sensitivity of the switch representation",
        "",
        f"`switch_feature_sensitivity.py` on **{meta['cohort']}/{meta['region']}** "
        f"({meta['quantifier']}, n={meta['n_samples']}, {meta['n_modules']} modules).",
        "",
        "## Scope, stated up front",
        "",
        meta["scope"],
        "",
        f"Gate passed: the baseline rebuild reproduces the published switch channel with "
        f"max |diff| = {meta['baseline_max_abs_diff_vs_feature_scores']:g}, so every "
        "comparison below is against the published quantity and not a lookalike.",
        "",
        "## 1-3. Pseudocount, expression filter, minor-isoform threshold",
        "",
        _md(summary[keep]),
        "",
        "`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt "
        "switch coordinate with the published one; `n_published_sig_retained` is how many "
        "of the published FDR-significant module-age associations survive. A sign flip in "
        "a switch coordinate is not itself a problem -- PC1's sign is arbitrary and "
        "sign-stabilised -- but a flip among *published-significant modules* would change "
        "the direction of a reported effect, so it is counted separately.",
        "",
        "## 4. Identifiability (transcript number)",
        "",
        "If the switch signal were a quantification artefact it should grow with the number "
        "of annotated isoforms, which is what makes a gene hard to quantify. It is reported "
        "by stratum rather than adjusted away:",
        "",
        _md(ident),
        "",
    ]
    if not quant.empty:
        lines += [
            "## 5. Quantification pipeline",
            "",
            "Neither cohort releases a second quantification, so the pipelines cannot be "
            "swapped within a cohort. This compares matched regions across cohorts that "
            "differ in quantifier, which **confounds quantifier with cohort** and is "
            "therefore an upper bound on the quantifier's contribution, not an isolated "
            "estimate of it.",
            "",
            _md(quant),
            "",
        ]
    (out_dir / "SWITCH_FEATURE_SENSITIVITY.md").write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cohort", choices=list(COHORTS), default="brainseq")
    ap.add_argument("--region", default="caudate")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--tol", type=float, default=0.0,
                    help="max |diff| allowed between the rebuilt and published switch "
                         "channel; 0 because the rebuild is exact by construction")
    args = ap.parse_args()
    run(args)


if __name__ == "__main__":
    main()
