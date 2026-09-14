"""Cross-cohort eigengene projection: does a module's frozen signature exist in the other cohort?

Stage 03's cross-cohort arm matches modules by gene-set overlap and compares each cohort's own
eigengene. Overlap is not the same test as projection: two cohorts can share genes while
weighting them differently, and a best-Jaccard match is granularity-confounded. This CLI asks
the WGCNA module-preservation question directly. For each trusted module in the **source**
cohort it freezes the module's eigengene weights, applies them unchanged to the **target**
cohort's expression of the same features, and tests two things:

1. **Preservation** -- do the target features still covary along the frozen axis? The statistic
   is the mean signed kME: the correlation of each target feature with the projected eigengene,
   signed by that feature's source weight. It is compared with a null that applies the same
   weight vector to size- and type-matched random target features.
2. **Aging** -- does the projected eigengene carry the source module's age association, with the
   same sign? Source and target effects use the same frozen-eigengene definition.

**Switch axes are oriented before projection.** A gene's switch coordinate is PC1 of its
within-gene composition, and PC1's sign is fixed per cohort by a pivot that need not agree
across cohorts. Applying a source weight to a target switch row with the opposite orientation
would silently invert that gene's contribution. Both cohorts share GENCODE transcript ids, so
each gene's switch loadings are compared over the shared transcripts: the cosine sets the sign,
and a gene whose |cosine| < 0.5 (a different axis in the other cohort) or with fewer than two
shared transcripts is dropped from the switch rows rather than guessed. Abundance rows need no
alignment.

Both directions (BrainSEQ -> GTEx and GTEx -> BrainSEQ) over the three matched region pairs,
for IsoGraph (switch + abundance rows) and classical WGCNA (abundance rows). No covariates enter
the age test, matching the published linear arm.

  python -m isograph_benchmark.real_data.eigengene_projection run --pair caudate --method isograph
  python -m isograph_benchmark.real_data.eigengene_projection aggregate
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.paths import ensure_dir, region_store, rel, stage_out
from isograph_benchmark.real_data.module_trust import (
    REGION_PAIRS,
    _load_production_modules,
    _trusted_set,
)

BUNDLES = {"brainseq": ("inputs", "bundles", "brainseq_v1"),
           "gtex": ("inputs", "bundles", "gtex_v11_brain")}
AGE_COL = {"brainseq": "Age", "gtex": "AGE"}
FEATURE_TYPES = {"isograph": ("abundance", "switch"), "wgcna": ("abundance",)}
_META = {"feature_id", "gene_id", "feature_type", "n_transcripts"}
MIN_FEATURES = 10
MIN_SHARED_TRANSCRIPTS = 2
MIN_ABS_COSINE = 0.5


def _out_dir() -> Path:
    return ensure_dir(stage_out("trust.stability", "eigengene_projection"))


def _bare(s: pd.Series) -> pd.Series:
    return s.astype(str).str.split(".").str[0]


# --------------------------------------------------------------------------- #
# Switch-axis orientation across cohorts
# --------------------------------------------------------------------------- #
def switch_loadings(cohort: str, region: str) -> pd.DataFrame:
    """Per-transcript switch loadings on the production feature construction."""
    from isograph.features.switch import gene_switch_loadings
    from isograph.io.artifacts import load_dataset_bundle
    from isograph_benchmark.real_data.run_models import _filter_expressed_transcripts

    bundle = load_dataset_bundle(rel(*BUNDLES[cohort], region))
    tc = np.asarray(bundle.matrices["transcript_counts"])
    tt = bundle.feature_tables["transcript"]
    tc, tt = _filter_expressed_transcripts(tc, tt)  # both production fits filter
    load = gene_switch_loadings(tc, tt)
    load["gene"] = _bare(load["gene_id"])
    load["transcript"] = _bare(load["transcript_id"])
    return load[["gene", "transcript", "loading"]]


def align_switch_axes(source: pd.DataFrame, target: pd.DataFrame,
                      min_shared: int = MIN_SHARED_TRANSCRIPTS,
                      min_abs_cos: float = MIN_ABS_COSINE) -> pd.DataFrame:
    """Per shared gene: cosine of switch loadings over shared transcripts, and the sign that
    orients the TARGET axis like the SOURCE one. `usable` is False when the axis cannot be
    oriented (too few shared transcripts, or the two cohorts' PC1 are different axes)."""
    m = source.merge(target, on=["gene", "transcript"], suffixes=("_s", "_t"))
    rows = []
    for gene, g in m.groupby("gene", sort=False):
        a, b = g["loading_s"].to_numpy(float), g["loading_t"].to_numpy(float)
        den = np.linalg.norm(a) * np.linalg.norm(b)
        cos = float(a @ b / den) if den > 0 else np.nan
        rows.append((gene, len(g), cos))
    t = pd.DataFrame(rows, columns=["gene", "n_shared_transcripts", "cosine"])
    t["sign"] = np.where(t["cosine"] < 0, -1.0, 1.0)
    t["usable"] = (t["n_shared_transcripts"] >= min_shared) & (t["cosine"].abs() >= min_abs_cos)
    return t


# --------------------------------------------------------------------------- #
# Projection
# --------------------------------------------------------------------------- #
def _zrows(x: np.ndarray) -> np.ndarray:
    c = x - np.nanmean(x, axis=1, keepdims=True)
    sd = np.nanstd(c, axis=1, keepdims=True)
    return np.nan_to_num(c / np.where(sd > 1e-12, sd, 1.0))


def frozen_weights(xs: np.ndarray) -> np.ndarray:
    """PC1 loadings of the (row-standardized) source feature matrix, oriented so the
    eigengene correlates positively with the mean of its features."""
    _, _, vt = np.linalg.svd(xs.T - xs.T.mean(axis=0), full_matrices=False)
    w = vt[0]
    score = xs.T @ w
    if np.corrcoef(score, xs.mean(axis=0))[0, 1] < 0:
        w = -w
    return w / np.linalg.norm(w)


def project(w: np.ndarray, xt: np.ndarray) -> np.ndarray:
    e = xt.T @ w
    sd = e.std()
    return (e - e.mean()) / (sd if sd > 1e-12 else 1.0)


def signed_kme(w: np.ndarray, xt: np.ndarray, e: np.ndarray) -> float:
    """Mean correlation of each target feature with the projected eigengene, signed by the
    feature's frozen weight. Features are row-standardized, so the correlation is the mean
    product divided by the sample count."""
    r = (xt @ e) / xt.shape[1]
    return float(np.mean(np.sign(w) * r))


def _age_r(e: np.ndarray, age: np.ndarray) -> tuple[float, float]:
    ok = np.isfinite(age)
    if ok.sum() < 10:
        return np.nan, np.nan
    r, p = stats.pearsonr(e[ok], age[ok])
    return float(r), float(p)


def age_vs_null(r: float, null: np.ndarray, min_null: int = 20) -> tuple[float, float]:
    """(z, two-sided empirical p) of an observed age correlation against same-weight random
    projections: how far the module's age association sits from what the cohort's shared
    structure gives any projection with these weights."""
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not np.isfinite(r) or len(null) < min_null:
        return np.nan, np.nan
    mu, sd = float(null.mean()), float(null.std())
    z = (r - mu) / sd if sd > 1e-12 else np.nan
    p = float((1 + np.sum(np.abs(null - mu) >= abs(r - mu))) / (1 + len(null)))
    return float(z), p


class Cohort:
    """Row-standardized production features and ages for one cohort x region."""

    def __init__(self, cohort: str, region: str):
        from isograph.io.artifacts import load_dataset_bundle

        fs = pd.read_parquet(region_store(cohort, region, "isograph_vae", "feature_scores.parquet"))
        self.samples = [c for c in fs.columns if c not in _META]
        self.info = pd.DataFrame({"gene": _bare(fs["gene_id"]).to_numpy(),
                                  "feature_type": fs["feature_type"].astype(str).to_numpy()})
        self.x = _zrows(fs[self.samples].to_numpy(float))
        self.index = {(g, t): i for i, (g, t) in enumerate(zip(self.info["gene"],
                                                               self.info["feature_type"]))}
        st = load_dataset_bundle(rel(*BUNDLES[cohort], region)).sample_table
        st = st.set_index(st["sample_id"].astype(str))
        self.age = pd.to_numeric(st.loc[self.samples, AGE_COL[cohort]], errors="coerce").to_numpy(float)
        self.cohort, self.region = cohort, region


def project_modules(src: Cohort, tgt: Cohort, modules: dict[str, set], trusted: set,
                    method: str, align: pd.DataFrame | None, n_perm: int, seed: int,
                    sig: float = 0.05) -> pd.DataFrame:
    types = FEATURE_TYPES[method]
    sign = {} if align is None else dict(zip(align.loc[align["usable"], "gene"],
                                             align.loc[align["usable"], "sign"]))
    # target universe for the null, by feature type; switch rows only where alignable
    pool = {t: [i for (g, ft), i in tgt.index.items()
                if ft == t and (t != "switch" or g in sign)] for t in types}
    rng = np.random.default_rng(seed)
    rows = []
    for mid in sorted(trusted):
        genes = {g.split(".")[0] for g in modules.get(mid, set())}
        feats, dropped = [], 0
        for g in sorted(genes):
            for t in types:
                if (g, t) not in src.index or (g, t) not in tgt.index:
                    continue
                if t == "switch" and g not in sign:
                    dropped += 1
                    continue
                feats.append((g, t))
        if len(feats) < MIN_FEATURES:
            rows.append({"module_id": mid, "n_features": len(feats), "status": "too_few_features",
                         "n_switch_unalignable": dropped})
            continue
        si = [src.index[f] for f in feats]
        ti = [tgt.index[f] for f in feats]
        flip = np.array([sign.get(g, 1.0) if t == "switch" else 1.0 for g, t in feats])
        xs, xt = src.x[si], tgt.x[ti] * flip[:, None]
        w = frozen_weights(xs)
        es, et = project(w, xs), project(w, xt)
        kme = signed_kme(w, xt, et)
        # Null: the same frozen weights on type-matched random features, in BOTH cohorts. It
        # serves two statistics. kME: a cohort with strong shared structure gives random feature
        # sets a high signed kME too. Age: that shared structure (degradation, RIN, ischemic time)
        # is itself age-correlated, so ANY weighted projection inherits one age sign -- the raw
        # projected age r then reports the cohort, not the module. Module-specific aging is the
        # observed age r relative to the same-weight random projections.
        n_by_type = pd.Series([t for _, t in feats]).value_counts().to_dict()
        src_pool = {t: [src.index[k] for k in tgt_keys if k in src.index]
                    for t, tgt_keys in (
                        (t, [key for key in tgt.index if key[1] == t and (t != "switch" or key[0] in sign)])
                        for t in types)}
        null_kme = np.empty(n_perm)
        null_age_t = np.empty(n_perm)
        null_age_s = np.empty(n_perm)
        for k in range(n_perm):
            idx_t = np.concatenate([rng.choice(pool[t], size=n, replace=False)
                                    for t, n in n_by_type.items()])
            xr = tgt.x[idx_t]
            er = project(w, xr)
            null_kme[k] = signed_kme(w, xr, er)
            null_age_t[k] = _age_r(er, tgt.age)[0]
            idx_s = np.concatenate([rng.choice(src_pool[t], size=n, replace=False)
                                    for t, n in n_by_type.items()])
            null_age_s[k] = _age_r(project(w, src.x[idx_s]), src.age)[0]
        native = np.linalg.svd(xt.T - xt.T.mean(axis=0), full_matrices=False)[0][:, 0]
        pres_r = abs(float(np.corrcoef(native, et)[0, 1]))
        r_s, p_s = _age_r(es, src.age)
        r_t, p_t = _age_r(et, tgt.age)
        z_s, pz_s = age_vs_null(r_s, null_age_s)
        z_t, pz_t = age_vs_null(r_t, null_age_t)
        rows.append({
            "module_id": mid, "status": "projected", "n_features": len(feats),
            "n_switch_features": int(sum(t == "switch" for _, t in feats)),
            "n_switch_unalignable": dropped,
            "signed_kme": kme, "kme_null_mean": float(null_kme.mean()),
            "kme_null_sd": float(null_kme.std()),
            "kme_perm_p": float((1 + np.sum(null_kme >= kme)) / (1 + n_perm)),
            "preservation_r_native_pc1": pres_r,
            # raw frozen-eigengene age correlations -- confounded by cohort-wide structure
            "age_r_source": r_s, "age_p_source": p_s,
            "age_r_target": r_t, "age_p_target": p_t,
            "raw_sign_match": bool(np.isfinite(r_s) and np.isfinite(r_t) and np.sign(r_s) == np.sign(r_t)),
            "raw_both_sig": bool(np.isfinite(p_s) and np.isfinite(p_t) and p_s < sig and p_t < sig),
            # module-specific aging: age r relative to same-weight random projections
            "age_null_mean_source": float(np.nanmean(null_age_s)),
            "age_null_mean_target": float(np.nanmean(null_age_t)),
            "age_z_source": z_s, "age_perm_p_source": pz_s,
            "age_z_target": z_t, "age_perm_p_target": pz_t,
            "sign_match": bool(np.isfinite(z_s) and np.isfinite(z_t) and np.sign(z_s) == np.sign(z_t)),
            "both_sig": bool(np.isfinite(pz_s) and np.isfinite(pz_t) and pz_s < sig and pz_t < sig),
        })
    out = pd.DataFrame(rows)
    ok = out["status"] == "projected"
    if ok.any():
        out.loc[ok, "kme_perm_q"] = stats.false_discovery_control(out.loc[ok, "kme_perm_p"], method="bh")
    return out


def run(pair: str, method: str, n_perm: int, seed: int) -> None:
    (bc, br), (gc, gr) = REGION_PAIRS[pair]
    cohorts = {("brainseq", br): Cohort(bc, br), ("gtex", gr): Cohort(gc, gr)}
    align = None
    if method == "isograph":
        cache = _out_dir() / f"switch_axis_alignment__{pair}.parquet"
        if cache.exists():
            align = pd.read_parquet(cache)
        else:
            align = align_switch_axes(switch_loadings(bc, br), switch_loadings(gc, gr))
            align.to_parquet(cache, index=False)
        print(f"[{pair}] switch axes: {len(align):,} shared genes, "
              f"{int(align['usable'].sum()):,} orientable, "
              f"{int((align['usable'] & (align['sign'] < 0)).sum()):,} flipped", flush=True)
    for (s_key, t_key) in ((("brainseq", br), ("gtex", gr)), (("gtex", gr), ("brainseq", br))):
        src, tgt = cohorts[s_key], cohorts[t_key]
        # alignment orients GTEx relative to BrainSEQ; the reverse direction uses the same sign
        modules = _load_production_modules(src.cohort, src.region, method)
        trusted = _trusted_set(src.cohort, src.region, method)
        res = project_modules(src, tgt, modules, trusted, method, align, n_perm, seed)
        direction = f"{src.cohort}_to_{tgt.cohort}"
        res.insert(0, "direction", direction)
        res.insert(0, "method", method)
        res.insert(0, "pair", pair)
        path = _out_dir() / f"eigengene_projection__{pair}__{method}__{direction}.parquet"
        res.to_parquet(path, index=False)
        p = res[res["status"] == "projected"]
        print(f"[{pair}/{method}/{direction}] {len(p)}/{len(res)} projected | preserved q<0.05 "
              f"{int((p['kme_perm_q'] < 0.05).sum()) if len(p) else 0} | sign match "
              f"{int(p['sign_match'].sum())} | both sig {int(p['both_sig'].sum())}", flush=True)


def summarize(df: pd.DataFrame, sig: float = 0.05) -> pd.DataFrame:
    rows = []
    for (method, direction), g in df[df["status"] == "projected"].groupby(["method", "direction"]):
        n = len(g)
        valid = g[g["age_z_source"].notna() & g["age_z_target"].notna()]
        nv = len(valid)
        n_match = int(valid["sign_match"].sum())
        src_sig = valid[valid["age_perm_p_source"] < sig]
        rho = stats.spearmanr(valid["age_z_source"], valid["age_z_target"]) if nv >= 5 else None
        rows.append({
            "method": method, "direction": direction, "n_modules": n,
            "n_preserved_q05": int((g["kme_perm_q"] < 0.05).sum()),
            "frac_preserved_q05": float((g["kme_perm_q"] < 0.05).mean()),
            "median_signed_kme": float(g["signed_kme"].median()),
            "median_kme_null": float(g["kme_null_mean"].median()),
            "median_preservation_r": float(g["preservation_r_native_pc1"].median()),
            # headline: module-specific aging, each cohort's age r centred on its own
            # same-weight random projections
            "n_age_testable": nv,
            "sign_match": n_match,
            "sign_match_binom_p": float(stats.binomtest(n_match, nv, 0.5, alternative="greater").pvalue)
            if nv else np.nan,
            "n_source_age_sig": len(src_sig),
            "sign_match_among_source_sig": int(src_sig["sign_match"].sum()),
            "n_both_sig": int(valid["both_sig"].sum()),
            "spearman_age_z": float(rho.statistic) if rho is not None else np.nan,
            "spearman_age_z_p": float(rho.pvalue) if rho is not None else np.nan,
            # raw projected age r, kept to show the cohort-wide confound it carries
            "raw_sign_match": int(g["raw_sign_match"].sum()),
            "raw_both_sig": int(g["raw_both_sig"].sum()),
            "median_age_null_target": float(g["age_null_mean_target"].median()),
        })
    return pd.DataFrame(rows)


def aggregate() -> None:
    d = _out_dir()
    parts = [pd.read_parquet(p) for p in sorted(d.glob("eigengene_projection__*.parquet"))]
    if not parts:
        raise SystemExit(f"no projection outputs under {d}")
    df = pd.concat(parts, ignore_index=True)
    df.to_parquet(d / "eigengene_projection_all.parquet", index=False)
    s = summarize(df)
    s.to_parquet(d / "eigengene_projection_summary.parquet", index=False)
    aligns = [pd.read_parquet(p).assign(pair=p.stem.split("__")[1])
              for p in sorted(d.glob("switch_axis_alignment__*.parquet"))]
    lines = [
        "# Cross-cohort eigengene projection",
        "",
        "Each trusted module's eigengene weights are frozen in the source cohort and applied "
        "unchanged to the target cohort's features. **Preservation** is the mean signed kME "
        "against a null that applies the same weights to size- and type-matched random target "
        "features (BH across modules). **Aging** is module-specific: in each cohort the frozen "
        "eigengene's age correlation is expressed as a z against the same-weight random "
        "projections, and the two cohorts' z are compared. Pooled over the three matched region "
        "pairs.",
        "",
        "**Why the raw age correlation is not the statistic.** Both cohorts carry cohort-wide "
        "structure (RNA quality, ischemic time, composition) that is itself age-correlated, so any "
        "weighted projection inherits one age sign regardless of the module. On the raw "
        "correlations the three BrainSEQ→GTEx pairs gave sign agreement of 44/44, 1/50 and 32/36 — "
        "a property of each target cohort, not of the modules. The raw counts are kept in the last "
        "columns to show the confound.",
        "",
        "| method | direction | modules | preserved q<0.05 | median signed kME (null) | "
        "median |r| vs native PC1 | age-z sign match | among source-age-sig | both sig | "
        "Spearman age z | raw sign match | raw both sig |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in s.itertuples(index=False):
        lines.append(
            f"| {r.method} | {r.direction} | {r.n_modules} | {r.n_preserved_q05} "
            f"({r.frac_preserved_q05:.2f}) | {r.median_signed_kme:.3f} ({r.median_kme_null:.3f}) | "
            f"{r.median_preservation_r:.3f} | {r.sign_match}/{r.n_age_testable} "
            f"(p = {r.sign_match_binom_p:.3g}) | {r.sign_match_among_source_sig}/{r.n_source_age_sig} | "
            f"{r.n_both_sig} | {r.spearman_age_z:.3f} (p = {r.spearman_age_z_p:.3g}) | "
            f"{r.raw_sign_match}/{r.n_modules} | {r.raw_both_sig} |")
    if aligns:
        a = pd.concat(aligns, ignore_index=True)
        lines += ["", "## Switch-axis orientation", "",
                  "| pair | shared genes | orientable | flipped | median |cosine| |",
                  "|---|---|---|---|---|"]
        for pair, g in a.groupby("pair"):
            lines.append(f"| {pair} | {len(g):,} | {int(g['usable'].sum()):,} | "
                         f"{int((g['usable'] & (g['sign'] < 0)).sum()):,} | "
                         f"{g['cosine'].abs().median():.3f} |")
    lines += ["",
              "**Scope.** Granularity still matters: WGCNA's larger modules average more features, "
              "which raises kME stability mechanically, so compare each method with its own null "
              "rather than the two methods' raw preservation rates. No covariates enter the age "
              "test, matching the published linear arm; cohort and quantifier remain confounded.",
              ""]
    (d / "EIGENGENE_PROJECTION.md").write_text("\n".join(lines))
    print(s.to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--pair", choices=list(REGION_PAIRS), required=True)
    r.add_argument("--method", choices=list(FEATURE_TYPES), default="isograph")
    r.add_argument("--n-perm", type=int, default=1000)
    r.add_argument("--seed", type=int, default=13)
    sub.add_parser("aggregate")
    args = ap.parse_args()
    if args.cmd == "run":
        run(args.pair, args.method, args.n_perm, args.seed)
    else:
        aggregate()


if __name__ == "__main__":
    main()
