"""Module context and held-out DTU evidence (PI review §7 item 12a)

The question
------------
Gene-wise DTU and IsoGraph answer complementary questions. satuRn asks which individual genes
show differential transcript usage; IsoGraph asks whether those gene-level effects are organized
into coordinated programs. ``isa_concordance`` shows the two agree. This CLI tests whether the
organization carries information that gene-wise DTU statistics do not, on held-out donors:

    Does the DTU evidence of a gene's module neighbours in one donor subset predict that gene's
    DTU evidence in an independent donor subset, after conditioning on its own discovery-subset
    evidence?

If yes, coordinated switching structure contains reproducible information beyond each gene's own
statistic. If no, the modules restate gene-wise DTU. Nothing here asks whether satuRn is adequate
at its own job, and results should not be framed that way.

Design (pre-registered 2026-09-18, before any result was seen)
---------------------------------------------------------------
* Unit: gene. Replicate: one split-half seed x one direction (A->B, B->A); 5 seeds x 2 = 10
  replicates per region. The halves are EXACTLY the stage-04 split halves
  (``stability._split_indices(n, SEED_BASE + k)`` on the same bundle and filter), so the
  IsoGraph partition used for discovery half A was fit on half A's samples alone. The
  replication half never touches the partition or the context.
* Per-gene DTU: satuRn (the test IsoformSwitchAnalyzeR v2 runs), continuous age_z with the
  lean covariates, staged exactly as ``isa_concordance``, run separately in each half. Gene
  evidence e = -log10(Simes p over the gene's tested isoforms); Simes rather than min p, so e is
  not inflated by isoform count.
* Regions: the six stage-04 split-half regions (BrainSEQ caudate / hippocampus / dlpfc; GTEx
  caudate_basal_ganglia / hippocampus / frontal_cortex_ba9), production resolution.
* Genes: satuRn-tested in both halves AND assigned to a discovery-half IsoGraph module that
  holds >= 5 such genes.
* Module context: leave-one-out mean of the module-mates' discovery evidence. A gene's own
  evidence never enters its context.
* Model: e_rep ~ ns(e_disc, df=4) + log n_tx + log mean count + minor-isoform usage + context.
  Statistic: partial correlation of context with e_rep given the rest, averaged over the 10
  replicates.
* Primary null: module labels permuted within technical strata (deciles of log mean count x
  tertiles of minor-isoform usage). This keeps each module's size AND technical make-up, so a
  module that is only a bin of well-measured genes cannot pass. 1,000 permutations, each
  replicate permuted independently and the replicate mean recomputed per permutation;
  one-sided p = (1 + #null >= obs) / (1 + N). BH across the six regions.
* Decision rule: SUPPORTED if BH q < 0.05 in >= 4 of 6 regions with at least one per cohort;
  PARTIAL if 1-3 regions; NOT SUPPORTED if none.

Direction of the known bias: module-mates share discovery-half sampling noise with the gene.
Given e_disc, a high context then says the gene's own e_disc is partly noise, which pushes the
context coefficient NEGATIVE. The test is conservative, not anti-conservative.

Secondary (reported, never the verdict)
---------------------------------------
* Unstratified (size-only) permutation p.
* Abundance co-expression comparator: context = mean discovery evidence of the gene's 50 most
  positively correlated genes on covariate-residualized log-CPM in the discovery half; the
  IsoGraph context is then tested GIVEN that comparator (joint model, same stratified null).
  Answers "would any co-expression neighbourhood do?". The committed WGCNA split-half
  partitions cannot serve here: ``stability_wgcna.R`` draws its halves with R's RNG, so its
  "half A" is not IsoGraph's half A and would leak replication samples into the grouping.
* Network context: |weight|-weighted mean of graph neighbours' discovery evidence on the saved
  half graph, to separate the module from the immediate neighbourhood.
* Sub-threshold enrichment (descriptive; the continuous test above is the inference): among
  genes that do not individually meet the discovery DTU criterion (BH q >= 0.05 over the genes
  tested in both halves), the fraction reaching p < 0.05 in the replication half, top vs bottom
  tertile of module context within that set, with an OR adjusted for the gene's own discovery
  evidence and the technical covariates. Both tertiles are selected by the same criterion, so
  regression to the mean acts on them alike.

Reporting (added 2026-09-18, after the first run; the verdict rule is unchanged): effect sizes
lead -- the mean partial r across the ten split directions with its median, range and number of
positive directions -- and the permutation q follows. The permutation null conditions on the
partition, so it does not measure split-to-split variation; the direction distribution does.
``giant-sensitivity`` (post hoc, supplement) recomputes the statistic without modules of >= 900
genes and without each split's largest module.

Outputs (``04_module_trust/_m/dtu_added_value/``):
  * ``halves/<cohort>__<region>__seed<k>__<half>/`` -- satuRn results, gene evidence and
    abundance neighbours per half (gitignored cache; regenerable with ``saturn``).
  * ``<cohort>__<region>__replicates.parquet`` / ``__null.parquet`` -- per-replicate statistics
    and permutation nulls, stamped with the partition fingerprint of every half used.
  * ``giant_module_sensitivity.parquet`` -- from ``giant-sensitivity`` (supplement).
  * ``region_summary.parquet``, ``summary.json``, ``DTU_ADDED_VALUE.md`` -- from ``summarize``.

Usage:
    python -m isograph_benchmark.real_data.module_dtu_added_value saturn \\
        --cohort brainseq --region caudate --seed 0 --cores 8
    python -m isograph_benchmark.real_data.module_dtu_added_value analyze \\
        --cohort brainseq --region caudate
    python -m isograph_benchmark.real_data.module_dtu_added_value giant-sensitivity
    python -m isograph_benchmark.real_data.module_dtu_added_value summarize
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from isograph_benchmark.paths import ensure_dir, rel, stage_out

REGIONS = {
    "brainseq": ("caudate", "hippocampus", "dlpfc"),
    "gtex": ("caudate_basal_ganglia", "hippocampus", "frontal_cortex_ba9"),
}
SEEDS = 5
HALVES = ("A", "B")
DIRECTIONS = (("A", "B"), ("B", "A"))  # (discovery, replication)
MIN_MODULE_GENES = 5
SPLINE_DF = 4
N_ABUNDANCE_NEIGHBOURS = 50
N_PERM = 1000
PERM_SEED = 13
ALPHA = 0.05
EXPR_BINS = 10
USAGE_BINS = 3
MIN_REGIONS_SUPPORTED = 4
_P_FLOOR = 1e-300


# --------------------------------------------------------------------------- #
# paths
# --------------------------------------------------------------------------- #
def _root() -> Path:
    return ensure_dir(stage_out("trust", "dtu_added_value"))


def _half_dir(cohort: str, region: str, seed: int, half: str) -> Path:
    return _root() / "halves" / f"{cohort}__{region}__seed{seed}__{half}"


def _strip_ver(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace(r"\.\d+$", "", regex=True)


# --------------------------------------------------------------------------- #
# per-gene evidence and covariates
# --------------------------------------------------------------------------- #
def simes_by_gene(res: pd.DataFrame) -> pd.DataFrame:
    """Gene-level Simes p over each gene's tested isoforms, plus min p and isoform count.

    Simes, not min p: min p shrinks with the number of isoforms tested, so a min-p evidence
    score would partly measure isoform count, which modules can group by.
    """
    r = res[["gene_id", "pval"]].copy()
    r["gene"] = _strip_ver(r["gene_id"])
    r["pval"] = pd.to_numeric(r["pval"], errors="coerce")
    r = r.dropna(subset=["pval"]).sort_values(["gene", "pval"], kind="mergesort")
    r["i"] = r.groupby("gene").cumcount() + 1
    r["n"] = r.groupby("gene")["pval"].transform("size")
    r["s"] = r["n"] * r["pval"] / r["i"]
    g = r.groupby("gene")
    out = pd.DataFrame({
        "n_tx": g["pval"].size(),
        "min_p": g["pval"].min(),
        "simes_p": g["s"].min().clip(upper=1.0),
    }).reset_index()
    out["evidence"] = -np.log10(out["simes_p"].clip(lower=_P_FLOOR))
    return out


def gene_covariates(counts: np.ndarray, tx_gene: pd.Series) -> pd.DataFrame:
    """Technical covariates per gene from the tested isoforms of one half.

    ``log_mean_count``: log1p of the mean per-sample gene total. ``minor_usage``: 1 minus the
    mean usage of the gene's dominant isoform (samples with a zero gene total skipped) -- how
    much room the gene has to switch, i.e. how estimable its usage is.
    """
    genes = _strip_ver(tx_gene).to_numpy()
    codes, uniq = pd.factorize(genes)
    agg = sparse.csr_matrix(
        (np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
    counts = np.asarray(counts, dtype=float)
    totals = np.asarray(agg @ counts)
    per_tx_total = totals[codes]
    with np.errstate(invalid="ignore", divide="ignore"):
        usage = np.where(per_tx_total > 0, counts / per_tx_total, np.nan)
    mean_usage = np.nanmean(usage, axis=1)
    top = pd.Series(mean_usage).groupby(codes).max().to_numpy()
    return pd.DataFrame({
        "gene": uniq,
        "log_mean_count": np.log1p(totals.mean(axis=1)),
        "minor_usage": 1.0 - top,
    })


def _covariate_matrix(sample_table: pd.DataFrame, covariates: list[str]) -> np.ndarray:
    """Numeric, mean-imputed discovery covariates with an intercept (constant columns dropped)."""
    cols = [np.ones(len(sample_table))]
    for c in covariates:
        if c not in sample_table.columns:
            continue
        x = pd.to_numeric(sample_table[c], errors="coerce").to_numpy(dtype=float)
        if np.all(np.isnan(x)):
            continue
        x = np.where(np.isnan(x), np.nanmean(x), x)
        if np.nanstd(x) > 0:
            cols.append(x)
    return np.column_stack(cols)


def abundance_neighbours(logcpm: np.ndarray, genes: np.ndarray, covars: np.ndarray,
                         k: int = N_ABUNDANCE_NEIGHBOURS, chunk: int = 2000) -> pd.DataFrame:
    """Each gene's k most positively correlated genes on covariate-residualized log-CPM.

    ``logcpm`` is genes x samples. Signed (positive correlation only), as the production
    WGCNA baselines are signed networks.
    """
    beta, *_ = np.linalg.lstsq(covars, logcpm.T, rcond=None)
    resid = logcpm.T - covars @ beta                     # samples x genes
    resid -= resid.mean(axis=0)
    sd = resid.std(axis=0)
    keep = sd > 0
    z = (resid[:, keep] / sd[keep]).astype(np.float32)
    g = genes[keep]
    n = z.shape[0]
    k = min(k, z.shape[1] - 1)
    rows = []
    for start in range(0, z.shape[1], chunk):
        block = (z[:, start:start + chunk].T @ z) / n    # chunk x genes
        for j in range(block.shape[0]):
            block[j, start + j] = -np.inf                # never your own neighbour
        top = np.argpartition(-block, k, axis=1)[:, :k]
        r = np.take_along_axis(block, top, axis=1)
        src = np.repeat(g[start:start + block.shape[0]], k)
        rows.append(pd.DataFrame({"gene": src, "neighbour": g[top.ravel()], "r": r.ravel()}))
    return pd.concat(rows, ignore_index=True)


# --------------------------------------------------------------------------- #
# step 1: satuRn + covariates + abundance neighbours, per half
# --------------------------------------------------------------------------- #
def split_indices(n: int, seed: int, half: str) -> np.ndarray:
    """The stage-04 split-half sample positions for (seed, half)."""
    from isograph_benchmark.real_data import stability
    a, b = stability._split_indices(n, stability.SEED_BASE + seed)
    return a if half == "A" else b


def _load_region(cohort: str, region: str) -> dict:
    """Bundle, filtered exactly as the split-half fits were (stability.fit_isograph)."""
    from isograph.io.artifacts import load_dataset_bundle
    from isograph_benchmark.real_data import stability
    from isograph_benchmark.real_data.run_models import filter_production_transcripts
    spec = stability.COHORTS[cohort]
    bundle = load_dataset_bundle(rel(*spec["bundle_root"], region))
    st = bundle.sample_table.reset_index(drop=True)
    tc, tt = bundle.matrices["transcript_counts"], bundle.feature_tables["transcript"]
    if spec["filter_transcripts"]:
        tc, tt = filter_production_transcripts(tc, tt)
    return {
        "sample_table": st, "tc": np.asarray(tc), "tt": tt.reset_index(drop=True),
        "gene_counts": np.asarray(bundle.matrices["gene_counts"]),
        "gene_ids": _strip_ver(bundle.feature_tables["gene"]["gene_id"]).to_numpy(),
        "covariates": list(spec["covariates"]),
    }


def run_saturn_half(cohort: str, region: str, seed: int, half: str, cores: int,
                    data: dict | None = None, force: bool = False) -> None:
    from isograph_benchmark.real_data import isa_concordance as isa
    out = _half_dir(cohort, region, seed, half)
    ev_path, nb_path = out / "gene_evidence.parquet", out / "abundance_neighbours.parquet"
    if ev_path.exists() and nb_path.exists() and not force:
        print(f"[{cohort}/{region}] seed{seed} {half}: exists, skipping", flush=True)
        return
    ensure_dir(out)
    data = data or _load_region(cohort, region)
    st = data["sample_table"]
    idx = split_indices(len(st), seed, half)
    st_h = st.iloc[idx].reset_index(drop=True)

    exposure = isa._build_exposure(st_h, cohort, "age")
    res = isa._stage_and_run_saturn(cohort, region, "age", out, st_h, data["tc"][:, idx],
                                    data["tt"], exposure, cores, 0)
    shutil.rmtree(out / "_saturn_inputs", ignore_errors=True)  # staged copy of the counts

    ev = simes_by_gene(res)
    tested = data["tt"]["transcript_id"].isin(set(res["isoform_id"])).to_numpy()
    cov = gene_covariates(data["tc"][tested][:, idx], data["tt"].loc[tested, "gene_id"])
    ev = ev.merge(cov, on="gene", how="left")
    ev["cohort"], ev["region"], ev["seed"], ev["half"] = cohort, region, seed, half
    ev["n_samples"] = len(idx)
    ev.to_parquet(ev_path, index=False, compression="zstd")

    # Abundance comparator, on the same half and the same (satuRn-tested) gene universe.
    gcounts = data["gene_counts"][:, idx].astype(float)
    lib = gcounts.sum(axis=0)
    in_u = np.isin(data["gene_ids"], ev["gene"].to_numpy())
    logcpm = np.log2(gcounts[in_u] / lib * 1e6 + 1.0)
    nb = abundance_neighbours(logcpm, data["gene_ids"][in_u],
                              _covariate_matrix(st_h, data["covariates"]))
    nb.to_parquet(nb_path, index=False, compression="zstd")
    print(f"[{cohort}/{region}] seed{seed} {half}: {len(ev)} genes, {len(idx)} samples, "
          f"median evidence {ev['evidence'].median():.2f} -> {out}", flush=True)


# --------------------------------------------------------------------------- #
# step 2: held-out context test
# --------------------------------------------------------------------------- #
def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def loo_context(labels: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Leave-one-out mean of ``values`` over each gene's module-mates."""
    codes, _ = pd.factorize(labels)
    sums = np.bincount(codes, weights=values)
    counts = np.bincount(codes)
    denom = counts[codes] - 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(denom > 0, (sums[codes] - values) / denom, np.nan)


def loo_context_matrix(label_matrix: np.ndarray, values: np.ndarray) -> np.ndarray:
    """``loo_context`` for every column of an (n_genes x n_perm) label-code matrix."""
    n, m = label_matrix.shape
    k = int(label_matrix.max()) + 1
    flat = label_matrix + (np.arange(m) * k)[None, :]
    sums = np.bincount(flat.ravel(), weights=np.repeat(values, m), minlength=k * m)
    counts = np.bincount(flat.ravel(), minlength=k * m)
    denom = counts[flat] - 1
    return (sums[flat] - values[:, None]) / denom


def technical_strata(log_mean_count: np.ndarray, minor_usage: np.ndarray,
                     expr_bins: int = EXPR_BINS, usage_bins: int = USAGE_BINS) -> np.ndarray:
    e = pd.qcut(pd.Series(log_mean_count).rank(method="first"), expr_bins, labels=False)
    u = pd.qcut(pd.Series(minor_usage).rank(method="first"), usage_bins, labels=False)
    return (e.to_numpy() * usage_bins + u.to_numpy()).astype(int)


def permute_within(labels: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Shuffle ``labels`` among positions of the same stratum (module sizes and each module's
    stratum make-up are preserved exactly)."""
    base = np.argsort(strata, kind="stable")
    shuffled = np.lexsort((rng.random(len(labels)), strata))
    out = np.empty_like(labels)
    out[base] = labels[shuffled]
    return out


def design_matrix(d: pd.DataFrame, extra: list[str] | None = None) -> np.ndarray:
    import patsy
    # Centred natural-spline basis + explicit intercept: an uncentred cr() basis already spans
    # the constant, and a rank-deficient design would silently distort the residualization.
    spline = np.asarray(patsy.dmatrix(f"cr(x, df={SPLINE_DF}, constraints='center')",
                                      {"x": d["e_disc"].to_numpy()}, return_type="matrix"))
    cols = [spline, np.log(d["n_tx"].to_numpy(dtype=float))[:, None],
            d["log_mean_count"].to_numpy(dtype=float)[:, None],
            d["minor_usage"].to_numpy(dtype=float)[:, None]]
    for c in extra or []:
        cols.append(d[c].to_numpy(dtype=float)[:, None])
    return np.hstack(cols)


def _basis(X: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Orthonormal basis of col(X), rank-revealing (QR on a rank-deficient X would return extra,
    arbitrary directions and residualize them away too)."""
    u, sv, _ = np.linalg.svd(X, full_matrices=False)
    return u[:, sv > tol * sv.max()]


def _residualizer(X: np.ndarray):
    q = _basis(X)
    return lambda v: v - q @ (q.T @ v)


def partial_r(resid_x: np.ndarray, resid_y: np.ndarray) -> np.ndarray:
    """Correlation of each column of ``resid_x`` (or a vector) with ``resid_y``."""
    rx = resid_x - resid_x.mean(axis=0)
    ry = resid_y - resid_y.mean()
    num = ry @ rx
    den = np.sqrt((rx ** 2).sum(axis=0) * (ry ** 2).sum())
    return num / den


def build_replicate_frame(disc: pd.DataFrame, rep: pd.DataFrame, partition: pd.DataFrame,
                          neighbours: pd.DataFrame, edges: pd.DataFrame | None) -> pd.DataFrame:
    """One discovery->replication frame: genes tested in both halves and module-assigned."""
    u = disc[["gene", "evidence", "simes_p", "n_tx", "log_mean_count", "minor_usage"]].rename(
        columns={"evidence": "e_disc", "simes_p": "p_disc"})
    u = u.merge(rep[["gene", "evidence", "simes_p"]].rename(
        columns={"evidence": "e_rep", "simes_p": "p_rep"}), on="gene", how="inner")
    u = u.dropna(subset=["e_disc", "e_rep", "log_mean_count", "minor_usage"])
    u["q_disc"] = bh(u["p_disc"].to_numpy())
    e_disc = u.set_index("gene")["e_disc"]

    nb = neighbours[neighbours["neighbour"].isin(e_disc.index)]
    ctx_abund = nb.assign(e=nb["neighbour"].map(e_disc)).groupby("gene")["e"].mean()
    u["ctx_abund"] = u["gene"].map(ctx_abund)

    if edges is not None and len(edges):
        e = edges.assign(source=_strip_ver(edges["source"]), target=_strip_ver(edges["target"]),
                         w=edges["weight"].abs())
        both = pd.concat([e[["source", "target", "w"]],
                          e.rename(columns={"source": "target", "target": "source"})[
                              ["source", "target", "w"]]], ignore_index=True)
        both = both[both["target"].isin(e_disc.index)]
        both["we"] = both["w"] * both["target"].map(e_disc)
        g = both.groupby("source")
        u["ctx_net"] = u["gene"].map(g["we"].sum() / g["w"].sum())
    else:
        u["ctx_net"] = np.nan

    p = partition.assign(gene=_strip_ver(partition["gene_id"]))[["gene", "module_id"]]
    d = u.merge(p, on="gene", how="inner")
    size = d["module_id"].map(d["module_id"].value_counts())
    d = d[size >= MIN_MODULE_GENES].reset_index(drop=True)
    d["ctx_module"] = loo_context(d["module_id"].to_numpy(), d["e_disc"].to_numpy())
    d.attrs["n_universe"] = len(u)
    return d


def replicate_stats(d: pd.DataFrame, n_perm: int, rng: np.random.Generator) -> tuple[dict, dict]:
    """Observed statistics and permutation nulls for one replicate frame."""
    y = d["e_rep"].to_numpy(dtype=float)
    e = d["e_disc"].to_numpy(dtype=float)
    labels, _ = pd.factorize(d["module_id"])
    strata = technical_strata(d["log_mean_count"].to_numpy(), d["minor_usage"].to_numpy())

    X = design_matrix(d)
    res = _residualizer(X)
    ry = res(y)
    obs = {"n_genes": len(d), "n_universe": int(d.attrs.get("n_universe", np.nan)),
           "n_modules": int(labels.max() + 1)}
    obs["r_module"] = float(partial_r(res(d["ctx_module"].to_numpy()), ry))

    # Abundance comparator: alone, and the module context GIVEN it (same stratified null).
    has_ab = d["ctx_abund"].notna().to_numpy()
    obs["n_missing_abund"] = int((~has_ab).sum())
    ab = d["ctx_abund"].fillna(d["e_disc"].mean())
    obs["r_abund"] = float(partial_r(res(ab.to_numpy()), ry))
    X2 = np.hstack([X, ab.to_numpy()[:, None]])
    res2 = _residualizer(X2)
    ry2 = res2(y)
    obs["r_module_given_abund"] = float(partial_r(res2(d["ctx_module"].to_numpy()), ry2))

    net = d["ctx_net"]
    obs["n_missing_net"] = int(net.isna().sum())
    if net.notna().sum() > 10:
        net = net.fillna(d["e_disc"].mean()).to_numpy()
        obs["r_net"] = float(partial_r(res(net), ry))
        res3 = _residualizer(np.hstack([X, net[:, None]]))
        obs["r_module_given_net"] = float(
            partial_r(res3(d["ctx_module"].to_numpy()), res3(y)))
    else:
        obs["r_net"] = obs["r_module_given_net"] = float("nan")

    # Effect in natural units: change in replication evidence per SD of module context.
    rc = res(d["ctx_module"].to_numpy())
    obs["beta_per_sd"] = float((rc @ ry) / (rc @ rc) * d["ctx_module"].std())

    obs.update(_subthreshold_enrichment(d, X))

    strat = np.column_stack([permute_within(labels, strata, rng) for _ in range(n_perm)])
    plain = np.column_stack([rng.permutation(labels) for _ in range(n_perm)])
    ctx_s = loo_context_matrix(strat, e)
    ctx_p = loo_context_matrix(plain, e)
    null = {
        "r_module_strat": partial_r(ctx_s - _proj(X, ctx_s), ry),
        "r_module_plain": partial_r(ctx_p - _proj(X, ctx_p), ry),
        "r_module_given_abund_strat": partial_r(ctx_s - _proj(X2, ctx_s), ry2),
    }
    return obs, null


def _proj(X: np.ndarray, V: np.ndarray) -> np.ndarray:
    q = _basis(X)
    return q @ (q.T @ V)


def _subthreshold_enrichment(d: pd.DataFrame, X: np.ndarray) -> dict:
    """Replication among genes below the discovery DTU criterion, by module-context tertile.

    Sub-threshold = BH q >= ALPHA in the discovery half. Tertiles are of module context within
    the sub-threshold set; replication = p < ALPHA in the replication half. Descriptive only.
    """
    sub = (d["q_disc"] >= ALPHA).to_numpy()
    out = {"n_disc_sig": int((~sub).sum()), "n_subthreshold": int(sub.sum())}
    if sub.sum() < 30:
        return out
    m = d.loc[sub]
    tert = pd.qcut(m["ctx_module"].rank(method="first"), 3, labels=False).to_numpy()
    hit = (m["p_rep"] < ALPHA).to_numpy()
    out["subthr_rate_top"] = float(hit[tert == 2].mean())
    out["subthr_rate_bottom"] = float(hit[tert == 0].mean())
    out["subthr_n_top"] = int((tert == 2).sum())
    out["subthr_hits_top"] = int(hit[tert == 2].sum())
    try:
        import statsmodels.api as sm
        keep = tert != 1
        Xm = np.hstack([X[sub][keep], (tert[keep] == 2).astype(float)[:, None]])
        fit = sm.Logit(hit[keep].astype(float), Xm).fit(disp=0)
        out["subthr_or_top_vs_bottom"] = float(np.exp(fit.params[-1]))
        out["subthr_or_p"] = float(fit.pvalues[-1])
    except Exception as exc:  # pragma: no cover - separation on tiny strata
        print(f"  sub-threshold logistic skipped: {exc}", flush=True)
    return out


def _load_half(cohort: str, region: str, seed: int, half: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = _half_dir(cohort, region, seed, half)
    ev, nb = d / "gene_evidence.parquet", d / "abundance_neighbours.parquet"
    if not (ev.exists() and nb.exists()):
        raise SystemExit(f"missing satuRn half {d}; run `saturn` first")
    return pd.read_parquet(ev), pd.read_parquet(nb)


def analyze_region(cohort: str, region: str, n_perm: int = N_PERM) -> pd.DataFrame:
    from isograph_benchmark.real_data import stability
    from isograph_benchmark.real_data.partition_provenance import partition_fingerprint
    rng = np.random.default_rng(PERM_SEED)
    rows, nulls = [], []
    for seed in range(SEEDS):
        halves = {h: _load_half(cohort, region, seed, h) for h in HALVES}
        for disc, rep in DIRECTIONS:
            pfile = stability._partitions_dir() / f"isograph__{cohort}__{region}__seed{seed}__{disc}.parquet"
            if not pfile.exists():
                raise SystemExit(f"missing split-half partition {pfile}")
            part = pd.read_parquet(pfile)
            efile = stability._edges_path(cohort, region, "isograph", seed, disc)
            edges = pd.read_parquet(efile) if efile.exists() else None
            d = build_replicate_frame(halves[disc][0], halves[rep][0], part, halves[disc][1], edges)
            obs, null = replicate_stats(d, n_perm, rng)
            obs.update({"cohort": cohort, "region": region, "seed": seed,
                        "discovery": disc, "replication": rep,
                        "partition_fingerprint": partition_fingerprint(part),
                        "leiden_resolution": stability.CANONICAL_LEIDEN_RESOLUTION,
                        "has_edges": edges is not None})
            rows.append(obs)
            nulls.append(pd.DataFrame({"seed": seed, "discovery": disc, "perm": np.arange(n_perm),
                                       **{k: v for k, v in null.items()}}))
            print(f"[{cohort}/{region}] seed{seed} {disc}->{rep}: {obs['n_genes']} genes / "
                  f"{obs['n_modules']} modules | r_module {obs['r_module']:+.3f} "
                  f"(null {null['r_module_strat'].mean():+.3f}) | given abund "
                  f"{obs['r_module_given_abund']:+.3f} | r_abund {obs['r_abund']:+.3f}", flush=True)
    rep_df = pd.DataFrame(rows)
    null_df = pd.concat(nulls, ignore_index=True)
    out = _root()
    rep_df.to_parquet(out / f"{cohort}__{region}__replicates.parquet", index=False, compression="zstd")
    null_df.to_parquet(out / f"{cohort}__{region}__null.parquet", index=False, compression="zstd")
    return rep_df


# --------------------------------------------------------------------------- #
# supplement (post hoc): giant modules
# --------------------------------------------------------------------------- #
GIANT_MODULE_GENES = 900  # the project's giant-module size criterion


def module_context_r(d: pd.DataFrame) -> float:
    """Observed partial r of module context, recomputing the leave-one-out context on ``d``."""
    d = d.reset_index(drop=True)
    ctx = loo_context(d["module_id"].to_numpy(), d["e_disc"].to_numpy())
    res = _residualizer(design_matrix(d))
    return float(partial_r(res(ctx), res(d["e_rep"].to_numpy(dtype=float))))


def giant_sensitivity_row(d: pd.DataFrame, partition: pd.DataFrame) -> dict:
    """The statistic with all genes, without giant modules, and without the largest module.

    Module size is the module's gene count in the whole discovery partition (the size the
    >= 900-gene criterion refers to), not its count among the tested genes.
    """
    size = partition["module_id"].value_counts()
    msize = d["module_id"].map(size)
    largest = size.idxmax()
    no_giant = d[msize < GIANT_MODULE_GENES]
    return {
        "r_all": module_context_r(d),
        "r_no_giant": module_context_r(no_giant) if no_giant["module_id"].nunique() > 1 else float("nan"),
        "r_no_largest": module_context_r(d[d["module_id"] != largest]),
        "n_giant_modules": int((size >= GIANT_MODULE_GENES).sum()),
        "largest_module_genes": int(size.max()),
        "frac_tested_in_giant": float((msize >= GIANT_MODULE_GENES).mean()),
        "frac_tested_in_largest": float((d["module_id"] == largest).mean()),
    }


def giant_sensitivity() -> pd.DataFrame:
    from isograph_benchmark.real_data import stability
    rows = []
    for cohort, regions in REGIONS.items():
        for region in regions:
            for seed in range(SEEDS):
                halves = {h: _load_half(cohort, region, seed, h) for h in HALVES}
                for disc, rep in DIRECTIONS:
                    part = pd.read_parquet(stability._partitions_dir()
                                           / f"isograph__{cohort}__{region}__seed{seed}__{disc}.parquet")
                    d = build_replicate_frame(halves[disc][0], halves[rep][0], part,
                                              halves[disc][1], None)
                    rows.append({"cohort": cohort, "region": region, "seed": seed,
                                 "discovery": disc, **giant_sensitivity_row(d, part)})
    g = pd.DataFrame(rows)
    g.to_parquet(_root() / "giant_module_sensitivity.parquet", index=False, compression="zstd")
    print(g.groupby(["cohort", "region"])[["r_all", "r_no_giant", "r_no_largest"]].mean()
          .round(3).to_string(), flush=True)
    return g


# --------------------------------------------------------------------------- #
# step 3: region summaries, decision rule, report
# --------------------------------------------------------------------------- #
def perm_p(obs: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= obs)) / (1 + len(null)))


def split_distribution(x: pd.Series, prefix: str) -> dict:
    """Mean, median, IQR, range and positive count of a statistic over the split directions."""
    x = x.astype(float)
    return {
        f"{prefix}": float(x.mean()),
        f"{prefix}_median": float(x.median()),
        f"{prefix}_q25": float(x.quantile(0.25)),
        f"{prefix}_q75": float(x.quantile(0.75)),
        f"{prefix}_min": float(x.min()),
        f"{prefix}_max": float(x.max()),
        f"{prefix}_n_positive": int((x > 0).sum()),
    }


def region_summary(rep: pd.DataFrame, null: pd.DataFrame) -> dict:
    """Split-direction distribution and permutation p for one region.

    The permutation statistic is the mean over split directions; its null is the mean over
    directions of each permutation index. The null conditions on the partitions, so it does not
    describe variation between splits -- the distribution columns do.
    """
    per_perm = null.groupby("perm")[["r_module_strat", "r_module_plain",
                                     "r_module_given_abund_strat"]].mean()
    s = per_perm["r_module_strat"]
    obs = rep["r_module"].mean()
    row = {
        "cohort": rep["cohort"].iloc[0], "region": rep["region"].iloc[0],
        "n_replicates": len(rep),
        "n_genes_median": float(rep["n_genes"].median()),
        "n_universe_median": float(rep["n_universe"].median()),
        "coverage": float((rep["n_genes"] / rep["n_universe"]).median()),
        "n_modules_median": float(rep["n_modules"].median()),
        **split_distribution(rep["r_module"], "r_module"),
        **split_distribution(rep["r_module_given_abund"], "r_module_given_abund"),
        **split_distribution(rep["r_abund"], "r_abund"),
        "null_mean": float(s.mean()), "null_sd": float(s.std()),
        "p_strat": perm_p(obs, s.to_numpy()),
        "p_plain": perm_p(obs, per_perm["r_module_plain"].to_numpy()),
        "p_given_abund_strat": perm_p(rep["r_module_given_abund"].mean(),
                                      per_perm["r_module_given_abund_strat"].to_numpy()),
        "r_net": float(rep["r_net"].mean()),
        "r_module_given_net": float(rep["r_module_given_net"].mean()),
        "beta_per_sd": float(rep["beta_per_sd"].mean()),
    }
    for c in ("subthr_rate_top", "subthr_rate_bottom", "n_subthreshold",
              "subthr_n_top", "subthr_hits_top"):
        if c in rep:
            row[c] = float(rep[c].mean())
    if "subthr_or_top_vs_bottom" in rep:
        row["subthr_or_median"] = float(rep["subthr_or_top_vs_bottom"].median())
        row["subthr_or_n_above_1"] = int((rep["subthr_or_top_vs_bottom"] > 1).sum())
    return row


def region_class(q_strat: float, q_given_abund: float) -> str:
    if q_strat >= ALPHA:
        return "not detected"
    if q_given_abund >= ALPHA:
        return "not independent of co-expression"
    return "robust to co-expression"


def verdict(summary: pd.DataFrame) -> str:
    """The pre-registered decision rule (unchanged since 2026-09-18)."""
    sig = summary["q_strat"] < ALPHA
    n_sig = int(sig.sum())
    cohorts = set(summary.loc[sig, "cohort"])
    if n_sig >= MIN_REGIONS_SUPPORTED and cohorts >= set(REGIONS):
        return "SUPPORTED"
    if n_sig >= 1:
        return "PARTIAL"
    return "NOT SUPPORTED"


def summarize() -> pd.DataFrame:
    out = _root()
    rows = []
    for cohort, regions in REGIONS.items():
        for region in regions:
            rp = out / f"{cohort}__{region}__replicates.parquet"
            npth = out / f"{cohort}__{region}__null.parquet"
            if not (rp.exists() and npth.exists()):
                raise SystemExit(f"missing analyze output for {cohort}/{region}")
            rows.append(region_summary(pd.read_parquet(rp), pd.read_parquet(npth)))
    s = pd.DataFrame(rows)
    s["q_strat"] = bh(s["p_strat"].to_numpy())
    s["q_given_abund_strat"] = bh(s["p_given_abund_strat"].to_numpy())
    s["region_class"] = [region_class(a, b) for a, b in zip(s["q_strat"], s["q_given_abund_strat"])]
    v = verdict(s)
    s.to_parquet(out / "region_summary.parquet", index=False, compression="zstd")

    giant_path = out / "giant_module_sensitivity.parquet"
    giant = pd.read_parquet(giant_path) if giant_path.exists() else None
    summary = {
        "question": ("Does the DTU evidence of a gene's module neighbours in one donor subset predict "
                     "that gene's DTU evidence in an independent donor subset, after conditioning on "
                     "its own discovery-subset evidence?"),
        "framing": ("Gene-wise DTU identifies genes whose transcript usage changes; IsoGraph identifies "
                    "coordinated structure among those changes. This tests whether that structure "
                    "carries information not represented by gene-wise DTU statistics alone."),
        "n_regions_robust_to_coexpression": int((s["region_class"] == "robust to co-expression").sum()),
        "n_regions_sig": int((s["q_strat"] < ALPHA).sum()),
        "preregistered_rule": (f"SUPPORTED if BH q_strat < {ALPHA} in >= {MIN_REGIONS_SUPPORTED}/6 "
                               "regions with at least one per cohort; PARTIAL if 1-3; NOT SUPPORTED if 0"),
        "preregistered_verdict": v,
        "n_perm": N_PERM, "seeds": SEEDS, "split_directions": SEEDS * len(DIRECTIONS),
        "min_module_genes": MIN_MODULE_GENES, "spline_df": SPLINE_DF,
        "coexpression_neighbours": N_ABUNDANCE_NEIGHBOURS,
        "strata": f"{EXPR_BINS} log-mean-count bins x {USAGE_BINS} minor-usage bins",
        "giant_sensitivity_present": giant is not None,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    (out / "DTU_ADDED_VALUE.md").write_text(_report(s, summary, giant))
    print(s[["cohort", "region", "r_module", "r_module_median", "r_module_min", "r_module_max",
             "r_module_n_positive", "q_strat", "r_module_given_abund", "q_given_abund_strat",
             "region_class"]].round(3).to_string(index=False))
    print(f"pre-registered verdict: {v}")
    return s


def _dist(r, prefix: str, n: int) -> str:
    g = lambda k: getattr(r, f"{prefix}{k}")
    return (f"{g(''):+.3f} (median {g('_median'):+.3f}; range {g('_min'):+.3f} to {g('_max'):+.3f}) "
            f"| {g('_n_positive')}/{n}")


def _report(s: pd.DataFrame, summary: dict, giant: pd.DataFrame | None) -> str:
    n = summary["split_directions"]
    label = lambda r: f"{'BrainSEQ' if r.cohort == 'brainseq' else 'GTEx'} {r.region}"
    lines = [
        "# Module context and held-out DTU evidence (split halves)",
        "",
        "*Generated by `isograph_benchmark.real_data.module_dtu_added_value summarize`; do not edit.*",
        "",
        f"**Framing.** {summary['framing']}",
        "",
        f"**Question.** {summary['question']}",
        "",
        "**Statistic.** Partial correlation of a gene's leave-one-out module context (the mean "
        "discovery-subset satuRn evidence of its module neighbours) with its satuRn evidence in the "
        "independent subset, given a df-4 spline of its own discovery evidence, log isoform count, log "
        f"mean count and minor-isoform usage. One value per split direction ({SEEDS} seeds x 2 "
        f"directions = {n}); the mean is the permutation statistic.",
        "",
        "## Main result: effect size, split consistency, then permutation q",
        "",
        f"| Analysis | Genes (coverage) | Modules | Partial r, mean (median; range) | Positive directions "
        f"| q | Given co-expression, mean (median; range) | Positive directions | q | Co-expression context alone | Pattern |",
        "|---|---:|---:|---|---:|---:|---|---:|---:|---:|---|",
    ]
    for r in s.itertuples(index=False):
        lines.append(
            f"| {label(r)} | {r.n_genes_median:.0f} ({r.coverage:.0%}) | {r.n_modules_median:.0f} | "
            f"{_dist(r, 'r_module', n)} | {r.q_strat:.2g} | {_dist(r, 'r_module_given_abund', n)} | "
            f"{r.q_given_abund_strat:.2g} | {r.r_abund:+.3f} | {r.region_class} |")

    robust = s[s["region_class"] == "robust to co-expression"]
    weak = s[s["region_class"] == "not independent of co-expression"]
    none = s[s["region_class"] == "not detected"]
    coexp_larger = robust[robust["r_abund"] > robust["r_module"]]
    names = lambda df: ", ".join(label(r) for r in df.itertuples())
    lines += [
        "",
        "**Heterogeneity across analyses is part of the result.** "
        f"Module context predicted held-out DTU evidence beyond co-expression in {len(robust)} of 6 "
        f"analyses ({names(robust)}). "
        + (f"In {names(weak)} the association was small and did not survive conditioning on "
           "co-expression. " if len(weak) else "")
        + (f"It was not detected in {names(none)}. " if len(none) else "")
        + (f"In {names(coexp_larger)}, co-expression context alone predicted more strongly than module "
           "context; there the module adds to co-expression rather than replacing it."
           if len(coexp_larger) else ""),
        "",
        "Permutation null: module labels shuffled within technical strata "
        f"({summary['strata']}), {summary['n_perm']} permutations, BH across the six analyses. It "
        "conditions on each split's partition, so it describes label exchangeability, not variation "
        "between splits; read the range and positive-direction count for that.",
        "",
        "## Sub-threshold enrichment (descriptive)",
        "",
        "Among genes that did not individually meet the discovery DTU criterion (BH q >= 0.05 over genes "
        "tested in both subsets), the fraction reaching p < 0.05 in the independent subset, for the top "
        "vs bottom tertile of module context within that set. Both tertiles pass through the same "
        "discovery criterion, so regression to the mean acts on them alike; the OR also conditions on "
        "each gene's own discovery evidence and the technical covariates. The continuous test above is "
        "the inference; this is its practical reading.",
        "",
        "| Analysis | Sub-threshold genes | Replicating: top vs bottom tertile | Adjusted OR, median | OR > 1 |",
        "|---|---:|---|---:|---:|",
    ]
    for r in s.itertuples(index=False):
        if not hasattr(r, "subthr_rate_top") or pd.isna(r.subthr_rate_top):
            continue
        lines.append(f"| {label(r)} | {r.n_subthreshold:.0f} | {r.subthr_rate_top:.3f} vs "
                     f"{r.subthr_rate_bottom:.3f} | {r.subthr_or_median:.2f} | {r.subthr_or_n_above_1}/{n} |")
    lines += [
        "",
        "## Limitations",
        "",
        "1. **Shared confounding is the major residual limitation.** A technical or biological factor "
        "present in both donor subsets replicates across a split as faithfully as coordinated switching "
        "does. The technical strata and covariates remove expression level and isoform-usage range, not "
        "every such factor. Cross-cohort module preservation (stage 04 step 03f) is the complementary "
        "evidence, because it does not share a cohort's confounds.",
        "2. Coverage: the test concerns genes in a discovery-subset module with >= 5 tested genes "
        "(coverage column); it says nothing about unassigned genes.",
        "3. Module neighbours share discovery-subset sampling noise with the gene; given the gene's own "
        "evidence this pushes the statistic negative, so the test is conservative.",
        "",
        "## Pre-registered decision rule",
        "",
        f"{summary['preregistered_rule']}. Outcome: **{summary['preregistered_verdict']}** "
        f"({summary['n_regions_sig']}/6 analyses at q < {ALPHA}). The rule was fixed before the first "
        "run; the reporting above was reorganized afterwards to lead with effect sizes.",
    ]
    if giant is not None:
        agg = giant.groupby(["cohort", "region"], sort=False)
        lines += [
            "",
            "## Supplement: giant modules (post hoc)",
            "",
            f"The statistic recomputed without modules of >= {GIANT_MODULE_GENES} genes (size in the whole "
            "discovery partition) and without each split's largest module. Mean over split directions, "
            "positive directions in brackets.",
            "",
            "| Analysis | All modules | Without giant modules | Without largest module | Giant modules per split, median | Tested genes in largest module |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        for (cohort, region), g in agg:
            name = f"{'BrainSEQ' if cohort == 'brainseq' else 'GTEx'} {region}"
            cell = lambda c: f"{g[c].mean():+.3f} [{int((g[c] > 0).sum())}/{g[c].notna().sum()}]"
            lines.append(f"| {name} | {cell('r_all')} | {cell('r_no_giant')} | {cell('r_no_largest')} | "
                         f"{g['n_giant_modules'].median():.1f} | {g['frac_tested_in_largest'].median():.0%} |")
    lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s1 = sub.add_parser("saturn", help="satuRn + covariates + co-expression neighbours per half")
    s1.add_argument("--cohort", choices=list(REGIONS), required=True)
    s1.add_argument("--region", required=True)
    s1.add_argument("--seed", type=int, action="append", dest="seeds",
                    help="Split seed(s); repeatable. Default: all.")
    s1.add_argument("--half", choices=HALVES, action="append", dest="halves")
    s1.add_argument("--cores", type=int, default=1)
    s1.add_argument("--force", action="store_true")
    s2 = sub.add_parser("analyze", help="held-out context test for one region")
    s2.add_argument("--cohort", choices=list(REGIONS), required=True)
    s2.add_argument("--region", required=True)
    s2.add_argument("--n-perm", type=int, default=N_PERM)
    sub.add_parser("giant-sensitivity", help="post-hoc supplement: statistic without giant modules")
    sub.add_parser("summarize", help="region summaries, decision rule, report")
    a = p.parse_args()

    if a.cmd in ("saturn", "analyze") and a.region not in REGIONS[a.cohort]:
        raise SystemExit(f"{a.cohort}/{a.region} has no stage-04 split halves "
                         f"(expected one of {REGIONS[a.cohort]})")
    if a.cmd == "saturn":
        data = _load_region(a.cohort, a.region)
        for seed in a.seeds or range(SEEDS):
            for half in a.halves or HALVES:
                run_saturn_half(a.cohort, a.region, seed, half, a.cores, data=data, force=a.force)
    elif a.cmd == "analyze":
        analyze_region(a.cohort, a.region, a.n_perm)
    elif a.cmd == "giant-sensitivity":
        giant_sensitivity()
    else:
        summarize()


if __name__ == "__main__":
    main()
