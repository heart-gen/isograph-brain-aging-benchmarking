"""Per-module trust funnel for real-data IsoGraph modules (see
``03_module_trust/docs/MODULE_TRUST_PLAN.md``).

Q1 (this module, ``stability`` command): which *production* modules are stable enough to
trust? A production module's gene set is scored for how tightly its genes stay co-clustered
across the split-half ensemble (`03_module_trust/_m/stability/partitions/`), relative to a
size-matched permutation null. Trusted = co-assignment density exceeds chance at BH-FDR<0.05.

Primary statistic — **co-assignment density** (the user-chosen Q1 gate): among a module's
genes that are assigned in a half-fit, the fraction of pairs that land in the same half-fit
module. Computed in closed form per half-fit as

    density_h(G) = sum_c C(|G ∩ module_c(h)|, 2) / C(|G assigned in h|, 2)

(no O(pairs) enumeration), averaged over the ensemble. Secondary, descriptive: median
best-match gene-set Jaccard across the ensemble.

No new model fits — pure partition arithmetic over artifacts already on disk.
"""

from __future__ import annotations

import argparse
from math import comb

import numpy as np
import pandas as pd
from scipy.stats import chi2, hypergeom, norm, spearmanr

from isograph_benchmark.paths import OUTPUT_DIRS, ensure_dir, rel, stage_out
from isograph_benchmark.real_data.stability import (
    COHORTS, SEED_BASE, _filter_expressed_transcripts, _split_indices,
)

# production output dir name <- split-half partition method tag
METHOD_DIRS = {"isograph": "isograph_vae", "wgcna": "wgcna_gene"}
# (cohort, region) -> production fit root
_STORE = OUTPUT_DIRS["modules"]

PROD_ROOTS = {
    ("brainseq", "caudate"): (*_STORE, "brainseq", "caudate", "_m"),
    ("brainseq", "hippocampus"): (*_STORE, "brainseq", "hippocampus", "_m"),
    ("brainseq", "dlpfc"): (*_STORE, "brainseq", "dlpfc", "_m"),
    ("gtex", "caudate_basal_ganglia"): (*_STORE, "gtex", "caudate_basal_ganglia", "_m"),
    ("gtex", "hippocampus"): (*_STORE, "gtex", "hippocampus", "_m"),
    ("gtex", "frontal_cortex_ba9"): (*_STORE, "gtex", "frontal_cortex_ba9", "_m"),
}

# cross-cohort caudate-matched pairs (BrainSEQ region, GTEx region) for Q3 replication
REGION_PAIRS = {
    "caudate": (("brainseq", "caudate"), ("gtex", "caudate_basal_ganglia")),
    "hippocampus": (("brainseq", "hippocampus"), ("gtex", "hippocampus")),
    "dlpfc_ba9": (("brainseq", "dlpfc"), ("gtex", "frontal_cortex_ba9")),
}


def _out_dir():
    return ensure_dir(stage_out("trust.stability", "module_trust"))


def _load_production_modules(cohort: str, region: str, method: str) -> dict[str, set]:
    root = PROD_ROOTS[(cohort, region)]
    path = rel(*root, METHOD_DIRS[method], "modules.parquet")
    if not path.exists():
        raise SystemExit(f"production modules missing: {path}")
    df = pd.read_parquet(path)
    df["gene_id"] = df["gene_id"].astype(str)
    return {m: set(g) for m, g in df.groupby("module_id")["gene_id"]}


def _load_halffit_maps(cohort: str, region: str, method: str) -> list[dict[str, str]]:
    """gene_id -> module label for each split-half partition of this method."""
    pdir = stage_out("trust.stability", "partitions")
    prefix = f"{method}__{cohort}__{region}__"
    maps = []
    for p in sorted(pdir.iterdir()):
        if not (p.name.startswith(prefix) and p.suffix == ".parquet"):
            continue
        df = pd.read_parquet(p)
        df["gene_id"] = df["gene_id"].astype(str)
        maps.append(dict(zip(df["gene_id"], df["module_id"].astype(str))))
    if not maps:
        raise SystemExit(f"no split-half partitions matching {prefix}* in {pdir}")
    return maps


def _density_in_half(genes, gene_map: dict[str, str]) -> float | None:
    """Closed-form within-set co-assignment density for one half-fit; None if <2 of the
    set's genes are assigned in this half (no pair to score)."""
    labels = [gene_map[g] for g in genes if g in gene_map]
    k = len(labels)
    if k < 2:
        return None
    counts = pd.Series(labels).value_counts().to_numpy()
    same_pairs = int(sum(comb(int(c), 2) for c in counts))
    return same_pairs / comb(k, 2)


def _mean_density(genes, halfmaps: list[dict[str, str]]) -> tuple[float, int]:
    vals = [d for d in (_density_in_half(genes, hm) for hm in halfmaps) if d is not None]
    if not vals:
        return float("nan"), 0
    return float(np.mean(vals)), len(vals)


def _best_match_jaccard(genes: set, half_modules: list[list[set]]) -> float:
    """Median over half-fits of the best gene-set Jaccard against any half-fit module."""
    best = []
    for mods in half_modules:
        j = 0.0
        for hm in mods:
            inter = len(genes & hm)
            if inter:
                j = max(j, inter / len(genes | hm))
        best.append(j)
    return float(np.median(best)) if best else float("nan")


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def stability(cohort: str, region: str, method: str, n_perm: int, seed: int,
              fdr: float) -> None:
    prod = _load_production_modules(cohort, region, method)
    halfmaps = _load_halffit_maps(cohort, region, method)
    half_modules = [
        [set(g for g, m in hm.items() if m == lbl) for lbl in set(hm.values())]
        for hm in halfmaps
    ]
    universe = np.array(sorted({g for hm in halfmaps for g in hm}))
    rng = np.random.default_rng(seed)
    print(f"[{cohort}/{region}/{method}] {len(prod)} production modules | "
          f"{len(halfmaps)} half-fits | universe {len(universe)} genes | "
          f"{n_perm} permutations", flush=True)

    rows = []
    for mid, genes in prod.items():
        obs, n_eff = _mean_density(genes, halfmaps)
        k = len(genes)
        # size-matched permutation null on the same statistic
        null = np.empty(n_perm)
        for b in range(n_perm):
            rand = set(universe[rng.choice(len(universe), size=k, replace=False)])
            null[b], _ = _mean_density(rand, halfmaps)
        valid = null[np.isfinite(null)]
        # one-sided p: how often does chance match/exceed the observed density
        p = (1 + int(np.sum(valid >= obs))) / (1 + len(valid)) if np.isfinite(obs) else 1.0
        rows.append({
            "cohort": cohort, "region": region, "method": method, "module_id": mid,
            "n_genes": k, "n_genes_assigned": n_eff,
            "coassign_density": obs, "null_mean": float(np.mean(valid)) if len(valid) else float("nan"),
            "perm_p": p, "best_match_jaccard": _best_match_jaccard(genes, half_modules),
        })

    out = pd.DataFrame(rows)
    out["fdr"] = _bh_fdr(out["perm_p"].to_numpy())
    out["trusted"] = out["fdr"] < fdr
    out = out.sort_values("coassign_density", ascending=False).reset_index(drop=True)
    path = _out_dir() / f"module_stability__{cohort}__{region}__{method}.parquet"
    out.to_parquet(path, index=False)

    n_trust = int(out["trusted"].sum())
    print(f"\n=== Q1 module stability: {n_trust} of {len(out)} modules trusted "
          f"(FDR<{fdr}) ===", flush=True)
    show = out[["module_id", "n_genes", "n_genes_assigned", "coassign_density",
                "null_mean", "perm_p", "fdr", "best_match_jaccard", "trusted"]]
    with pd.option_context("display.float_format", lambda x: f"{x:.3f}",
                           "display.max_rows", None):
        print(show.to_string(index=False), flush=True)
    print(f"\nwrote {path}", flush=True)


def _module_drivers(cohort: str, region: str, k: int, method: str = "isograph") -> dict[str, set]:
    """Top-k driver transcripts per production module (top |r| in the explain-module
    transcript_polarity_table, i.e. isoforms whose usage tracks the module).

    Only the IsoGraph backends produce a ``module_interpret`` tree; the classical WGCNA
    baseline has no transcript-level driver notion, so it yields an empty mapping and the
    driver-Jaccard column reads 0 rather than silently borrowing IsoGraph's drivers.
    """
    root = PROD_ROOTS[(cohort, region)]
    interp = rel(*root, METHOD_DIRS[method], "module_interpret")
    out: dict[str, set] = {}
    if not interp.exists():
        return out
    for mdir in interp.iterdir():
        tp = mdir / "transcript_polarity_table.parquet"
        if not (mdir.is_dir() and tp.exists()):
            continue
        df = pd.read_parquet(tp, columns=["transcript_id", "r"]).dropna(subset=["r"])
        if df.empty:
            continue
        top = df.reindex(df["r"].abs().sort_values(ascending=False).index).head(k)
        out[mdir.name] = set(top["transcript_id"].astype(str))
    return out


def _module_genes(cohort: str, region: str, method: str = "isograph") -> dict[str, set]:
    """gene set per production module (quantifier-robust matching key)."""
    return _load_production_modules(cohort, region, method)


AGE_MODELS = ("linear", "spline")


def _age_effects(cohort: str, region: str, method: str = "isograph",
                 model: str = "linear") -> pd.DataFrame:
    """Per-module Age effect and p-value, under one of the two fitted age models.

    ``linear`` is the published covariate-free Pearson r (``age_linear.parquet``).
    ``spline`` is the covariate-adjusted df=3 natural cubic spline the rest of the
    manuscript uses: the omnibus F-test p, with the direction taken from the projected
    late-vs-early contrast (effect at the 90th age percentile minus the 10th), so a
    non-monotone trajectory still gets a defined sign.
    """
    root = PROD_ROOTS[(cohort, region)]
    if model == "linear":
        df = pd.read_parquet(rel(*root, METHOD_DIRS[method], "age_linear.parquet"))
        df["module_id"] = df["module_id"].astype(str)
        return df.set_index("module_id")[["effect", "pvalue", "fdr"]]
    if model != "spline":
        raise SystemExit(f"unknown age model {model!r}; choose from {AGE_MODELS}")

    df = pd.read_parquet(rel(*root, METHOD_DIRS[method], "age_spline.parquet"))
    df["module_id"] = df["module_id"].astype(str)
    rows = []
    for mid, g in df.groupby("module_id"):
        g = g.sort_values("age_prob")
        rows.append({
            "module_id": mid,
            "effect": float(g["effect"].iloc[-1] - g["effect"].iloc[0]),
            "pvalue": float(g["pvalue_ftest"].iloc[0]),
            "fdr": float(g["fdr_ftest"].iloc[0]),
        })
    return pd.DataFrame(rows).set_index("module_id")[["effect", "pvalue", "fdr"]]


def _trusted_set(cohort: str, region: str, method: str) -> set:
    path = _out_dir() / f"module_stability__{cohort}__{region}__{method}.parquet"
    if not path.exists():
        raise SystemExit(f"run Q1 `stability` for {cohort}/{region}/{method} first ({path})")
    df = pd.read_parquet(path)
    return set(df.loc[df["trusted"], "module_id"].astype(str))


def _crosscohort_rows(pair: str, method: str, k: int, sig: float,
                      model: str = "linear") -> list[dict]:
    """One row per trusted discovery (BrainSEQ) module: its best gene-Jaccard match in the
    replication (GTEx) cohort, with both cohorts' module-Age effect/p and driver overlap.
    Shared by per-pair ``replication`` and pooled ``replication_pooled``."""
    (bc, br), (gc, gr) = REGION_PAIRS[pair]
    bs_trust = _trusted_set(bc, br, method)
    gt_trust = _trusted_set(gc, gr, method)
    # Every artifact below must come from `method`'s own fit.  Reading IsoGraph's modules
    # and Age effects while selecting module ids from WGCNA's trusted list looks plausible —
    # both backends label modules M000, M001, ... — but silently produces IsoGraph rows under
    # a WGCNA heading rather than raising.
    bs_genes, gt_genes = _module_genes(bc, br, method), _module_genes(gc, gr, method)
    bs_drv, gt_drv = _module_drivers(bc, br, k, method), _module_drivers(gc, gr, k, method)
    bs_age = _age_effects(bc, br, method, model)
    gt_age = _age_effects(gc, gr, method, model)
    print(f"[{pair}] {bc}/{br} ({len(bs_trust)} trusted) <-> {gc}/{gr} "
          f"({len(gt_trust)} trusted) | gene-set matching (driver Jaccard reported "
          f"separately; isoform drivers are quantifier-sensitive)", flush=True)

    rows = []
    for m in sorted(bs_trust):
        gm = bs_genes.get(m, set())
        # best-matching GTEx module by GENE-set Jaccard (genes are quantifier-robust;
        # top driver isoforms are not, so they cannot key the cross-cohort match)
        best, bestj = None, 0.0
        for g, gg in gt_genes.items():
            j = len(gm & gg) / len(gm | gg) if (gm | gg) else 0.0
            if j > bestj:
                best, bestj = g, j
        dm, dg = bs_drv.get(m, set()), gt_drv.get(best, set())
        drv_j = len(dm & dg) / len(dm | dg) if (best and (dm | dg)) else 0.0
        eb = bs_age.loc[m] if m in bs_age.index else None
        eg = gt_age.loc[best] if (best is not None and best in gt_age.index) else None
        sign_match = bool(eb is not None and eg is not None
                          and np.sign(eb["effect"]) == np.sign(eg["effect"]))
        both_sig = bool(eb is not None and eg is not None
                        and eb["pvalue"] < sig and eg["pvalue"] < sig)
        rows.append({
            "pair": pair, "bs_module": m, "gtex_match": best, "gene_jaccard": bestj,
            "driver_jaccard": drv_j, "match_trusted": best in gt_trust if best else False,
            "age_effect_bs": float(eb["effect"]) if eb is not None else np.nan,
            "age_p_bs": float(eb["pvalue"]) if eb is not None else np.nan,
            "age_effect_gtex": float(eg["effect"]) if eg is not None else np.nan,
            "age_p_gtex": float(eg["pvalue"]) if eg is not None else np.nan,
            "sign_match": sign_match, "both_sig": both_sig,
            "replicates": bool(sign_match and both_sig and bestj > 0),
        })
    return rows


def replication(pair: str, method: str, k: int, sig: float, model: str = "linear") -> None:
    """Q3 cross-cohort aging replication: match trusted modules across cohorts by driver-
    transcript overlap, then test Age-effect sign/significance concordance."""
    rows = _crosscohort_rows(pair, method, k, sig, model)
    out = pd.DataFrame(rows).sort_values("gene_jaccard", ascending=False).reset_index(drop=True)
    # The spline arm writes to its own prefix: trust_funnel_figure.R prefix-globs
    # "module_aging_replication__" and assemble_supp_tables.py reads those names exactly, so
    # reusing the prefix would silently mix the two age models in the figures and tables.
    stem = "module_aging_replication" if model == "linear" else f"module_aging_replication_{model}"
    path = _out_dir() / f"{stem}__{pair}__{method}.parquet"
    out.to_parquet(path, index=False)

    n_match = int((out["gene_jaccard"] > 0).sum())
    n_rep = int(out["replicates"].sum())
    print(f"\n=== Q3 cross-cohort aging replication ({pair}) ===", flush=True)
    print(f"trusted BrainSEQ modules: {len(out)} | gene-matched to GTEx: {n_match} | "
          f"aging replicates (sign+both p<{sig}): {n_rep}", flush=True)
    with pd.option_context("display.float_format", lambda x: f"{x:.3f}",
                           "display.max_rows", None):
        print(out[["bs_module", "gtex_match", "gene_jaccard", "driver_jaccard",
                   "match_trusted", "age_effect_bs", "age_effect_gtex", "sign_match",
                   "both_sig", "replicates"]].to_string(index=False), flush=True)
    print(f"\nwrote {path}", flush=True)


def _stouffer_z(disc_sign: np.ndarray, rep_eff: np.ndarray, rep_p: np.ndarray) -> float:
    """Directional Stouffer Z: for each discovery-significant module, convert the *replication*
    cohort's two-sided Age p-value into a one-sided p in the discovery effect direction, map to
    z = Phi^-1(1 - p_one), and combine as sum(z)/sqrt(K). Continuous per-module evidence, so it
    does not floor at the binomial sign-test's 0.5^K (the wall that sank the per-pair arm)."""
    if len(disc_sign) == 0:
        return float("nan")
    p2 = np.clip(rep_p, 1e-300, 1.0)
    same = np.sign(rep_eff) == disc_sign
    p_one = np.where(same, p2 / 2.0, 1.0 - p2 / 2.0)
    z = norm.isf(p_one)  # Phi^-1(1 - p_one)
    return float(z.sum() / np.sqrt(len(z)))


def replication_pooled(method: str, k: int, sig: float, min_jaccard: float,
                       n_perm: int, seed: int) -> None:
    """Pooled cross-cohort aging replication across ALL homologous region pairs.

    The per-pair ``replication`` arm is structurally underpowered: with only ~3 reproducible
    matched modules per pair, the binomial sign-test floor is 0.5^3 = 0.125, so even perfect
    concordance cannot reach p<0.05. Pooling the matched modules from every region pair into a
    single test breaks that floor in three complementary, increasingly assumption-light ways:

      (1) Stouffer directional meta-Z  -- HEADLINE. Among discovery (BrainSEQ) age-significant
          modules pooled across pairs, combine each module's *replication*-cohort directional z.
          Continuous evidence per module, so K=3 strong replications (p~0.01 each, in-direction)
          give Z~4 -> p~3e-5; the sign floor is gone. Reported with an analytic normal p AND a
          label-permutation p (shuffling which GTEx module pairs with each BrainSEQ module),
          which is robust to any within-region eigengene correlation that would violate the
          analytic independence assumption.
      (2) Sign-concordance permutation test -- shuffles the GTEx->BrainSEQ pairing, preserving
          each cohort's marginal sign distribution, so it controls for any global directional
          bias in aging (a bias that would inflate a naive binomial-vs-0.5 sign test).
      (3) Spearman magnitude concordance -- correlates discovery vs replication Age effect sizes
          across pooled pairs; uses magnitude, not just sign, with a permutation null.

    Pooling is valid because modules across region pairs are independent and gene sets within a
    pair are largely disjoint; the permutation nulls in (1)-(3) make no parametric independence
    claim. No new model fits -- reuses the trusted sets, production modules and Age effects."""
    frames, used, skipped = [], [], []
    for p in REGION_PAIRS:
        try:
            frames.append(pd.DataFrame(_crosscohort_rows(p, method, k, sig)))
            used.append(p)
        except SystemExit as e:  # missing Q1 trusted set / inputs for this pair -> skip
            print(f"[skip {p}] {e}", flush=True)
            skipped.append(p)
    if not frames:
        raise SystemExit("no region pairs have the required Q1 trusted sets + inputs yet")
    if len(used) < 2:
        print(f"\n!! WARNING: only {len(used)} region pair available ({used}); pooling needs "
              f">=2-3 pairs to clear the binomial sign floor. Run Q1 `stability` for the "
              f"missing pairs ({skipped}) for the powered result.", flush=True)
    allm = pd.concat(frames, ignore_index=True)
    rep = allm[(allm["gene_jaccard"] >= min_jaccard)
               & allm["age_effect_bs"].notna()
               & allm["age_effect_gtex"].notna()].reset_index(drop=True)
    path = _out_dir() / f"module_aging_replication_pooled__{method}.parquet"
    rep.to_parquet(path, index=False)

    n = len(rep)
    if n == 0:
        raise SystemExit(f"no reproducible matched pairs (gene Jaccard >= {min_jaccard}) "
                         f"pooled across {list(REGION_PAIRS)}")
    eb = rep["age_effect_bs"].to_numpy()
    eg = rep["age_effect_gtex"].to_numpy()
    pg = rep["age_p_gtex"].to_numpy()
    rng = np.random.default_rng(seed)

    # (1) Stouffer directional meta-Z among discovery (BrainSEQ) age-significant modules
    disc = rep[rep["age_p_bs"] < sig]
    d = np.sign(disc["age_effect_bs"].to_numpy())
    deg = disc["age_effect_gtex"].to_numpy()
    dpg = disc["age_p_gtex"].to_numpy()
    K = len(disc)
    Zs = _stouffer_z(d, deg, dpg)
    p_stouffer_norm = float(norm.sf(Zs)) if K else float("nan")
    if K:
        z_perm = np.array([_stouffer_z(d, deg[i := rng.permutation(K)], dpg[i])
                           for _ in range(n_perm)])
        p_stouffer_perm = float((1 + np.sum(z_perm >= Zs)) / (n_perm + 1))
    else:
        p_stouffer_perm = float("nan")

    # (2) sign concordance with label-permutation null (bias-robust)
    conc = np.sign(eb) == np.sign(eg)
    c_obs = float(conc.mean())
    c_perm = np.array([(np.sign(eb) == np.sign(rng.permutation(eg))).mean()
                       for _ in range(n_perm)])
    p_sign_perm = float((1 + np.sum(c_perm >= c_obs)) / (n_perm + 1))
    n_conc = int(conc.sum())
    p_binom = float(sum(comb(n, i) for i in range(n_conc, n + 1)) / 2 ** n)  # reference only

    # (3) magnitude concordance: Spearman of discovery vs replication Age effect
    rho = float(spearmanr(eb, eg).statistic)
    rho_perm = np.array([spearmanr(eb, rng.permutation(eg)).statistic for _ in range(n_perm)])
    p_rho = float((1 + np.sum(rho_perm >= rho)) / (n_perm + 1))

    stats = {
        "method": method, "n_pairs_reproducible": n, "min_jaccard": min_jaccard,
        "n_perm": n_perm, "sig": sig,
        "stouffer_K": K, "stouffer_Z": Zs,
        "stouffer_p_normal": p_stouffer_norm, "stouffer_p_perm": p_stouffer_perm,
        "sign_concordant": n_conc, "sign_frac": c_obs,
        "sign_p_perm": p_sign_perm, "sign_p_binom_vs_0.5": p_binom,
        "spearman_rho": rho, "spearman_p_perm": p_rho,
    }
    spath = _out_dir() / f"module_aging_replication_pooled__{method}__stats.json"
    import json
    spath.write_text(json.dumps(stats, indent=2))

    print(f"\n=== Q3 POOLED cross-cohort aging replication ({method}) ===", flush=True)
    print(f"region pairs pooled: {list(REGION_PAIRS)}", flush=True)
    print(f"reproducible matched pairs (gene J>={min_jaccard}): {n}", flush=True)
    print(f"\n(1) Stouffer directional meta-Z [HEADLINE] over K={K} discovery-significant "
          f"(BrainSEQ p<{sig}) modules:", flush=True)
    print(f"    Z={Zs:.3f} | analytic one-sided p={p_stouffer_norm:.3g} | "
          f"permutation p={p_stouffer_perm:.3g}", flush=True)
    print(f"(2) sign concordance: {n_conc}/{n} ({100*c_obs:.0f}%) | "
          f"permutation p={p_sign_perm:.3g} | (binomial-vs-0.5 ref p={p_binom:.3g})", flush=True)
    print(f"(3) magnitude (Spearman) discovery vs replication Age effect: rho={rho:.3f} | "
          f"permutation p={p_rho:.3g}", flush=True)
    print(f"\nwrote {path}\nwrote {spath}", flush=True)


def _composition_unique_genes(cohort: str, region: str) -> set:
    """Genes flagged DTU-without-DGE (switch significant, abundance not) -- invisible to
    an abundance method like WGCNA."""
    root = PROD_ROOTS[(cohort, region)]
    path = rel(*root, "isograph_vae", "composition_unique", "genes.parquet")
    if not path.exists():
        return set()
    df = pd.read_parquet(path)
    return set(df.loc[df["category"] == "composition_unique", "gene_id"].astype(str))


def _wgcna_age_genes(cohort: str, region: str, fdr: float) -> set:
    """Union of genes in age-associated WGCNA modules (the abundance-based age signal)."""
    root = PROD_ROOTS[(cohort, region)]
    mods = pd.read_parquet(rel(*root, "wgcna_gene", "modules.parquet"))
    mods["gene_id"] = mods["gene_id"].astype(str)
    age = pd.read_parquet(rel(*root, "wgcna_gene", "age_linear.parquet"))
    sig = set(age.loc[age["fdr"] < fdr, "module_id"].astype(str))
    return set(mods.loc[mods["module_id"].astype(str).isin(sig), "gene_id"])


def _structure_flags(cohort: str, region: str) -> pd.DataFrame:
    root = PROD_ROOTS[(cohort, region)]
    path = rel(*root, "isograph_vae", "module_interpret", "structure_annotations.parquet")
    df = pd.read_parquet(path)
    df["transcript_id"] = df["transcript_id"].astype(str)
    return df.set_index("transcript_id")


def complementarity(cohort: str, region: str, method: str, k: int, fdr: float) -> None:
    """Q4: for trusted modules, quantify the biology WGCNA cannot see -- DTU-without-DGE
    gene fraction, (non-)overlap with age-associated WGCNA modules, and the splicing-level
    structural switches of their driver transcripts."""
    stab = pd.read_parquet(
        _out_dir() / f"module_stability__{cohort}__{region}__{method}.parquet")
    trusted = stab.loc[stab["trusted"], "module_id"].astype(str).tolist()
    genes = _module_genes(cohort, region, method)
    drivers = _module_drivers(cohort, region, k, method)
    dtu = _composition_unique_genes(cohort, region)
    wgcna_age = _wgcna_age_genes(cohort, region, fdr)
    age = _age_effects(cohort, region, method)
    struct = _structure_flags(cohort, region)
    sflags = ["cds_changed", "utr_changed", "biotype_switch", "coding_status_change"]
    print(f"[{cohort}/{region}/{method}] {len(trusted)} trusted modules | "
          f"{len(dtu)} DTU-without-DGE genes | {len(wgcna_age)} genes in age-WGCNA modules",
          flush=True)

    rows = []
    for m in trusted:
        g = genes.get(m, set())
        drv = [t for t in drivers.get(m, set()) if t in struct.index]
        srow = {f"drv_{f}": (float(struct.loc[drv, f].mean()) if drv else float("nan"))
                for f in sflags}
        eb = age.loc[m] if m in age.index else None
        rows.append({
            "cohort": cohort, "region": region, "module_id": m, "n_genes": len(g),
            "n_dtu_without_dge": len(g & dtu),
            "frac_dtu_without_dge": len(g & dtu) / len(g) if g else 0.0,
            "frac_in_wgcna_age_modules": len(g & wgcna_age) / len(g) if g else 0.0,
            "age_effect": float(eb["effect"]) if eb is not None else np.nan,
            "age_fdr": float(eb["fdr"]) if eb is not None else np.nan,
            "age_sig": bool(eb is not None and eb["fdr"] < fdr),
            **srow,
        })
    out = pd.DataFrame(rows).sort_values(
        ["age_sig", "n_dtu_without_dge"], ascending=False).reset_index(drop=True)
    path = _out_dir() / f"module_complementarity__{cohort}__{region}__{method}.parquet"
    out.to_parquet(path, index=False)

    n_dtu_mod = int((out["n_dtu_without_dge"] > 0).sum())
    headline = out[out["age_sig"]]
    print(f"\n=== Q4 WGCNA-complementarity ({cohort}/{region}) ===", flush=True)
    print(f"trusted modules carrying DTU-without-DGE genes: {n_dtu_mod}/{len(out)} | "
          f"total DTU-without-DGE genes in trusted modules: {int(out['n_dtu_without_dge'].sum())}",
          flush=True)
    print(f"trusted & age-associated modules (the headline set): {len(headline)}", flush=True)
    with pd.option_context("display.float_format", lambda x: f"{x:.3f}",
                           "display.max_rows", None):
        cols = ["module_id", "n_genes", "n_dtu_without_dge", "frac_dtu_without_dge",
                "frac_in_wgcna_age_modules", "age_effect", "age_sig",
                "drv_cds_changed", "drv_biotype_switch"]
        print(out[cols].head(15).to_string(index=False), flush=True)
    print(f"\nwrote {path}", flush=True)
    print("NOTE: per-module GO enrichment not yet wired (raw GO files in "
          "inputs/go_annotations/); add as a follow-on.", flush=True)


def _meta_dir():
    return ensure_dir(stage_out("trust.stability", "modules_meta"))


def meta(cohort: str, region: str, method: str, k: int) -> None:
    """Step 0 (post-hoc, no re-fit): for each split-half partition, reconstruct the
    module eigengenes + Age effect (faithfully, via the same feature construction the fit
    used) and the per-module top-k driver transcripts (from switch-axis loadings). The
    module assignment is taken from the saved partition; everything else is a deterministic
    function of the bundle + the seeded split, so no VAE re-run is needed.
    """
    from isograph.features.channels import gene_feature_channels, make_feature_scores
    from isograph.features.residualize import build_design_matrix, residualize_rows
    from isograph.features.switch import gene_switch_loadings
    from isograph.io.artifacts import load_dataset_bundle
    from isograph.models.base import compute_trait_associations

    spec = COHORTS[cohort]
    age_col = spec["age_col"]
    bundle = load_dataset_bundle(rel(*spec["bundle_root"], region))
    sample_table = bundle.sample_table.reset_index(drop=True)
    tc = bundle.matrices["transcript_counts"]
    tt = bundle.feature_tables["transcript"]
    if spec["filter_transcripts"]:
        tc, tt = _filter_expressed_transcripts(tc, tt)
    else:
        tc = np.asarray(tc)
    del bundle
    n = sample_table.shape[0]

    pdir = stage_out("trust.stability", "partitions")
    prefix = f"{method}__{cohort}__{region}__"
    parts = sorted(p for p in pdir.iterdir()
                   if p.name.startswith(prefix) and p.suffix == ".parquet")
    print(f"[{cohort}/{region}/{method}] reconstructing meta for {len(parts)} half-fits",
          flush=True)

    out_rows = []
    load_rows = []  # per-(seed,half,module,gene) switch-axis importance for the Q2 reframe
    for p in parts:
        # parse seed + half from <method>__<cohort>__<region>__seed<k>__<half>.parquet
        stem = p.stem[len(prefix):]
        seedtok, half = stem.split("__")
        seed = int(seedtok.replace("seed", ""))
        a_idx, b_idx = _split_indices(n, SEED_BASE + seed)
        idx = a_idx if half == "A" else b_idx
        st = sample_table.iloc[idx].reset_index(drop=True)
        part = pd.read_parquet(p)[["gene_id", "module_id"]]
        part["gene_id"] = part["gene_id"].astype(str)

        switch_matrix, feature_info = gene_feature_channels(tc[:, idx], tt)
        if switch_matrix.size:
            design = build_design_matrix(st, spec["covariates"])
            switch_matrix = residualize_rows(switch_matrix, design)
        fs = make_feature_scores(switch_matrix, feature_info, st)
        assoc, _ = compute_trait_associations(part, fs, st, [age_col])
        age = assoc[assoc["trait"] == age_col].set_index("module_id")

        load = gene_switch_loadings(tc[:, idx], tt)
        load["abs"] = load["loading"].abs()
        gene2mod = dict(zip(part["gene_id"], part["module_id"].astype(str)))
        load["module_id"] = load["gene_id"].map(gene2mod)
        for mid, grp in load.dropna(subset=["module_id"]).groupby("module_id"):
            grp = grp.sort_values("abs", ascending=False)
            drivers = grp["transcript_id"].head(k).tolist()
            a = age.loc[mid] if mid in age.index else None
            out_rows.append({
                "cohort": cohort, "region": region, "method": method,
                "seed": seed, "half": half, "module_id": str(mid),
                "age_effect": float(a["effect"]) if a is not None else np.nan,
                "age_pvalue": float(a["pvalue"]) if a is not None else np.nan,
                "drivers": ";".join(drivers),
            })
            # per-gene switch-axis importance (max |loading| over the gene's transcripts).
            # |loading| because the SVD axis sign is arbitrary per half-fit (not comparable).
            gimp = grp.groupby("gene_id")["abs"].max()
            load_rows.append(pd.DataFrame({
                "seed": seed, "half": half, "module_id": str(mid),
                "gene_id": gimp.index.astype(str), "importance": gimp.to_numpy(),
            }))
        print(f"  {p.name}: {part['module_id'].nunique()} modules", flush=True)

    out = pd.DataFrame(out_rows)
    path = _meta_dir() / f"modules_meta__{cohort}__{region}__{method}.parquet"
    out.to_parquet(path, index=False)
    lpath = _meta_dir() / f"modules_meta_loadings__{cohort}__{region}__{method}.parquet"
    pd.concat(load_rows, ignore_index=True).to_parquet(lpath, index=False)
    n_sig = int((out["age_pvalue"] < 0.05).sum())
    print(f"\n=== Step-0 meta: {len(out)} (module,half) rows | "
          f"age p<0.05 in {n_sig} ({100*n_sig/max(len(out),1):.0f}%) ===", flush=True)
    print(f"wrote {path}\nwrote {lpath}", flush=True)


def _halffit_module_genes(cohort: str, region: str, method: str) -> dict:
    """(seed, half) -> {module_id: set(gene_id)} for each split-half partition. Module
    labels match the meta parquet (both read module_id straight from the partition)."""
    pdir = stage_out("trust.stability", "partitions")
    prefix = f"{method}__{cohort}__{region}__"
    out: dict = {}
    for p in sorted(pdir.iterdir()):
        if not (p.name.startswith(prefix) and p.suffix == ".parquet"):
            continue
        seedtok, half = p.stem[len(prefix):].split("__")
        seed = int(seedtok.replace("seed", ""))
        df = pd.read_parquet(p)[["gene_id", "module_id"]]
        df["gene_id"] = df["gene_id"].astype(str)
        out[(seed, half)] = {str(m): set(g) for m, g in df.groupby("module_id")["gene_id"]}
    return out


def _transcript_gene_map(cohort: str, region: str) -> dict[str, str]:
    """transcript_id -> gene_id from the bundle's transcripts.parquet (cheap; no VAE load)."""
    spec = COHORTS[cohort]
    tt = pd.read_parquet(rel(*spec["bundle_root"], region, "transcripts.parquet"),
                         columns=["transcript_id", "gene_id"])
    return dict(zip(tt["transcript_id"].astype(str), tt["gene_id"].astype(str)))


def within(cohort: str, region: str, method: str, sig: float, min_jaccard: float) -> None:
    """Within-cohort Q3 (aging-sign concordance) + Q2 (driver reproducibility). For each
    seed the two independent halves are matched module-to-module by best gene-set Jaccard;
    a pair is *reproducible* when that Jaccard >= ``min_jaccard``. Q3: among reproducible
    pairs, do the two halves agree on the sign of the module-Age effect (overall, and among
    pairs where both halves are age-nominal)? Q2 (driver reproducibility) is reported at TWO
    levels: (i) the raw cross-half top-k driver-*transcript* Jaccard (descriptive; weak by
    construction -- isoform identity is quantifier-sensitive), and (ii) the REFRAMED
    driver-*gene* test, which is the level that actually reproduces. (ii) is made non-circular
    by conditioning on the matched pair's *shared* genes: among genes present in both halves'
    modules, are the same shared genes flagged as drivers more than chance (per-pair
    hypergeometric), and pooled across reproducible pairs (Fisher)? Pure post-hoc arithmetic
    over the meta parquet + partitions + the bundle transcript->gene map."""
    mpath = _meta_dir() / f"modules_meta__{cohort}__{region}__{method}.parquet"
    if not mpath.exists():
        raise SystemExit(f"run `meta` for {cohort}/{region}/{method} first ({mpath})")
    meta_df = pd.read_parquet(mpath)
    meta_df["module_id"] = meta_df["module_id"].astype(str)
    # (seed, half, module_id) -> (age_effect, age_pvalue, driver set)
    rec: dict = {}
    for r in meta_df.itertuples(index=False):
        drv = set(str(r.drivers).split(";")) if isinstance(r.drivers, str) and r.drivers else set()
        rec[(int(r.seed), r.half, str(r.module_id))] = (r.age_effect, r.age_pvalue, drv)

    genes = _halffit_module_genes(cohort, region, method)
    t2g = _transcript_gene_map(cohort, region)
    # per-(seed,half,module) gene->switch-importance vectors for the loading-correlation reframe
    lpath = _meta_dir() / f"modules_meta_loadings__{cohort}__{region}__{method}.parquet"
    imp: dict = {}
    if lpath.exists():
        ldf = pd.read_parquet(lpath)
        ldf["module_id"] = ldf["module_id"].astype(str)
        ldf["gene_id"] = ldf["gene_id"].astype(str)
        for (s, h, m), grp in ldf.groupby(["seed", "half", "module_id"]):
            imp[(int(s), h, m)] = grp.set_index("gene_id")["importance"]
    seeds = sorted({s for (s, _) in genes})
    print(f"[{cohort}/{region}/{method}] within-cohort A<->B matching over {len(seeds)} seeds "
          f"| reproducible = gene Jaccard >= {min_jaccard}", flush=True)

    rows = []
    for seed in seeds:
        ga, gb = genes.get((seed, "A")), genes.get((seed, "B"))
        if not ga or not gb:
            continue
        for ma, gset_a in ga.items():
            best, bestj, gset_best = None, 0.0, set()
            for mb, gset_b in gb.items():
                j = len(gset_a & gset_b) / len(gset_a | gset_b) if (gset_a | gset_b) else 0.0
                if j > bestj:
                    best, bestj, gset_best = mb, j, gset_b
            if best is None:
                continue
            ea, pa, da = rec.get((seed, "A", ma), (np.nan, np.nan, set()))
            eb, pb, db = rec.get((seed, "B", best), (np.nan, np.nan, set()))
            drv_j = len(da & db) / len(da | db) if (da | db) else np.nan
            # REFRAMED Q2: driver *genes* (top-k transcripts mapped to genes), tested against
            # the matched pair's SHARED-gene background so the matching itself can't inflate it.
            dga = {t2g[t] for t in da if t in t2g}
            dgb = {t2g[t] for t in db if t in t2g}
            drv_gene_j = len(dga & dgb) / len(dga | dgb) if (dga | dgb) else np.nan
            shared = gset_a & gset_best
            DA, DB = dga & shared, dgb & shared  # driver genes drawn from the shared background
            x, N, n1, n2 = len(DA & DB), len(shared), len(DA), len(DB)
            # P(overlap >= x) for n2 draws from N with n1 successes; 1.0 if any margin empty
            hyp_p = float(hypergeom.sf(x - 1, N, n1, n2)) if (N and n1 and n2) else np.nan
            # PRIMARY reframe: do the SHARED genes rank as drivers consistently across halves?
            # Spearman of per-gene switch-importance over the shared genes (continuous, no
            # top-k cutoff, conditioned on shared genes so it is not inflated by the matching).
            load_rho, load_p, n_shared_imp = np.nan, np.nan, 0
            ia, ib = imp.get((seed, "A", ma)), imp.get((seed, "B", best))
            if ia is not None and ib is not None:
                common = sorted(shared & set(ia.index) & set(ib.index))
                n_shared_imp = len(common)
                if n_shared_imp >= 5:
                    sr = spearmanr(ia.loc[common].to_numpy(), ib.loc[common].to_numpy())
                    load_rho, load_p = float(sr.statistic), float(sr.pvalue)
            rows.append({
                "cohort": cohort, "region": region, "method": method, "seed": seed,
                "module_a": ma, "module_b": best, "gene_jaccard": bestj,
                "reproducible": bestj >= min_jaccard,
                "age_effect_a": ea, "age_effect_b": eb,
                "age_p_a": pa, "age_p_b": pb,
                "sign_concordant": bool(np.sign(ea) == np.sign(eb))
                                   if (pd.notna(ea) and pd.notna(eb)) else False,
                "both_age_sig": bool(pd.notna(pa) and pd.notna(pb) and pa < sig and pb < sig),
                "driver_jaccard": drv_j, "driver_gene_jaccard": drv_gene_j,
                "n_shared_genes": N, "driver_gene_overlap": x,
                "driver_gene_hyp_p": hyp_p,
                "driver_load_rho": load_rho, "driver_load_p": load_p,
                "n_shared_importance": n_shared_imp,
            })
    out = pd.DataFrame(rows).sort_values(["seed", "gene_jaccard"],
                                         ascending=[True, False]).reset_index(drop=True)
    path = _out_dir() / f"within_cohort__{cohort}__{region}__{method}.parquet"
    out.to_parquet(path, index=False)

    repro = out[out["reproducible"]]
    both = repro[repro["both_age_sig"]]
    n_conc_repro = int(repro["sign_concordant"].sum())
    n_conc_both = int(both["sign_concordant"].sum())
    # binomial sign-test p (one-sided, H0: p=0.5) on the both-age-sig pairs
    from math import comb as _comb
    nb = len(both)
    sign_p = (sum(_comb(nb, i) for i in range(n_conc_both, nb + 1)) / 2 ** nb
              if nb else float("nan"))
    print(f"\n=== within-cohort Q3 aging concordance ({cohort}/{region}) ===", flush=True)
    print(f"A<->B matched pairs: {len(out)} | reproducible (J>={min_jaccard}): {len(repro)} "
          f"| median gene Jaccard (repro): {repro['gene_jaccard'].median():.3f}", flush=True)
    print(f"sign-concordant among reproducible: {n_conc_repro}/{len(repro)} "
          f"({100*n_conc_repro/max(len(repro),1):.0f}%)", flush=True)
    print(f"sign-concordant among reproducible & both-halves-age-nominal: "
          f"{n_conc_both}/{nb} ({100*n_conc_both/max(nb,1):.0f}%) | sign-test p={sign_p:.3g}",
          flush=True)
    print(f"\n=== Q2 driver reproducibility ({cohort}/{region}) ===", flush=True)
    dj = repro["driver_jaccard"].dropna()
    print(f"(i) raw top-k driver-TRANSCRIPT Jaccard [descriptive, weak by construction]: "
          f"median={dj.median():.3f} mean={dj.mean():.3f} | pairs with any shared driver: "
          f"{int((dj > 0).sum())}/{len(dj)}", flush=True)
    # (ii) reframed driver-GENE level, conditioned on the shared-gene background
    dgj = repro["driver_gene_jaccard"].dropna()
    hp = repro["driver_gene_hyp_p"].dropna()
    n_pair_sig = int((hp < sig).sum())
    # Fisher pooled across reproducible pairs (independent matched pairs): -2 sum ln p ~ chi2_2k
    hp_pos = hp[hp > 0]
    fisher_stat = float(-2.0 * np.log(hp_pos).sum()) if len(hp_pos) else float("nan")
    fisher_p = float(chi2.sf(fisher_stat, 2 * len(hp_pos))) if len(hp_pos) else float("nan")
    # (iii) PRIMARY reframe: shared-gene switch-importance rank correlation across halves
    rr = repro.dropna(subset=["driver_load_rho"])
    if len(rr):
        rho_med = rr["driver_load_rho"].median()
        n_rho_pos = int((rr["driver_load_rho"] > 0).sum())
        rp_pos = rr.loc[rr["driver_load_p"] > 0, "driver_load_p"]
        # one-sided (positive-concordance) p per pair, then Fisher-pool across pairs
        rp_one = np.where(rr["driver_load_rho"] > 0, rr["driver_load_p"] / 2,
                          1 - rr["driver_load_p"] / 2)
        rp_one = rp_one[rp_one > 0]
        f_stat = float(-2.0 * np.log(rp_one).sum()) if len(rp_one) else float("nan")
        f_p = float(chi2.sf(f_stat, 2 * len(rp_one))) if len(rp_one) else float("nan")
        f_p_str = "<1e-300 (underflow)" if f_p == 0 else f"{f_p:.3g}"
        rho_lo, rho_hi = rr["driver_load_rho"].min(), rr["driver_load_rho"].max()
        print(f"(iii) PRIMARY reframe -- shared-gene switch-importance rank concordance "
              f"(Spearman over each matched pair's shared genes, >=5 genes):", flush=True)
        print(f"      EFFECT SIZE (headline): median rho={rho_med:.3f} [{rho_lo:.3f}-{rho_hi:.3f}]"
              f" | rho>0 in {n_rho_pos}/{len(rr)} pairs", flush=True)
        print(f"      significance: pooled one-sided Fisher p={f_p_str} (note: per-pair p is "
              f"driven by the 100s of shared genes, so rho is the meaningful quantity)", flush=True)
        print(f"      caveat: reproducible switch-importance may be partly gene-intrinsic "
              f"(stable isoform structure), not only module/aging-specific", flush=True)
    else:
        miss = "loadings parquet absent" if not lpath.exists() else "no pair has >=5 shared genes"
        print(f"(iii) PRIMARY reframe -- shared-gene importance correlation: not computed "
              f"({miss})", flush=True)
    print(f"(ii) driver-GENE set test (k-limited): driver-gene Jaccard median={dgj.median():.3f} "
          f"| shared-background hypergeometric POOLED Fisher p={fisher_p:.3g} "
          f"({n_pair_sig}/{len(hp)} pairs p<{sig})", flush=True)
    print(f"\nwrote {path}", flush=True)


def replication_model_contrast(method: str = "isograph") -> None:
    """Decompose the linear-vs-spline gap in the Q3 concordance count.

    Panel C of the trust funnel reports the covariate-free linear (Pearson) arm, while the
    covariate-adjusted spline is the age model used elsewhere in the manuscript.  The two
    give very different counts, and a reader is entitled to know *which component* moves.
    Splitting the count into its two conjuncts answers that: `sign_match` is whether the two
    cohorts agree on direction, `both_sig` is whether both clear the p cutoff.  Reporting
    only the final count would leave it ambiguous whether the models disagree about the
    direction of aging or merely about how much of it survives covariate adjustment.

    Reads the per-pair tables both `replication --model {linear,spline}` runs already wrote;
    it re-fits nothing, so it cannot drift from the numbers those tables carry.
    """
    out_dir = _out_dir()
    rows = []
    for model in AGE_MODELS:
        tag = "" if model == "linear" else "_spline"
        frames = []
        for path in sorted(out_dir.iterdir()):
            name = path.name
            if not name.startswith(f"module_aging_replication{tag}__"):
                continue
            if "_pooled__" in name or not name.endswith(f"__{method}.parquet"):
                continue
            # linear's prefix is a strict prefix of nothing else, but spell the guard out:
            # "module_aging_replication__" must not swallow "..._spline__".
            if model == "linear" and "_spline__" in name:
                continue
            frames.append(pd.read_parquet(path))
        if not frames:
            continue
        t = pd.concat(frames, ignore_index=True)
        rows.append({
            "method": method,
            "model": model,
            "n_pairs": int(len(t)),
            "sign_match": int(t["sign_match"].sum()),
            "both_sig": int(t["both_sig"].sum()),
            "concordant": int(t["replicates"].sum()),
            "discovery_sig": int((t["age_p_bs"] < 0.05).sum()),
            "replication_sig": int((t["age_p_gtex"] < 0.05).sum()),
        })
    if not rows:
        raise SystemExit(f"no replication tables found for method={method!r}; run "
                         "`replication --model linear` and `--model spline` first")
    out = pd.DataFrame(rows)
    path = out_dir / f"replication_model_contrast__{method}.parquet"
    out.to_parquet(path, index=False, compression="zstd")
    print(out.to_string(index=False), flush=True)
    if len(out) == 2:
        lin = out[out["model"] == "linear"].iloc[0]
        spl = out[out["model"] == "spline"].iloc[0]
        print(f"\nsign_match: {lin.sign_match} -> {spl.sign_match} "
              f"(delta {spl.sign_match - lin.sign_match:+d})", flush=True)
        print(f"both_sig:   {lin.both_sig} -> {spl.both_sig} "
              f"(delta {spl.both_sig - lin.both_sig:+d})", flush=True)
        print(f"discovery-cohort modules at p<0.05:   {lin.discovery_sig} -> "
              f"{spl.discovery_sig}", flush=True)
        print(f"replication-cohort modules at p<0.05: {lin.replication_sig} -> "
              f"{spl.replication_sig}", flush=True)
    print(f"\nwrote {path}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    st = sub.add_parser("stability", help="Q1: trust gate from co-assignment density")
    st.add_argument("--cohort", default="brainseq")
    st.add_argument("--region", default="caudate")
    st.add_argument("--method", default="isograph", choices=list(METHOD_DIRS))
    st.add_argument("--n-perm", type=int, default=1000)
    st.add_argument("--seed", type=int, default=13)
    st.add_argument("--fdr", type=float, default=0.05)
    rp = sub.add_parser("replication", help="Q3: cross-cohort aging replication of trusted modules")
    rp.add_argument("--pair", default="caudate", choices=list(REGION_PAIRS))
    rp.add_argument("--method", default="isograph", choices=list(METHOD_DIRS))
    rp.add_argument("--k", type=int, default=5, help="top-k driver transcripts for matching")
    rp.add_argument("--sig", type=float, default=0.05, help="Age p-value cutoff for 'both_sig'")
    rp.add_argument("--model", default="linear", choices=list(AGE_MODELS),
                    help="age model: linear = published covariate-free Pearson; "
                         "spline = covariate-adjusted df=3 spline F-test (written to "
                         "module_aging_replication_spline__*)")
    pl = sub.add_parser("replication-pooled",
                        help="Q3 pooled: cross-cohort aging replication pooled over all region "
                             "pairs (Stouffer meta-Z + permutation tests; breaks the n=3 floor)")
    pl.add_argument("--method", default="isograph", choices=list(METHOD_DIRS))
    pl.add_argument("--k", type=int, default=5, help="top-k driver transcripts for driver Jaccard")
    pl.add_argument("--sig", type=float, default=0.05, help="Age p cutoff for discovery-significant")
    pl.add_argument("--min-jaccard", type=float, default=0.25,
                    help="gene-set Jaccard for a reproducible cross-cohort match")
    pl.add_argument("--n-perm", type=int, default=10000)
    pl.add_argument("--seed", type=int, default=13)
    cm = sub.add_parser("complementarity", help="Q4: DTU-without-DGE / WGCNA-complementarity")
    cm.add_argument("--cohort", default="brainseq")
    cm.add_argument("--region", default="caudate")
    cm.add_argument("--method", default="isograph", choices=list(METHOD_DIRS))
    cm.add_argument("--k", type=int, default=5, help="top-k driver transcripts for structure")
    cm.add_argument("--fdr", type=float, default=0.05)
    mt = sub.add_parser("meta", help="Step 0: per-half module eigengene Age effect + drivers")
    mt.add_argument("--cohort", default="brainseq")
    mt.add_argument("--region", default="caudate")
    mt.add_argument("--method", default="isograph", choices=list(METHOD_DIRS))
    mt.add_argument("--k", type=int, default=5)
    wn = sub.add_parser("within", help="within-cohort Q3 aging concordance + Q2 driver reproducibility")
    wn.add_argument("--cohort", default="brainseq")
    wn.add_argument("--region", default="caudate")
    wn.add_argument("--method", default="isograph", choices=list(METHOD_DIRS))
    wn.add_argument("--sig", type=float, default=0.05)
    wn.add_argument("--min-jaccard", type=float, default=0.25)
    mc = sub.add_parser("replication-model-contrast",
                        help="decompose the linear-vs-spline Q3 count into sign agreement "
                             "vs both-cohort significance")
    mc.add_argument("--method", default="isograph", choices=list(METHOD_DIRS))
    args = ap.parse_args()
    if args.cmd == "stability":
        stability(args.cohort, args.region, args.method, args.n_perm, args.seed, args.fdr)
    elif args.cmd == "replication":
        replication(args.pair, args.method, args.k, args.sig, args.model)
    elif args.cmd == "replication-pooled":
        replication_pooled(args.method, args.k, args.sig, args.min_jaccard,
                           args.n_perm, args.seed)
    elif args.cmd == "complementarity":
        complementarity(args.cohort, args.region, args.method, args.k, args.fdr)
    elif args.cmd == "meta":
        meta(args.cohort, args.region, args.method, args.k)
    elif args.cmd == "within":
        within(args.cohort, args.region, args.method, args.sig, args.min_jaccard)
    elif args.cmd == "replication-model-contrast":
        replication_model_contrast(args.method)


if __name__ == "__main__":
    main()
