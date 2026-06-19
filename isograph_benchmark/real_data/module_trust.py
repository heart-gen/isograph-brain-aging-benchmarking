"""Per-module trust funnel for real-data IsoGraph modules (see
``real_data/stability/MODULE_TRUST_PLAN.md``).

Q1 (this module, ``stability`` command): which *production* modules are stable enough to
trust? A production module's gene set is scored for how tightly its genes stay co-clustered
across the split-half ensemble (`real_data/stability/_m/partitions/`), relative to a
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

from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.stability import (
    COHORTS, SEED_BASE, _filter_expressed_transcripts, _split_indices,
)

# production output dir name <- split-half partition method tag
METHOD_DIRS = {"isograph": "isograph_vae", "wgcna": "wgcna_gene"}
# (cohort, region) -> production fit root
PROD_ROOTS = {
    ("brainseq", "caudate"): ("real_data", "brainseq", "caudate", "_m"),
    ("brainseq", "hippocampus"): ("real_data", "brainseq", "hippocampus", "_m"),
    ("brainseq", "dlpfc"): ("real_data", "brainseq", "dlpfc", "_m"),
    ("gtex", "caudate_basal_ganglia"): ("real_data", "gtex", "caudate_basal_ganglia", "_m"),
    ("gtex", "hippocampus"): ("real_data", "gtex", "hippocampus", "_m"),
    ("gtex", "frontal_cortex_ba9"): ("real_data", "gtex", "frontal_cortex_ba9", "_m"),
}

# cross-cohort caudate-matched pairs (BrainSEQ region, GTEx region) for Q3 replication
REGION_PAIRS = {
    "caudate": (("brainseq", "caudate"), ("gtex", "caudate_basal_ganglia")),
    "hippocampus": (("brainseq", "hippocampus"), ("gtex", "hippocampus")),
    "dlpfc_ba9": (("brainseq", "dlpfc"), ("gtex", "frontal_cortex_ba9")),
}


def _out_dir():
    return ensure_dir(rel("real_data", "stability", "_m", "module_trust"))


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
    pdir = rel("real_data", "stability", "_m", "partitions")
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


def _module_drivers(cohort: str, region: str, k: int) -> dict[str, set]:
    """Top-k driver transcripts per production module (top |r| in the explain-module
    transcript_polarity_table, i.e. isoforms whose usage tracks the module)."""
    root = PROD_ROOTS[(cohort, region)]
    interp = rel(*root, "isograph_vae", "module_interpret")
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


def _module_genes(cohort: str, region: str) -> dict[str, set]:
    """gene set per production module (quantifier-robust matching key)."""
    return _load_production_modules(cohort, region, "isograph")


def _age_effects(cohort: str, region: str) -> pd.DataFrame:
    root = PROD_ROOTS[(cohort, region)]
    df = pd.read_parquet(rel(*root, "isograph_vae", "age_linear.parquet"))
    df["module_id"] = df["module_id"].astype(str)
    return df.set_index("module_id")[["effect", "pvalue", "fdr"]]


def _trusted_set(cohort: str, region: str, method: str) -> set:
    path = _out_dir() / f"module_stability__{cohort}__{region}__{method}.parquet"
    if not path.exists():
        raise SystemExit(f"run Q1 `stability` for {cohort}/{region}/{method} first ({path})")
    df = pd.read_parquet(path)
    return set(df.loc[df["trusted"], "module_id"].astype(str))


def replication(pair: str, method: str, k: int, sig: float) -> None:
    """Q3 cross-cohort aging replication: match trusted modules across cohorts by driver-
    transcript overlap, then test Age-effect sign/significance concordance."""
    (bc, br), (gc, gr) = REGION_PAIRS[pair]
    bs_trust = _trusted_set(bc, br, method)
    gt_trust = _trusted_set(gc, gr, method)
    bs_genes, gt_genes = _module_genes(bc, br), _module_genes(gc, gr)
    bs_drv, gt_drv = _module_drivers(bc, br, k), _module_drivers(gc, gr, k)
    bs_age, gt_age = _age_effects(bc, br), _age_effects(gc, gr)
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
    out = pd.DataFrame(rows).sort_values("gene_jaccard", ascending=False).reset_index(drop=True)
    path = _out_dir() / f"module_aging_replication__{pair}__{method}.parquet"
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
    genes = _module_genes(cohort, region)
    drivers = _module_drivers(cohort, region, k)
    dtu = _composition_unique_genes(cohort, region)
    wgcna_age = _wgcna_age_genes(cohort, region, fdr)
    age = _age_effects(cohort, region)
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
    return ensure_dir(rel("real_data", "stability", "_m", "modules_meta"))


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

    pdir = rel("real_data", "stability", "_m", "partitions")
    prefix = f"{method}__{cohort}__{region}__"
    parts = sorted(p for p in pdir.iterdir()
                   if p.name.startswith(prefix) and p.suffix == ".parquet")
    print(f"[{cohort}/{region}/{method}] reconstructing meta for {len(parts)} half-fits",
          flush=True)

    out_rows = []
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
            drivers = grp.sort_values("abs", ascending=False)["transcript_id"].head(k).tolist()
            a = age.loc[mid] if mid in age.index else None
            out_rows.append({
                "cohort": cohort, "region": region, "method": method,
                "seed": seed, "half": half, "module_id": str(mid),
                "age_effect": float(a["effect"]) if a is not None else np.nan,
                "age_pvalue": float(a["pvalue"]) if a is not None else np.nan,
                "drivers": ";".join(drivers),
            })
        print(f"  {p.name}: {part['module_id'].nunique()} modules", flush=True)

    out = pd.DataFrame(out_rows)
    path = _meta_dir() / f"modules_meta__{cohort}__{region}__{method}.parquet"
    out.to_parquet(path, index=False)
    n_sig = int((out["age_pvalue"] < 0.05).sum())
    print(f"\n=== Step-0 meta: {len(out)} (module,half) rows | "
          f"age p<0.05 in {n_sig} ({100*n_sig/max(len(out),1):.0f}%) ===", flush=True)
    print(f"wrote {path}", flush=True)


def _halffit_module_genes(cohort: str, region: str, method: str) -> dict:
    """(seed, half) -> {module_id: set(gene_id)} for each split-half partition. Module
    labels match the meta parquet (both read module_id straight from the partition)."""
    pdir = rel("real_data", "stability", "_m", "partitions")
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


def within(cohort: str, region: str, method: str, sig: float, min_jaccard: float) -> None:
    """Within-cohort Q3 (aging-sign concordance) + Q2 (driver reproducibility). For each
    seed the two independent halves are matched module-to-module by best gene-set Jaccard;
    a pair is *reproducible* when that Jaccard >= ``min_jaccard``. Q3: among reproducible
    pairs, do the two halves agree on the sign of the module-Age effect (overall, and among
    pairs where both halves are age-nominal)? Q2: cross-half top-k driver-transcript Jaccard
    on the same matched pairs. Pure post-hoc arithmetic over the meta parquet + partitions."""
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
    seeds = sorted({s for (s, _) in genes})
    print(f"[{cohort}/{region}/{method}] within-cohort A<->B matching over {len(seeds)} seeds "
          f"| reproducible = gene Jaccard >= {min_jaccard}", flush=True)

    rows = []
    for seed in seeds:
        ga, gb = genes.get((seed, "A")), genes.get((seed, "B"))
        if not ga or not gb:
            continue
        for ma, gset_a in ga.items():
            best, bestj = None, 0.0
            for mb, gset_b in gb.items():
                j = len(gset_a & gset_b) / len(gset_a | gset_b) if (gset_a | gset_b) else 0.0
                if j > bestj:
                    best, bestj = mb, j
            if best is None:
                continue
            ea, pa, da = rec.get((seed, "A", ma), (np.nan, np.nan, set()))
            eb, pb, db = rec.get((seed, "B", best), (np.nan, np.nan, set()))
            drv_j = len(da & db) / len(da | db) if (da | db) else np.nan
            rows.append({
                "cohort": cohort, "region": region, "method": method, "seed": seed,
                "module_a": ma, "module_b": best, "gene_jaccard": bestj,
                "reproducible": bestj >= min_jaccard,
                "age_effect_a": ea, "age_effect_b": eb,
                "age_p_a": pa, "age_p_b": pb,
                "sign_concordant": bool(np.sign(ea) == np.sign(eb))
                                   if (pd.notna(ea) and pd.notna(eb)) else False,
                "both_age_sig": bool(pd.notna(pa) and pd.notna(pb) and pa < sig and pb < sig),
                "driver_jaccard": drv_j,
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
    print(f"\n=== Q2 driver-transcript reproducibility ({cohort}/{region}) ===", flush=True)
    dj = repro["driver_jaccard"].dropna()
    print(f"cross-half top-k driver Jaccard over reproducible pairs: "
          f"median={dj.median():.3f} mean={dj.mean():.3f} | pairs with any shared driver: "
          f"{int((dj > 0).sum())}/{len(dj)}", flush=True)
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
    st.add_argument("--seed", type=int, default=0)
    st.add_argument("--fdr", type=float, default=0.05)
    rp = sub.add_parser("replication", help="Q3: cross-cohort aging replication of trusted modules")
    rp.add_argument("--pair", default="caudate", choices=list(REGION_PAIRS))
    rp.add_argument("--method", default="isograph", choices=list(METHOD_DIRS))
    rp.add_argument("--k", type=int, default=5, help="top-k driver transcripts for matching")
    rp.add_argument("--sig", type=float, default=0.05, help="Age p-value cutoff for 'both_sig'")
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
    args = ap.parse_args()
    if args.cmd == "stability":
        stability(args.cohort, args.region, args.method, args.n_perm, args.seed, args.fdr)
    elif args.cmd == "replication":
        replication(args.pair, args.method, args.k, args.sig)
    elif args.cmd == "complementarity":
        complementarity(args.cohort, args.region, args.method, args.k, args.fdr)
    elif args.cmd == "meta":
        meta(args.cohort, args.region, args.method, args.k)
    elif args.cmd == "within":
        within(args.cohort, args.region, args.method, args.sig, args.min_jaccard)


if __name__ == "__main__":
    main()
