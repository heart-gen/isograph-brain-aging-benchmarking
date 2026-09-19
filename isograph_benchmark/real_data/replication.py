"""Cross-cohort replication: BrainSEQ vs GTEx for matched brain regions.

Two independent postmortem cohorts (BrainSEQ Salmon counts; GTEx v11 RSEM counts)
are profiled with the same pipeline. For the three brain regions present in both,
this module asks how reproducible the discovered structure is, with two
complementary analyses run per method (``isograph_vae`` and ``wgcna_gene``):

1. **Module preservation** — for each module in one cohort, the best-matching
   module in the other cohort (highest Jaccard over the shared expressed genes),
   in both directions. Significance is assessed against a label-permutation null
   that preserves module sizes, giving an empirical p-value and z-score per
   module. A module counts as *preserved* when its observed best Jaccard exceeds
   the null at p < 0.05.

2. **Age-effect concordance** — for module pairs matched by best Jaccard, whether
   the eigengene–age association (the ``effect`` in ``age_linear.parquet``, a
   Pearson correlation with age) agrees in sign across cohorts (binomial test vs
   0.5) and in magnitude (Spearman across matched pairs). This tests whether the
   aging signal — not just the gene grouping — reproduces.

Gene IDs are versioned Ensembl IDs that match directly across cohorts (same
GENCODE annotation), so no ID harmonization is needed; analyses are restricted to
genes expressed (bundle-filtered) in both cohorts for a region.

Usage::

    python -m isograph_benchmark.real_data.replication                 # both methods
    python -m isograph_benchmark.real_data.replication --methods isograph_vae

Outputs (under ``04_module_trust/_m/replication/``):
    <method>_module_match.parquet   per-source-module best match + age concordance
    replication_summary.parquet     per method × region × direction summary
    replication_summary.json        compact headline summary
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

from isograph_benchmark.paths import ensure_dir, region_store, stage_out

# (brainseq_region, gtex_region, label) for the three matched brain regions.
REGION_PAIRS = [
    ("caudate", "caudate_basal_ganglia", "caudate"),
    ("hippocampus", "hippocampus", "hippocampus"),
    ("dlpfc", "frontal_cortex_ba9", "dlpfc_ba9"),
]

DEFAULT_METHODS = ["isograph_vae", "wgcna_gene"]
N_PERMUTATIONS = 1000
MIN_OVERLAP_FOR_MATCH = 5      # min shared genes for a match to seed age concordance
PRESERVED_P = 0.05
SEED = 13


def _modules_path(cohort: str, region: str, method: str):
    return region_store(cohort, region, method, "modules.parquet")


def _age_linear_path(cohort: str, region: str, method: str):
    return region_store(cohort, region, method, "age_linear.parquet")


def _load_modules(cohort: str, region: str, method: str) -> pd.DataFrame | None:
    p = _modules_path(cohort, region, method)
    if not p.exists():
        return None
    df = pd.read_parquet(p)
    return df[["gene_id", "module_id"]].dropna()


def _masks(modules: pd.DataFrame, gene_index: dict[str, int]) -> tuple[np.ndarray, list[str]]:
    """Boolean (n_modules × n_shared_genes) membership matrix over shared genes."""
    mod_ids = sorted(modules["module_id"].unique())
    mask = np.zeros((len(mod_ids), len(gene_index)), dtype=bool)
    pos = {m: i for i, m in enumerate(mod_ids)}
    for gene, mod in zip(modules["gene_id"], modules["module_id"]):
        gi = gene_index.get(gene)
        if gi is not None:
            mask[pos[mod], gi] = True
    keep = mask.any(axis=1)
    return mask[keep], [m for m, k in zip(mod_ids, keep) if k]


def _best_jaccard(src_mask: np.ndarray, tgt_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Best-match Jaccard of each source module against all target modules.

    Returns (best_jaccard, best_target_idx, best_overlap) per source module.
    """
    inter = src_mask.astype(np.int32) @ tgt_mask.astype(np.int32).T   # (nS, nT)
    src_sz = src_mask.sum(1)[:, None]
    tgt_sz = tgt_mask.sum(1)[None, :]
    union = src_sz + tgt_sz - inter
    jacc = np.where(union > 0, inter / union, 0.0)
    best_idx = jacc.argmax(1)
    rows = np.arange(jacc.shape[0])
    return jacc[rows, best_idx], best_idx, inter[rows, best_idx]


def _permutation_null(src_mask: np.ndarray, tgt_mask: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Null distribution of each source module's best Jaccard under permuted target labels.

    Returns (N_PERMUTATIONS × n_source) array. Shuffling target columns preserves
    target module sizes while destroying gene-level correspondence.
    """
    n_genes = tgt_mask.shape[1]
    null = np.empty((N_PERMUTATIONS, src_mask.shape[0]), dtype=float)
    for k in range(N_PERMUTATIONS):
        perm = rng.permutation(n_genes)
        best, _, _ = _best_jaccard(src_mask, tgt_mask[:, perm])
        null[k] = best
    return null


def _direction(
    src_mod: pd.DataFrame, tgt_mod: pd.DataFrame,
    src_age: pd.DataFrame | None, tgt_age: pd.DataFrame | None,
    shared_genes: list[str], method: str, region: str, direction: str,
    rng: np.random.Generator,
) -> pd.DataFrame:
    gene_index = {g: i for i, g in enumerate(shared_genes)}
    src_mask, src_ids = _masks(src_mod, gene_index)
    tgt_mask, tgt_ids = _masks(tgt_mod, gene_index)
    if src_mask.shape[0] == 0 or tgt_mask.shape[0] == 0:
        return pd.DataFrame()

    best_j, best_idx, best_ov = _best_jaccard(src_mask, tgt_mask)
    null = _permutation_null(src_mask, tgt_mask, rng)
    null_mean = null.mean(0)
    null_sd = null.std(0)
    perm_p = (1.0 + (null >= best_j[None, :]).sum(0)) / (N_PERMUTATIONS + 1)
    perm_z = np.where(null_sd > 0, (best_j - null_mean) / null_sd, np.nan)

    src_eff = (src_age.set_index("module_id")["effect"].to_dict() if src_age is not None else {})
    src_fdr = (src_age.set_index("module_id")["fdr"].to_dict() if src_age is not None else {})
    tgt_eff = (tgt_age.set_index("module_id")["effect"].to_dict() if tgt_age is not None else {})

    rows = []
    for i, sid in enumerate(src_ids):
        tid = tgt_ids[best_idx[i]]
        se = src_eff.get(sid, np.nan)
        te = tgt_eff.get(tid, np.nan)
        sign_conc = (
            bool(np.sign(se) == np.sign(te))
            if np.isfinite(se) and np.isfinite(te) and se != 0 and te != 0
            else None
        )
        rows.append({
            "method": method, "region": region, "direction": direction,
            "source_module": sid, "n_source_genes": int(src_mask[i].sum()),
            "best_target_module": tid, "n_overlap": int(best_ov[i]),
            "jaccard": float(best_j[i]),
            "null_mean_jaccard": float(null_mean[i]), "perm_p": float(perm_p[i]),
            "perm_z": float(perm_z[i]) if np.isfinite(perm_z[i]) else np.nan,
            "preserved": bool(perm_p[i] < PRESERVED_P),
            "source_age_effect": float(se) if np.isfinite(se) else np.nan,
            "source_age_fdr": float(src_fdr.get(sid, np.nan)),
            "target_age_effect": float(te) if np.isfinite(te) else np.nan,
            "age_sign_concordant": sign_conc,
        })
    return pd.DataFrame(rows)


def run_method(method: str) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    parts = []
    for bs_region, gtex_region, label in REGION_PAIRS:
        bs_mod = _load_modules("brainseq", bs_region, method)
        gt_mod = _load_modules("gtex", gtex_region, method)
        if bs_mod is None or gt_mod is None:
            missing = "brainseq" if bs_mod is None else "gtex"
            print(f"  [{method}] {label}: missing {missing} modules — skipping")
            continue
        shared = sorted(set(bs_mod["gene_id"]) & set(gt_mod["gene_id"]))
        if len(shared) < MIN_OVERLAP_FOR_MATCH:
            print(f"  [{method}] {label}: only {len(shared)} shared genes — skipping")
            continue
        bs_age = _read_optional(_age_linear_path("brainseq", bs_region, method))
        gt_age = _read_optional(_age_linear_path("gtex", gtex_region, method))
        print(f"  [{method}] {label}: {len(shared)} shared genes | "
              f"BrainSEQ {bs_mod['module_id'].nunique()} mods, GTEx {gt_mod['module_id'].nunique()} mods")
        parts.append(_direction(bs_mod, gt_mod, bs_age, gt_age, shared, method, label, "brainseq_to_gtex", rng))
        parts.append(_direction(gt_mod, bs_mod, gt_age, bs_age, shared, method, label, "gtex_to_brainseq", rng))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def _read_optional(path):
    return pd.read_parquet(path) if path.exists() else None


def summarize(match: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (method, region, direction), g in match.groupby(["method", "region", "direction"]):
        matched = g[g["n_overlap"] >= MIN_OVERLAP_FOR_MATCH]
        conc = matched["age_sign_concordant"].dropna()
        sign_rate = float(conc.mean()) if len(conc) else np.nan
        sign_p = (binomtest(int(conc.sum()), len(conc)).pvalue if len(conc) else np.nan)
        pair = matched.dropna(subset=["source_age_effect", "target_age_effect"])
        if len(pair) >= 4:
            rho, rho_p = spearmanr(pair["source_age_effect"], pair["target_age_effect"])
        else:
            rho, rho_p = np.nan, np.nan
        rows.append({
            "method": method, "region": region, "direction": direction,
            "n_modules": int(len(g)),
            "mean_best_jaccard": float(g["jaccard"].mean()),
            "median_best_jaccard": float(g["jaccard"].median()),
            "n_preserved": int(g["preserved"].sum()),
            "frac_preserved": float(g["preserved"].mean()),
            "n_age_matched": int(len(conc)),
            "age_sign_concordance": sign_rate,
            "age_sign_binom_p": float(sign_p) if np.isfinite(sign_p) else np.nan,
            "age_effect_spearman": float(rho) if np.isfinite(rho) else np.nan,
            "age_effect_spearman_p": float(rho_p) if np.isfinite(rho_p) else np.nan,
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="BrainSEQ vs GTEx cross-cohort replication.")
    ap.add_argument("--methods", nargs="+", default=DEFAULT_METHODS,
                    help=f"Methods to run (default: {DEFAULT_METHODS}).")
    args = ap.parse_args()

    out_dir = stage_out("trust.replication")
    ensure_dir(out_dir)

    all_match = []
    for method in args.methods:
        print(f"Running replication for {method} ...")
        m = run_method(method)
        if m.empty:
            print(f"  {method}: no results (modules not available?)")
            continue
        m.to_parquet(out_dir / f"{method}_module_match.parquet", index=False, compression="zstd")
        all_match.append(m)

    if not all_match:
        raise SystemExit("No replication results produced — are the module fits present?")

    match = pd.concat(all_match, ignore_index=True)
    summary = summarize(match)
    summary.to_parquet(out_dir / "replication_summary.parquet", index=False, compression="zstd")

    headline = {
        "region_pairs": [{"brainseq": b, "gtex": g, "label": l} for b, g, l in REGION_PAIRS],
        "n_permutations": N_PERMUTATIONS,
        "methods": sorted(match["method"].unique().tolist()),
        "by_method_region": [],
    }
    for (method, region), g in summary.groupby(["method", "region"]):
        headline["by_method_region"].append({
            "method": method, "region": region,
            "mean_best_jaccard": round(float(g["mean_best_jaccard"].mean()), 3),
            "frac_preserved": round(float(g["frac_preserved"].mean()), 3),
            "age_sign_concordance": round(float(g["age_sign_concordance"].mean()), 3),
            "age_effect_spearman": round(float(g["age_effect_spearman"].mean()), 3),
        })
    (out_dir / "replication_summary.json").write_text(json.dumps(headline, indent=2))

    print(f"\nWrote {out_dir}/replication_summary.parquet ({len(summary)} rows) + .json")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
