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

# production output dir name <- split-half partition method tag
METHOD_DIRS = {"isograph": "isograph_vae", "wgcna": "wgcna_gene"}
# (cohort, region) -> production fit root
PROD_ROOTS = {
    ("brainseq", "caudate"): ("real_data", "brainseq", "caudate", "_m"),
    ("brainseq", "hippocampus"): ("real_data", "brainseq", "hippocampus", "_m"),
    ("brainseq", "dlpfc"): ("real_data", "brainseq", "dlpfc", "_m"),
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
    args = ap.parse_args()
    if args.cmd == "stability":
        stability(args.cohort, args.region, args.method, args.n_perm, args.seed, args.fdr)


if __name__ == "__main__":
    main()
