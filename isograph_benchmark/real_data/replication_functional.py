"""Functional preservation of the cross-cohort matched modules.

The BrainSEQ<->GTEx matched pairs share very few genes (median gene Jaccard ~0.02), which
on its own reads as a failure to reproduce.  But gene-level overlap is the strictest
possible criterion: two cohorts can recover the *same biological program* through partly
different genes.  This asks whether they do, at three levels of organisation:

* **pathway** — Jaccard of the two modules' enriched GO:BP term sets,
* **cell type** — correlation of their cell-type marker enrichment profiles,
* **structure** — correlation of their switch coding-consequence class-rate profiles.

The control is the piece ``replication_go.py`` lacks: **size-matched** random pairs.  Module
similarity rises with module size for all three measures, and the matched partner of a large
module is itself usually large, so an unmatched-random null would be beaten by size alone.
Each BrainSEQ module is therefore re-paired with a random GTEx module drawn from the same
gene-count decile.

Pairs come from ``module_trust``'s cross-cohort table, so this speaks to exactly the same
modules as the replication permutation test (``replication_permutation.py``) rather than
``replication.py``'s separate matching.

Decision rule, pre-registered: only if the age-concordant pairs are *significantly* more
functionally similar than size-matched random pairs may the manuscript describe the cohorts
as recovering related biological programs.  Otherwise the language stays at "matched modules
with concordant age effects", and "replicated programs" is not used.

Not every measure is computable for every method.  ``celltype_composition`` is only run for
the IsoGraph artifact tree, and ``switch_consequence`` skips any region with no
phenotype-significant switch genes -- which is the case for three of the six regions here,
so ``structure_r`` has no data on either side of some pairs.  A missing input yields NaN
rather than an error, so the report carries an explicit **Data coverage** section naming
which inputs were absent; an "n/a" row must never be read as a null result.

Outputs (under ``03_module_trust/_m/replication/``):
    functional_preservation__{method}__{model}.parquet    one row per matched pair
    functional_preservation__{method}__{model}__stats.json
    FUNCTIONAL_PRESERVATION__{method}__{model}.md
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.paths import ensure_dir, region_store, rel, stage_out
from isograph_benchmark.real_data.module_trust import (
    METHOD_DIRS, PROD_ROOTS, REGION_PAIRS, _out_dir,
)

SEED = 13
N_PERMUTATIONS = 1000
N_BOOT = 2000
N_DECILES = 10

# module_trust method tag -> module_enrichment file prefix
_ENRICH_PREFIX = {"isograph": "isograph", "wgcna": "wgcna"}

_CONSEQUENCE_CLASSES = (
    "first_exon_changed", "last_exon_changed", "internal_exon_difference",
    "cds_changed", "utr_changed", "biotype_switch", "coding_status_change",
    "coding_consequence", "nmd_switch",
)

MEASURES = ("go_jaccard", "celltype_r", "structure_r")


# --------------------------------------------------------------------------- #
# Per-module profiles
# --------------------------------------------------------------------------- #
def _gene_stem(s: pd.Series) -> pd.Series:
    """Ensembl gene id without its version suffix."""
    return s.astype(str).str.split(".", n=1).str[0]


def _go_sets(cohort: str, region: str, method: str) -> dict[str, set]:
    """module_id -> set of enriched GO:BP term ids (empty set when none passed)."""
    base = region_store(cohort, region, "module_enrichment")
    prefix = _ENRICH_PREFIX[method]
    go_path, mod_path = base / f"{prefix}_module_go.parquet", base / f"{prefix}_modules.parquet"
    if not (go_path.exists() and mod_path.exists()):
        return {}
    sets: dict[str, set] = {
        str(m): set() for m in pd.read_parquet(mod_path)["module_id"].unique()
    }
    go = pd.read_parquet(go_path)
    for mid, grp in go.groupby("module_id"):
        sets[str(mid)] = set(grp["term_id"])
    return sets


def _celltype_profiles(cohort: str, region: str, method: str) -> dict[str, pd.Series]:
    """module_id -> log2 marker enrichment per cell type.

    A +0.5 offset keeps the ratio defined for modules with no observed markers of a type,
    which is common for the smaller modules.
    """
    path = rel(*PROD_ROOTS[(cohort, region)], METHOD_DIRS[method],
               "celltype_composition", "marker_enrichment.parquet")
    if not path.exists():
        return {}
    df = pd.read_parquet(path)
    df["value"] = np.log2((df["observed"] + 0.5) / (df["expected"] + 0.5))
    wide = df.pivot_table(index="module_id", columns="cell_type", values="value")
    return {str(m): row for m, row in wide.iterrows()}


def _structure_profiles(cohort: str, region: str, method: str) -> dict[str, pd.Series]:
    """module_id -> rate of each switch coding-consequence class over the module's genes."""
    root = PROD_ROOTS[(cohort, region)]
    pc = rel(*root, METHOD_DIRS[method], "switch_consequence", "pair_consequence.parquet")
    mods = rel(*root, METHOD_DIRS[method], "modules.parquet")
    if not (pc.exists() and mods.exists()):
        return {}
    pairs = pd.read_parquet(pc)
    modules = pd.read_parquet(mods)
    # switch_consequence strips the Ensembl version suffix ("ENSG…") while modules.parquet
    # keeps it ("ENSG….10"), so the join has to be on the stem or it silently matches nothing.
    gene_to_module = dict(zip(_gene_stem(modules["gene_id"]),
                              modules["module_id"].astype(str)))
    pairs["module_id"] = _gene_stem(pairs["gene"]).map(gene_to_module)
    pairs = pairs.dropna(subset=["module_id"])
    cols = [c for c in _CONSEQUENCE_CLASSES if c in pairs.columns]
    if pairs.empty or not cols:
        return {}
    rates = pairs.groupby("module_id")[cols].mean()
    return {str(m): row for m, row in rates.iterrows()}


def _probe(cohort: str, region: str, method: str) -> dict[str, bool]:
    """Which of the three profile inputs exist for this (cohort, region, method).

    Absent inputs propagate to NaN measures rather than raising, so record them: the
    difference between "no signal" and "no data" is not recoverable from the numbers.
    """
    base = region_store(cohort, region, "module_enrichment")
    prefix = _ENRICH_PREFIX[method]
    root = PROD_ROOTS[(cohort, region)]
    return {
        "go_jaccard": (base / f"{prefix}_module_go.parquet").exists()
        and (base / f"{prefix}_modules.parquet").exists(),
        "celltype_r": rel(*root, METHOD_DIRS[method], "celltype_composition",
                          "marker_enrichment.parquet").exists(),
        "structure_r": rel(*root, METHOD_DIRS[method], "switch_consequence",
                           "pair_consequence.parquet").exists(),
    }


def _module_sizes(cohort: str, region: str, method: str) -> pd.Series:
    df = pd.read_parquet(rel(*PROD_ROOTS[(cohort, region)], METHOD_DIRS[method],
                             "modules.parquet"))
    return df.groupby(df["module_id"].astype(str)).size()


# --------------------------------------------------------------------------- #
# Similarity measures
# --------------------------------------------------------------------------- #
def _jaccard(a: set, b: set) -> float:
    u = a | b
    return len(a & b) / len(u) if u else np.nan


def _profile_r(a: pd.Series | None, b: pd.Series | None) -> float:
    """Pearson r over the cell types / consequence classes both profiles define."""
    if a is None or b is None:
        return np.nan
    shared = [c for c in a.index if c in b.index]
    if len(shared) < 3:
        return np.nan
    x = pd.to_numeric(a[shared], errors="coerce").to_numpy(float)
    y = pd.to_numeric(b[shared], errors="coerce").to_numpy(float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3 or np.std(x[mask]) == 0 or np.std(y[mask]) == 0:
        return np.nan
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def _similarity(bs_module: str, gt_module: str | None, profiles: dict) -> dict[str, float]:
    if gt_module is None or (isinstance(gt_module, float) and not np.isfinite(gt_module)):
        return dict.fromkeys(MEASURES, np.nan)
    return {
        "go_jaccard": _jaccard(profiles["bs_go"].get(bs_module, set()),
                               profiles["gt_go"].get(gt_module, set())),
        "celltype_r": _profile_r(profiles["bs_ct"].get(bs_module),
                                 profiles["gt_ct"].get(gt_module)),
        "structure_r": _profile_r(profiles["bs_st"].get(bs_module),
                                  profiles["gt_st"].get(gt_module)),
    }


# --------------------------------------------------------------------------- #
# Size-matched null
# --------------------------------------------------------------------------- #
def _decile_pools(sizes: pd.Series) -> dict[int, list[str]]:
    """Gene-count decile -> the GTEx modules in it (the draw pool for a matched module)."""
    if sizes.empty:
        return {}
    ranks = sizes.rank(method="average", pct=True)
    bins = np.clip((ranks * N_DECILES).astype(int), 0, N_DECILES - 1)
    pools: dict[int, list[str]] = {}
    for mid, b in zip(sizes.index, bins):
        pools.setdefault(int(b), []).append(str(mid))
    return pools


def _decile_of(size: float, sizes: pd.Series) -> int:
    if sizes.empty:
        return 0
    pct = float((sizes <= size).mean())
    return int(np.clip(int(pct * N_DECILES), 0, N_DECILES - 1))


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def _matched_table(method: str, model: str) -> pd.DataFrame:
    stem = "module_aging_replication" if model == "linear" else f"module_aging_replication_{model}"
    frames = []
    for pair in REGION_PAIRS:
        path = _out_dir() / f"{stem}__{pair}__{method}.parquet"
        if not path.exists():
            raise SystemExit(
                f"missing {path}; run 03_module_trust/_h/09.module_trust_replication.sh"
                + (f" --model {model}" if model != "linear" else ""))
        frames.append(pd.read_parquet(path))
    return pd.concat(frames, ignore_index=True)


def run(method: str, model: str, n_perm: int, seed: int, n_boot: int) -> None:
    matched = _matched_table(method, model)
    rng = np.random.default_rng(seed)

    rows: list[dict] = []
    null_by_pair: dict[str, np.ndarray] = {}
    coverage: dict[str, dict] = {}

    for pair, sub in matched.groupby("pair"):
        (bc, br), (gc, gr) = REGION_PAIRS[pair]
        coverage[pair] = {
            f"{bc}/{br}": _probe(bc, br, method),
            f"{gc}/{gr}": _probe(gc, gr, method),
        }
        profiles = {
            "bs_go": _go_sets(bc, br, method), "gt_go": _go_sets(gc, gr, method),
            "bs_ct": _celltype_profiles(bc, br, method),
            "gt_ct": _celltype_profiles(gc, gr, method),
            "bs_st": _structure_profiles(bc, br, method),
            "gt_st": _structure_profiles(gc, gr, method),
        }
        gt_sizes = _module_sizes(gc, gr, method)
        bs_sizes = _module_sizes(bc, br, method)
        pools = _decile_pools(gt_sizes)

        obs_rows = []
        for r in sub.itertuples():
            gt = None if pd.isna(r.gtex_match) else str(r.gtex_match)
            sim = _similarity(str(r.bs_module), gt, profiles)
            obs_rows.append({
                "method": method, "model": model, "pair": pair,
                "bs_module": str(r.bs_module), "gtex_match": gt,
                "bs_n_genes": int(bs_sizes.get(str(r.bs_module), 0)),
                "gtex_n_genes": int(gt_sizes.get(gt, 0)) if gt else 0,
                "gene_jaccard": float(r.gene_jaccard),
                "replicates": bool(r.replicates), "both_sig": bool(r.both_sig),
                **sim,
            })
        rows.extend(obs_rows)

        # Size-matched random re-pairing: each BrainSEQ module keeps its own decile.
        draws = {m: np.full(n_perm, np.nan) for m in MEASURES}
        candidates = [
            (str(r["bs_module"]),
             pools.get(_decile_of(gt_sizes.get(str(r["gtex_match"]), np.nan)
                                  if r["gtex_match"] else np.nan, gt_sizes), []))
            for r in obs_rows
        ]
        for b in range(n_perm):
            vals = {m: [] for m in MEASURES}
            for bs_module, pool in candidates:
                if not pool:
                    continue
                pick = pool[int(rng.integers(0, len(pool)))]
                sim = _similarity(bs_module, pick, profiles)
                for m in MEASURES:
                    if np.isfinite(sim[m]):
                        vals[m].append(sim[m])
            for m in MEASURES:
                if vals[m]:
                    draws[m][b] = float(np.mean(vals[m]))
        null_by_pair[pair] = draws

    out = pd.DataFrame(rows)
    outdir = ensure_dir(stage_out("trust.replication"))
    # method AND model in every filename: the array runs both methods concurrently and the
    # linear/spline runs are separate results, so a fixed name silently loses one of them.
    stem = f"{method}__{model}"
    out.to_parquet(outdir / f"functional_preservation__{stem}.parquet", index=False)

    payload = _summarise(out, null_by_pair, method, model, n_perm, seed, n_boot)
    payload["coverage"] = coverage
    (outdir / f"functional_preservation__{stem}__stats.json").write_text(
        json.dumps(payload, indent=2))
    _write_report(out, payload, outdir, stem)
    print(f"[{method}/{model}] {len(out)} matched pairs -> {outdir}")
    for m in MEASURES:
        s = payload["measures"][m]
        print(f"  {m}: matched mean {s['mean_matched']} vs size-matched null "
              f"{s['null_mean']} (p_emp={s['p_emp']}); concordant vs discordant "
              f"MWU p={s['mwu_p_concordant_vs_not']}")


def _mean(v) -> float | None:
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(v.mean()) if v.size else None


def _boot_diff_ci(a, b, n_boot: int, seed: int) -> tuple[float | None, float | None]:
    a = np.asarray(a, float)[np.isfinite(a)]
    b = np.asarray(b, float)[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return None, None
    rng = np.random.default_rng(seed)
    diffs = np.empty(n_boot)
    for i in range(n_boot):
        diffs[i] = (a[rng.integers(0, a.size, a.size)].mean()
                    - b[rng.integers(0, b.size, b.size)].mean())
    lo, hi = np.quantile(diffs, [0.025, 0.975])
    return float(lo), float(hi)


def _summarise(out: pd.DataFrame, null_by_pair: dict, method: str, model: str,
               n_perm: int, seed: int, n_boot: int) -> dict:
    payload = {
        "method": method, "model": model, "n_pairs": int(len(out)),
        "n_concordant": int(out["replicates"].sum()),
        "n_permutations": n_perm, "n_boot": n_boot, "seed": seed,
        "median_gene_jaccard": float(out["gene_jaccard"].median()),
        "measures": {},
    }
    # Pool the per-pair size-matched nulls by averaging the three region pairs' draw means,
    # so the null statistic has the same shape as the observed overall mean.
    for m in MEASURES:
        stacks = [d[m] for d in null_by_pair.values() if m in d and np.isfinite(d[m]).any()]
        null = (np.nanmean(np.vstack(stacks), axis=0) if stacks else np.empty(0))
        null = null[np.isfinite(null)]
        obs = _mean(out[m])
        conc = out.loc[out["replicates"], m]
        disc = out.loc[~out["replicates"], m]
        mwu_p = None
        if conc.notna().sum() >= 3 and disc.notna().sum() >= 3:
            mwu_p = float(stats.mannwhitneyu(conc.dropna(), disc.dropna(),
                                             alternative="greater").pvalue)
        lo, hi = _boot_diff_ci(conc, disc, n_boot, seed)
        payload["measures"][m] = {
            "n_finite": int(np.isfinite(pd.to_numeric(out[m], errors="coerce")).sum()),
            "mean_matched": obs,
            "mean_concordant": _mean(conc),
            "mean_discordant": _mean(disc),
            "null_mean": float(null.mean()) if null.size else None,
            "null_sd": float(null.std(ddof=1)) if null.size > 1 else None,
            "p_emp": (float((1 + int((null >= obs).sum())) / (null.size + 1))
                      if null.size and obs is not None else None),
            "mwu_p_concordant_vs_not": mwu_p,
            "diff_ci_low": lo, "diff_ci_high": hi,
        }
    return payload


def _write_report(out: pd.DataFrame, payload: dict, outdir, stem: str) -> None:
    m = payload["measures"]
    lines = [
        "# Functional preservation of cross-cohort matched modules", "",
        f"Method `{payload['method']}`, age model `{payload['model']}`, "
        f"{payload['n_pairs']} matched BrainSEQ<->GTEx module pairs "
        f"({payload['n_concordant']} with concordant age effects). Median gene Jaccard "
        f"{payload['median_gene_jaccard']:.4f} — the gene-level overlap these measures are "
        f"asked to look past.", "",
        "The null re-pairs each BrainSEQ module with a random GTEx module **from the same "
        "gene-count decile**, because all three similarity measures increase with module "
        "size. `p_emp` is the size-matched permutation p for the overall matched mean; the "
        "Mann-Whitney column asks the separate question of whether the age-concordant pairs "
        "are more similar than the other matched pairs.", "",
        "| measure | n finite | matched mean | size-matched null | p_emp | concordant | "
        "discordant | diff 95% CI | MWU p |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    def f(v, spec="{:.4f}"):
        return spec.format(v) if isinstance(v, (int, float)) and np.isfinite(v) else "n/a"

    for name in MEASURES:
        s = m[name]
        ci = ("n/a" if s["diff_ci_low"] is None
              else f"{s['diff_ci_low']:.4f} to {s['diff_ci_high']:.4f}")
        lines.append(
            f"| `{name}` | {s['n_finite']} | {f(s['mean_matched'])} | {f(s['null_mean'])} | "
            f"{f(s['p_emp'], '{:.4g}')} | {f(s['mean_concordant'])} | "
            f"{f(s['mean_discordant'])} | {ci} | {f(s['mwu_p_concordant_vs_not'], '{:.4g}')} |")

    # Which measures had no input at all, and where.  Without this an "n/a" row reads as a
    # negative result when it is actually an absent upstream analysis.
    cov = payload.get("coverage", {})
    missing: dict[str, list[str]] = {}
    for pair, sides in cov.items():
        for side, avail in sides.items():
            for meas, present in avail.items():
                if not present:
                    missing.setdefault(meas, []).append(f"{side} ({pair})")
    lines += ["", "## Data coverage", ""]
    if missing:
        lines.append("A measure is NaN wherever its upstream analysis was never produced "
                     "for that region. **These are absent inputs, not null results.**", )
        lines += ["", "| measure | regions with no input |", "|---|---|"]
        for meas in MEASURES:
            if meas in missing:
                lines.append(f"| `{meas}` | {', '.join(sorted(set(missing[meas])))} |")
        lines += ["", "`celltype_composition` is only run for the IsoGraph artifact tree, "
                      "and `switch_consequence` skips regions with no phenotype-significant "
                      "switch genes; both are upstream properties, not failures of this "
                      "analysis. A pair contributes to a measure only when **both** sides "
                      "have it."]
    else:
        lines.append("All three measures had inputs present on both sides of every pair.")

    sig = [n for n in MEASURES
           if m[n]["n_finite"] > 0
           if (m[n]["p_emp"] is not None and m[n]["p_emp"] < 0.05
               and m[n]["mean_matched"] is not None and m[n]["null_mean"] is not None
               and m[n]["mean_matched"] > m[n]["null_mean"])]
    untested = [n for n in MEASURES if m[n]["n_finite"] == 0]
    lines += ["", "## Verdict", ""]
    if sig:
        lines.append(
            "Matched pairs exceed size-matched random pairs on "
            + ", ".join(f"`{s}`" for s in sig)
            + ". The manuscript may state that the cohorts recover related biological "
              "programs at a higher level of organisation than gene identity, citing these "
              "measures specifically — not as a general claim of replication.")
    else:
        lines.append(
            "No measure with data exceeds its size-matched null. Per the pre-registered "
            "decision rule, the language stays at \"matched modules with concordant age "
            "effects\"; \"replicated programs\" must not be used.")
    if untested:
        lines += ["", "**Untested, not negative:** "
                  + ", ".join(f"`{u}`" for u in untested)
                  + " had no computable pair (see *Data coverage*). The decision rule was "
                    "applied only to the measures that had data; these remain open."]
    (outdir / f"FUNCTIONAL_PRESERVATION__{stem}.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--method", default="isograph", choices=list(METHOD_DIRS))
    p.add_argument("--model", default="linear", choices=["linear", "spline"],
                   help="which module_trust replication table to read")
    p.add_argument("--n-perm", type=int, default=N_PERMUTATIONS)
    p.add_argument("--n-boot", type=int, default=N_BOOT)
    p.add_argument("--seed", type=int, default=SEED)
    args = p.parse_args()
    run(args.method, args.model, args.n_perm, args.seed, args.n_boot)


if __name__ == "__main__":
    main()
