"""Module-level colocalization convergence — do GWAS coloc hits concentrate in modules?

The gene-level coloc (`coloc_prep` -> `coloc_clpp.R`) answers "does this switch gene's
sQTL share a causal variant with the GWAS?". It never asks the module-level question,
and `clpp_results.tsv` carries no module column at all. `scz_age_projection.layer_convergence`
does ask it, but only for SCZ. This CLI generalises that layer to every trait, so the
four aging traits (AD, PD, LBD, ALS) — which are the paper's actual subject — get the
same module-level treatment SCZ already has.

WHY NOT A MODULE-EIGENSWITCH QTL. The obvious version of "module-level genetics" is to
map a QTL for the module summary and colocalize that. It is not worth doing here: a
module eigenswitch is a genome-wide trait with no cis window, and at BrainSEQ n=232 the
genome-wide-significant floor is r^2 ~ 0.135 — a single variant would have to explain
>13% of the variance of an average over hundreds of genes. If it did hit, it would be
one member gene's cis-sQTL leaking into the summary, i.e. a noisier restatement of the
gene-level result. This aggregation asks the same biological question using evidence
that already exists, and costs no new data.

MODULE ASSIGNMENT IS RE-DERIVED, NOT READ FROM THE COLOC OUTPUT. For a bundle gene
source (`aging`), `coloc_prep._resolve_switch_genes` assigns `module_id` via
`drop_duplicates("gene")` — an ARBITRARY contributing analysis's label for a gene that
recurs in three or more. Using that column would be the same class of error as joining
a stale enrichment table (see partition_provenance.py). Instead, genes are re-mapped
onto each source's own `modules.parquet`, exactly as `layer_convergence` does, and the
`source` column records which partition a row belongs to.

GENES VS LOCI. Several colocalizing genes can sit under ONE GWAS peak — gene-dense
regions (17q21/MAPT, APOE) are the classic way a "module converges on disease risk"
claim turns out to be one locus. Both counts are reported, and the leave-one-locus-out
sensitivity drops each module's largest contributing locus and refits. A module whose
signal does not survive that is a single-locus result wearing a module label.
(Note: `scz_age_projection` names its column `n_coloc_loci` but computes
`gene.nunique()` — it counts GENES. This CLI reports the two separately.)

Enrichment is hypergeometric per module against the trait's own tested pool (the genes
that actually received a CLPP test in that trait and are in the source's partition),
BH-corrected across modules within (trait, source). A size-matched permutation gives a
global concentration p, because per-module counts are small and the question "is the
overall distribution more concentrated than chance" is better powered than any single
module.

Writes under 05_genetic_anchoring/_m/module_coloc_convergence/:
  convergence.parquet      one row per (trait, source, module): tested/coloc counts,
                           distinct loci, hypergeometric p + FDR, MAGMA trait P,
                           go_invisible count, leave-one-locus-out worst-case p.
  global.parquet           one row per (trait, source): pool sizes, concentration
                           statistic, permutation p, hypergeometric enrichment of
                           coloc genes in MAGMA-anchored modules.
  MODULE_COLOC_CONVERGENCE.md  the writeup.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.paths import ensure_dir, region_store, stage_out

_COLOC_ROOT = stage_out("anchoring", "coloc")
_MAGMA = stage_out("anchoring.gwas", "magma_results_combined.parquet")
_OUT = stage_out("anchoring", "module_coloc_convergence")
_SEED = 13
_N_PERM = 20000
_MIN_POOL_GENES = 3   # a module with fewer genes in the tested pool cannot be informative

# trait -> coloc output tag. The aging traits share the `aging` bundle gene source;
# SCZ is carried in both its disease-cohort and aging-bundle forms.
TRAIT_TAGS = {
    "AD": "aging__ad",
    "PD": "aging__pd",
    "LBD": "aging__lbd",
    "ALS": "aging__als",
    "SCZ": "aging__scz",
}

# source -> (cohort, region, MAGMA module-name prefix). Mirrors
# scz_age_projection._AGING_SOURCES so the aging traits are reported on exactly the
# partitions SCZ is already reported on.
SOURCES = {
    "gtex_caudate_bg": ("gtex", "caudate_basal_ganglia", "gtex__caudate_basal_ganglia__"),
    "brainseq_caudate": ("brainseq", "caudate", "brainseq__caudate__"),
}


def _bare(s: pd.Series) -> pd.Series:
    return s.astype(str).str.split(".", n=1).str[0]


def load_clpp(tag: str) -> pd.DataFrame:
    """Gene-level coloc results for one trait: gene, LOCUS_ID, colocalized."""
    f = _COLOC_ROOT / tag / "coloc" / "clpp_results.tsv"
    if not f.exists():
        return pd.DataFrame()
    d = pd.read_csv(f, sep="\t")
    d["gene"] = _bare(d["gene"])
    d["colocalized"] = d["colocalized"].astype(str).str.upper().eq("TRUE")
    return d


def load_source_modules(source: str) -> pd.DataFrame:
    """gene -> module_id for one source's own partition (never the coloc table's label)."""
    cohort, region, _ = SOURCES[source]
    m = pd.read_parquet(region_store(cohort, region, "isograph_vae", "modules.parquet"))
    return pd.DataFrame({"gene": _bare(m["gene_id"]), "module_id": m["module_id"].astype(str)})


def magma_p(source: str, trait: str) -> dict[str, float]:
    """module_id -> MAGMA gene-set P for this trait, for the anchored-module contrast."""
    if not _MAGMA.exists():
        return {}
    m = pd.read_parquet(_MAGMA)
    m = m[(m["trait"] == trait) & (m["backend"] == "isograph_vae")]
    prefix = SOURCES[source][2]
    hit = m[m["FULL_NAME"].astype(str).str.startswith(prefix)]
    return {str(n)[len(prefix):]: float(p) for n, p in zip(hit["FULL_NAME"], hit["P"])}


def _hyper(k: int, n: int, K: int, N: int) -> float:
    """P(X >= k) for k coloc genes among n module genes, K coloc in a pool of N."""
    if n == 0 or K == 0 or N == 0:
        return np.nan
    return float(stats.hypergeom.sf(k - 1, N, K, n))


def convergence_for(trait: str, source: str, n_perm: int, rng: np.random.Generator,
                    min_pool_genes: int = _MIN_POOL_GENES) -> tuple[pd.DataFrame, dict]:
    clpp = load_clpp(TRAIT_TAGS[trait])
    if clpp.empty:
        return pd.DataFrame(), {}
    mods = load_source_modules(source)

    # Pool = genes that actually received a CLPP test AND live in this partition. Genes
    # tested but absent from the partition cannot contribute to any module and must not
    # inflate the denominator.
    tested = clpp.merge(mods, on="gene", how="inner")
    if tested.empty:
        return pd.DataFrame(), {}
    pool_genes = tested.drop_duplicates("gene")[["gene", "module_id"]]
    hit_rows = tested[tested["colocalized"]]
    hit_genes = set(hit_rows["gene"])
    N, K = len(pool_genes), len(hit_genes)

    go_inv = _go_invisible_map(TRAIT_TAGS[trait])
    pmap = magma_p(source, trait)
    # locus per hit gene, for the genes-vs-loci distinction and the LOO sensitivity
    gene_locus = (hit_rows.drop_duplicates("gene").set_index("gene")["LOCUS_ID"]
                  if "LOCUS_ID" in hit_rows.columns else pd.Series(dtype=object))

    rows = []
    for mid, grp in pool_genes.groupby("module_id"):
        genes = set(grp["gene"])
        hits = genes & hit_genes
        n = len(genes)
        loci = {gene_locus.get(g) for g in hits} - {None}
        # Leave-one-locus-out: drop this module's largest contributing locus (and the
        # matching genes from the pool's hit count) and refit. Reports the WORST case.
        loo_p = np.nan
        if loci:
            worst = 0.0
            for drop in loci:
                keep = {g for g in hits if gene_locus.get(g) != drop}
                dropped_all = {g for g in hit_genes if gene_locus.get(g) == drop}
                worst = max(worst, _hyper(len(keep), n, K - len(dropped_all), N))
            loo_p = float(worst)
        rows.append({
            "trait": trait, "source": source, "module_id": mid,
            "n_module_genes_in_pool": n,
            "n_coloc_genes": len(hits),
            "n_coloc_loci": len(loci),
            "n_go_invisible": int(sum(go_inv.get(g, False) for g in hits)),
            "hyperg_p": _hyper(len(hits), n, K, N),
            "loo_worst_p": loo_p,
            "magma_p": pmap.get(mid, np.nan),
            "coloc_genes": ", ".join(sorted(hits)),
        })
    conv = pd.DataFrame(rows)
    conv = conv[conv["n_coloc_genes"] > 0].copy()
    if not conv.empty:
        # A module contributing 1-2 genes to the tested pool cannot be informative: if its
        # single gene colocalizes the hypergeometric is "significant" (p = K/N) purely
        # because the module is tiny in the pool, and it consumes an FDR slot that then
        # penalises the real tests. LBD/brainseq_caudate M012 is exactly this — one gene
        # in the pool (SNCA), one hit, p=0.033. Such rows are kept for transparency,
        # flagged, and excluded from the multiple-testing correction.
        conv["testable"] = conv["n_module_genes_in_pool"] >= min_pool_genes
        conv["hyperg_fdr"] = np.nan
        t = conv["testable"].to_numpy()
        if t.any():
            conv.loc[t, "hyperg_fdr"] = _bh(conv.loc[t, "hyperg_p"].to_numpy())
        conv = conv.sort_values(["testable", "n_coloc_genes", "hyperg_p"],
                                ascending=[False, False, True])

    glob = _global_stats(pool_genes, hit_genes, pmap, trait, source, n_perm, rng)
    return conv, glob


def _global_stats(pool: pd.DataFrame, hit_genes: set, pmap: dict, trait: str,
                  source: str, n_perm: int, rng: np.random.Generator) -> dict:
    """Concentration permutation + anchored-module hypergeometric.

    Concentration statistic = sum of squared per-module hit counts, which is large when
    hits pile into few modules and small when they spread. The null redraws the same
    number of hit genes uniformly from the pool, so module SIZES are held fixed by
    construction — a big module attracting hits by size alone is not evidence.
    """
    labels = pool["module_id"].to_numpy()
    idx_hits = np.flatnonzero(pool["gene"].isin(hit_genes).to_numpy())
    N, K = len(pool), len(idx_hits)
    if K == 0:
        return {"trait": trait, "source": source, "n_pool": N, "n_coloc_genes": 0}
    codes, _ = pd.factorize(labels)
    n_mod = codes.max() + 1
    obs = float((np.bincount(codes[idx_hits], minlength=n_mod) ** 2).sum())
    null = np.empty(n_perm)
    for i in range(n_perm):
        pick = rng.choice(N, size=K, replace=False)
        null[i] = (np.bincount(codes[pick], minlength=n_mod) ** 2).sum()
    conc_p = float((np.sum(null >= obs) + 1) / (n_perm + 1))

    # Are coloc genes concentrated in MAGMA-anchored modules (trait P < 0.05)?
    #
    # THE DENOMINATOR IS THE WHOLE RESULT. Both denominators are emitted because they
    # disagree and the disagreement is the finding. `scz_age_projection` uses ALL module
    # genes, against which SCZ coloc genes look strongly enriched in anchored modules
    # (15/31 = 48% vs a 25% background, P=0.004). But a gene can only colocalize if it
    # was TESTED, i.e. if it sat under a GWAS peak with a QTL credible set — and anchored
    # modules are MAGMA-enriched for that same GWAS, so their genes enter the tested pool
    # preferentially. Conditioning on the tested pool, the anchored background is already
    # 45% and 48% is null (P=0.40). MAGMA anchoring and coloc testing select on the same
    # signal; the published enrichment is largely that shared ascertainment.
    anchored = {m for m, p in pmap.items() if np.isfinite(p) and p < 0.05}
    anch_mask = pool["module_id"].isin(anchored).to_numpy()
    Kk = int(anch_mask.sum())
    x = int(anch_mask[idx_hits].sum())
    all_mods = load_source_modules(source)
    N_all = int(all_mods["gene"].nunique())
    K_all = int(all_mods.loc[all_mods["module_id"].isin(anchored), "gene"].nunique())
    return {
        "trait": trait, "source": source, "n_pool": N, "n_coloc_genes": K,
        "n_modules_with_coloc": int(len(set(labels[idx_hits]))),
        "concentration_obs": obs, "concentration_null_mean": float(null.mean()),
        "concentration_p": conc_p,
        "n_anchored_modules": len(anchored),
        "frac_coloc_anchored": x / K if K else np.nan,
        # correct denominator: only genes that could have been a hit
        "frac_pool_anchored": Kk / N if N else np.nan,
        "anchored_hyperg_p": _hyper(x, K, Kk, N),
        # published denominator: every module gene, including untested ones
        "frac_allgenes_anchored": K_all / N_all if N_all else np.nan,
        "anchored_hyperg_p_allgenes_denom": _hyper(x, K, K_all, N_all),
    }


def _go_invisible_map(tag: str) -> dict[str, bool]:
    f = _COLOC_ROOT / tag / "candidate_genes.tsv"
    if not f.exists():
        return {}
    d = pd.read_csv(f, sep="\t")
    return dict(zip(_bare(d["gene"]), d["go_invisible"].astype(str).str.lower() == "true"))


def _bh(p: np.ndarray) -> np.ndarray:
    ok = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    if not ok.any():
        return out
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order] * len(q) / (np.arange(len(q)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    res = np.empty(len(q))
    res[order] = ranked
    out[ok] = res
    return out


def run(traits: list[str], n_perm: int = _N_PERM,
        min_pool_genes: int = _MIN_POOL_GENES) -> pd.DataFrame:
    rng = np.random.default_rng(_SEED)
    convs, globs = [], []
    for trait in traits:
        for source in SOURCES:
            c, g = convergence_for(trait, source, n_perm, rng, min_pool_genes)
            if not c.empty:
                convs.append(c)
            if g:
                globs.append(g)
            n = 0 if c.empty else len(c)
            print(f"[{trait}/{source}] modules with coloc: {n}", flush=True)
    out = ensure_dir(_OUT)
    conv = pd.concat(convs, ignore_index=True) if convs else pd.DataFrame()
    glob = pd.DataFrame(globs)
    conv.to_parquet(out / "convergence.parquet", index=False, compression="zstd")
    glob.to_parquet(out / "global.parquet", index=False, compression="zstd")
    _write_report(out, conv, glob)
    return conv


def _md(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    rows = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    rows += ["| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |"
             for r in df.itertuples(index=False)]
    return "\n".join(rows)


def _write_report(out: Path, conv: pd.DataFrame, glob: pd.DataFrame) -> None:
    L = [
        "# Module-level colocalization convergence",
        "",
        "Do GWAS colocalizing switch genes concentrate in particular IsoGraph co-switch "
        "modules, beyond what module sizes alone predict? Generalises the SCZ-only layer "
        "in `scz_age_projection.py` to every trait, so the four aging traits are reported "
        "on the same partitions and with the same statistics.",
        "",
        "Reproduce: `python -m isograph_benchmark.real_data.module_coloc_convergence`",
        "",
        "## Reading",
        "",
        "- **Per-module counts are small.** With single-digit colocalizing genes per trait, "
        "the per-module hypergeometric is underpowered; the global concentration "
        "permutation is the better-powered test and should lead.",
        "- **`n_coloc_genes` vs `n_coloc_loci`.** Several genes can sit under one GWAS "
        "peak. A module whose genes collapse to one locus is a single-locus result, not a "
        "converging program. `loo_worst_p` drops each contributing locus in turn and "
        "reports the worst case — read it before believing any module row.",
        "- Module assignment is re-derived from each source's own `modules.parquet`; the "
        "`module_id` column in the coloc outputs is an arbitrary contributing analysis's "
        "label for bundle gene sources and is deliberately not used.",
        "- The permutation holds module sizes fixed by construction, so a large module "
        "collecting hits in proportion to its size is not evidence.",
        "",
    ]
    if not glob.empty:
        show = glob.copy()
        for c in ("concentration_obs", "concentration_null_mean", "frac_pool_anchored",
                  "frac_coloc_anchored"):
            if c in show:
                show[c] = show[c].round(3)
        for c in ("concentration_p", "anchored_hyperg_p"):
            if c in show:
                show[c] = show[c].apply(lambda p: f"{p:.3g}" if pd.notna(p) else "NA")
        L += ["## Global concentration per (trait, source)", "", _md(show), ""]
    if not conv.empty:
        top = conv[conv["n_coloc_genes"] >= 1].copy()
        for c in ("hyperg_p", "hyperg_fdr", "loo_worst_p", "magma_p"):
            top[c] = top[c].apply(lambda p: f"{p:.3g}" if pd.notna(p) else "NA")
        cols = ["trait", "source", "module_id", "testable", "n_module_genes_in_pool",
                "n_coloc_genes", "n_coloc_loci", "n_go_invisible", "hyperg_p",
                "hyperg_fdr", "loo_worst_p", "magma_p"]
        L += ["## Modules carrying colocalizing switch genes", "", _md(top[cols]), ""]
    (out / "MODULE_COLOC_CONVERGENCE.md").write_text("\n".join(L) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="Module-level coloc convergence per trait.")
    p.add_argument("--traits", nargs="+", default=list(TRAIT_TAGS),
                   choices=list(TRAIT_TAGS))
    p.add_argument("--n-perm", type=int, default=_N_PERM)
    p.add_argument("--min-pool-genes", type=int, default=_MIN_POOL_GENES,
                   help="modules with fewer genes in the tested pool are reported but "
                        "flagged untestable and excluded from the FDR")
    args = p.parse_args()
    run(args.traits, args.n_perm, args.min_pool_genes)


if __name__ == "__main__":
    main()
