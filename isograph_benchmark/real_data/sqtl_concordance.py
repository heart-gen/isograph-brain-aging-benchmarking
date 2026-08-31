"""Direction concordance between GTEx sQTLs and IsoGraph isoform switches.

The companion `qtl_anchoring.py` asks *whether* co-switch module genes carry
splicing QTLs (a matched enrichment). This module asks the stronger, orthogonal
question: when a gene has an sQTL, does the genetic control of splicing point in
the *same direction* as IsoGraph's inferred isoform switch? A directional match
is much harder to produce by a size/detectability confound than mere enrichment,
so it hardens the genetic-anchoring headline.

The sign problem, and why the test is relative. A LeafCutter sQTL slope is signed
relative to the variant's (arbitrary) alt allele; an IsoGraph transcript switch is
signed relative to the module's (arbitrary) switch-axis orientation. No *absolute*
per-transcript direction comparison is therefore meaningful. What is meaningful is
the *within-gene relative* structure: transcripts that IsoGraph places on opposite
sides of the switch should, if the switch is genetically driven, have their introns
pushed in opposite directions by the same sQTL variant. We measure that as a
per-gene rank correlation between two within-gene direction vectors, and because a
joint global flip of either vector is meaningless we aggregate the *magnitude*
|rho| and calibrate it against a within-gene transcript-label permutation null.

Construction, per gene with an sQTL and >= MIN_TX shared transcripts:
  * IsoGraph direction d_t: signed usage change of transcript t along the module
    switch axis, from the per-module high_vs_low table (mean_high - mean_low).
  * Genetic direction g_t: for the gene's lead sQTL variant (smallest nominal p),
    the net LeafCutter slope over the introns that transcript t splices out. An
    intron maps to a transcript when the transcript contains it as an intron
    (exon_end, next_exon_start) -- verified to join ~80% of GTEx brain sQTL introns
    to the GENCODE v47 cache. slope > 0 = alt allele raises that intron's excision
    ratio, favouring transcripts that use the intron, so g_t = sign(sum of slopes).
  * Concordance rho_g = Spearman(d_t, g_t) across the gene's shared transcripts.

Aggregation per module set (same sets as qtl_anchoring: all / pheno_sig /
go_invisible / go_visible): mean |rho_g|, and a permutation p from shuffling the
transcript->g_t assignment within each gene (seeded, deterministic). On-thesis
result: phenotype-associated and GO-invisible module sets show mean |rho| above
the null, i.e. sQTL and switch directions are coupled, not merely co-located.

Reads only saved artifacts + the copied GTEx xQTL catalog. Genes matched on
unversioned Ensembl id. Writes under <artifact-parent>/_m/:
  sqtl_concordance.parquet  -- one row per module_set: n_genes, mean_abs_rho,
      mean_signed_rho, frac genes with |rho|>=0.5, permutation p, null mean.
  sqtl_concordance_pergene.parquet -- one row per (module_set, gene): rho, n_tx.
  SQTL_CONCORDANCE.md  -- Manubot summary.
  sqtl_concordance.json -- run parameters.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.interpret_modules import DEFAULT_GTF_CACHE
from isograph_benchmark.real_data.qtl_anchoring import (
    _GTEX_TISSUE,
    _bare,
    build_gene_sets,
    resolve_tissue,
)
from isograph_benchmark.real_data.sweep_leiden import _artifact_dir

DEFAULT_XQTL_DIR = rel("inputs", "raw", "gtex_v11", "xqtl")
_QVAL = 0.05
_MIN_TX = 3          # genes need >= this many shared transcripts for a stable rho
_N_PERM = 2000       # within-gene label permutations for the null
_SEED = 13


def transcript_introns(gtf_cache: Path) -> pd.DataFrame:
    """One row per (transcript, intron): chrom, istart=exon_end, iend=next_exon_start.

    Intron coordinates follow the LeafCutter phenotype convention (donor exon end,
    acceptor exon start) so they join the GTEx sQTL phenotype_id directly.
    """
    g = pd.read_parquet(gtf_cache, columns=["transcript_id", "gene_id", "chrom", "feature", "start", "end"])
    ex = g[g["feature"] == "exon"].sort_values(["transcript_id", "start"])
    ex["next_start"] = ex.groupby("transcript_id")["start"].shift(-1)
    ex["next_tx"] = ex.groupby("transcript_id")["transcript_id"].shift(-1)
    intr = ex[ex["next_tx"].notna()].copy()
    intr["istart"] = intr["end"].astype("int64")
    intr["iend"] = intr["next_start"].astype("int64")
    intr["gene"] = _bare(intr["gene_id"])
    return intr[["transcript_id", "gene", "chrom", "istart", "iend"]]


def load_lead_sqtl(xqtl_dir: Path, tissue: str) -> pd.DataFrame:
    """Lead sQTL per intron: smallest nominal-p variant, with its signed slope.

    Returns one row per (gene, intron): chrom, istart, iend, variant_id, slope.
    """
    path = xqtl_dir / f"{tissue}.v11.sQTLs.signif_pairs.parquet"
    d = pd.read_parquet(path, columns=["phenotype_id", "variant_id", "slope", "pval_nominal"])
    parts = d["phenotype_id"].str.split(":", expand=True)
    d["chrom"] = parts[0]
    d["istart"] = parts[1].astype("int64")
    d["iend"] = parts[2].astype("int64")
    d["gene"] = _bare(parts[4])
    # lead variant per intron
    idx = d.groupby("phenotype_id")["pval_nominal"].idxmin()
    lead = d.loc[idx, ["gene", "chrom", "istart", "iend", "variant_id", "slope"]].reset_index(drop=True)
    return lead


def genetic_direction(gene: str, lead: pd.DataFrame, introns: pd.DataFrame) -> pd.Series:
    """Net sQTL slope per transcript for a gene's single lead variant.

    Picks the variant that leads the most introns for the gene (ties -> the one
    with the largest summed |slope|), then sums that variant's intron slopes onto
    each transcript that splices those introns. Returns Series indexed by
    transcript_id (bare-gene transcripts), empty if no mappable intron.
    """
    gl = lead[lead["gene"] == gene]
    if gl.empty:
        return pd.Series(dtype=float)
    # choose one variant so all transcript directions share a common allele reference
    order = (gl.assign(abss=gl["slope"].abs())
               .groupby("variant_id")
               .agg(n=("istart", "size"), s=("abss", "sum"))
               .sort_values(["n", "s"], ascending=False))
    variant = order.index[0]
    gv = gl[gl["variant_id"] == variant]
    gi = introns[introns["gene"] == gene]
    joined = gv.merge(gi, on=["gene", "chrom", "istart", "iend"], how="inner")
    if joined.empty:
        return pd.Series(dtype=float)
    return joined.groupby("transcript_id")["slope"].sum()


def _isograph_direction(mod_interpret_dir: Path, module_id: str) -> pd.Series:
    """Signed switch-axis usage change per transcript for a module (mean_high-mean_low)."""
    path = mod_interpret_dir / module_id / "high_vs_low_table.parquet"
    if not path.exists():
        return pd.Series(dtype=float)
    hl = pd.read_parquet(path, columns=["feature_id", "delta", "tstat"])
    hl = hl[hl["feature_id"].astype(str).str.startswith("ENST") & hl["delta"].notna()]
    return hl.set_index("feature_id")["delta"]


def _gene_of(transcript_ids: pd.Index, introns: pd.DataFrame) -> dict[str, str]:
    m = introns.drop_duplicates("transcript_id").set_index("transcript_id")["gene"]
    return m.to_dict()


def per_gene_concordance(iso_dir: Path, gene_sets: dict[str, set[str]],
                         lead: pd.DataFrame, introns: pd.DataFrame,
                         rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute rho_g for every gene reachable in any module set, then tag by set."""
    modules = pd.read_parquet(iso_dir / "modules.parquet")
    modules["gene"] = _bare(modules["gene_id"])
    modules["module_id"] = modules["module_id"].astype(str)
    interpret_dir = iso_dir / "module_interpret"
    tx_gene = _gene_of(introns["transcript_id"], introns)

    # a gene's IsoGraph direction comes from whichever interpreted module it sits in
    gene_rho: dict[str, tuple[float, int]] = {}
    for module_id, grp in modules.groupby("module_id"):
        d_all = _isograph_direction(interpret_dir, module_id)
        if d_all.empty:
            continue
        d_gene = pd.Series(d_all.index.map(lambda t: tx_gene.get(t, None)), index=d_all.index)
        for gene in set(grp["gene"]):
            if gene in gene_rho:
                continue
            d_t = d_all[d_gene == gene]
            g_t = genetic_direction(gene, lead, introns)
            shared = d_t.index.intersection(g_t.index)
            if len(shared) < _MIN_TX:
                continue
            dv = d_t.loc[shared].to_numpy(float)
            gv = g_t.loc[shared].to_numpy(float)
            if np.ptp(dv) == 0 or np.ptp(gv) == 0:
                continue
            rho = spearmanr(dv, gv).statistic
            if np.isnan(rho):
                continue
            gene_rho[gene] = (float(rho), int(len(shared)))

    rho_df = pd.DataFrame(
        [{"gene": g, "rho": r, "n_tx": n} for g, (r, n) in gene_rho.items()]
    )

    rows, pergene_rows = [], []
    for set_name, genes in gene_sets.items():
        sub = rho_df[rho_df["gene"].isin(genes)] if not rho_df.empty else rho_df
        n = len(sub)
        if n == 0:
            rows.append({"module_set": set_name, "n_genes": 0, "mean_abs_rho": np.nan,
                         "mean_signed_rho": np.nan, "frac_strong": np.nan,
                         "perm_mean_abs_rho": np.nan, "pvalue": np.nan})
            continue
        obs = float(sub["rho"].abs().mean())
        perm_means = _perm_null(sub, rng)
        pval = float((np.sum(np.asarray(perm_means) >= obs) + 1) / (len(perm_means) + 1))
        rows.append({"module_set": set_name, "n_genes": n,
                     "mean_abs_rho": round(obs, 4),
                     "mean_signed_rho": round(float(sub["rho"].mean()), 4),
                     "frac_strong": round(float((sub["rho"].abs() >= 0.5).mean()), 4),
                     "perm_mean_abs_rho": round(float(np.mean(perm_means)), 4),
                     "pvalue": pval})
        for _, r in sub.iterrows():
            pergene_rows.append({"module_set": set_name, "gene": r["gene"],
                                 "rho": round(float(r["rho"]), 4), "n_tx": int(r["n_tx"])})
    return pd.DataFrame(rows), pd.DataFrame(pergene_rows)


def _null_abs_rho(n: int, rng: np.random.Generator) -> np.ndarray:
    """|Spearman rho| for _N_PERM random rankings of length n, vectorised.

    Under the null the two within-gene direction vectors are unrelated, so rho is
    the Spearman of arange(n) against a random permutation. For distinct integer
    ranks that is the closed form rho = 1 - 6*sum(d^2)/(n*(n^2-1)), which avoids the
    per-call cost of scipy.spearmanr entirely.
    """
    if n < 2:
        return np.zeros(_N_PERM)
    perms = np.argsort(rng.random((_N_PERM, n)), axis=1)
    a = np.arange(n)
    d2 = np.sum((a[None, :] - perms) ** 2, axis=1)
    return np.abs(1.0 - 6.0 * d2 / (n * (n * n - 1)))


def _perm_null(sub: pd.DataFrame, rng: np.random.Generator) -> list[float]:
    """Null distribution of mean |rho| for a gene set.

    Each gene contributes an independent null |rho| draw at its own n_tx (so the
    small-n_tx genes that mechanically inflate |rho| are matched in the null), then
    the per-permutation mean is taken over the set's genes.
    """
    n_tx = sub["n_tx"].to_numpy(int)
    mat = np.empty((len(n_tx), _N_PERM))
    for j, n in enumerate(n_tx):
        mat[j] = _null_abs_rho(int(n), rng)
    return mat.mean(axis=0).tolist()


def run(analysis: str, region: str | None, variant: str, fdr: float, xqtl_dir: Path,
        tissue: str | None, gtf_cache: Path) -> pd.DataFrame:
    iso_dir = _artifact_dir(analysis, region, variant)  # the isograph_vae modules dir
    tissue = tissue or resolve_tissue(analysis, region)
    enrich_path = iso_dir.parent / "module_enrichment" / "isograph_modules.parquet"
    gene_sets = build_gene_sets(iso_dir, enrich_path, fdr)

    introns = transcript_introns(gtf_cache)
    lead = load_lead_sqtl(xqtl_dir, tissue)
    rng = np.random.default_rng(_SEED)
    summary, pergene = per_gene_concordance(iso_dir, gene_sets, lead, introns, rng)
    summary.insert(0, "tissue", tissue)
    summary.insert(0, "region", region or "")
    summary.insert(0, "analysis", analysis)

    out_dir = ensure_dir(iso_dir.parent)
    summary.to_parquet(out_dir / "sqtl_concordance.parquet", index=False, compression="zstd")
    if not pergene.empty:
        pergene.to_parquet(out_dir / "sqtl_concordance_pergene.parquet", index=False, compression="zstd")
    (out_dir / "sqtl_concordance.json").write_text(json.dumps(
        {"analysis": analysis, "region": region, "tissue": tissue, "variant": variant,
         "fdr": fdr, "min_tx": _MIN_TX, "n_perm": _N_PERM, "seed": _SEED,
         "n_lead_sqtl_introns": int(len(lead)),
         "module_set_sizes": {k: len(v) for k, v in gene_sets.items()}}, indent=2))
    _write_report(out_dir, analysis, region, tissue, summary)
    return summary


def _write_report(out_dir: Path, analysis: str, region: str | None, tissue: str,
                  summary: pd.DataFrame) -> None:
    show = summary[["module_set", "n_genes", "mean_abs_rho", "perm_mean_abs_rho",
                    "frac_strong", "mean_signed_rho", "pvalue"]].copy()
    show["pvalue"] = show["pvalue"].apply(lambda p: f"{p:.3g}" if pd.notna(p) else "NA")
    cols = list(show.columns)
    tbl = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    tbl += ["| " + " | ".join(str(v) for v in r) + " |" for r in show.itertuples(index=False)]
    lines = [
        f"# sQTL direction concordance — IsoGraph switches vs GTEx {tissue} sQTLs "
        f"({analysis}" + (f"/{region}" if region else "") + ")",
        "",
        "Per gene with a lead sQTL and >= {} shared transcripts, the within-gene "
        "Spearman correlation (rho) between IsoGraph's switch-axis usage change "
        "(mean_high - mean_low) and the sQTL lead variant's net intron slope per "
        "transcript. Aggregated as mean |rho| per module set; `pvalue` is a within-"
        "gene rank-permutation null (seed {}, {} draws).".format(_MIN_TX, _SEED, _N_PERM),
        "",
        "Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance "
        f"--analysis {analysis}" + (f" --region {region}" if region else "") + "`",
        "",
        "## Concordance by module set",
        "",
        "\n".join(tbl),
        "",
        "## Reading",
        "",
        "- mean |rho| above the permutation null would indicate the sQTL and the "
        "isoform switch move the *same* transcripts in a coupled direction. Note this "
        "within-gene test is underpowered: a single lead sQTL variant often tags "
        "introns with near-constant per-transcript direction, so per-cohort |rho| "
        "hugs the null; the pooled meta and the diagnosis live in "
        "`05_genetic_anchoring/_m/sqtl_concordance_meta/`.",
        "- The statistic is flip-invariant by construction: the sQTL allele reference "
        "and the module switch-axis orientation are both arbitrary, so only the "
        "*relative* within-gene direction structure is testable. The disease-anchored "
        "directional test is colocalization (`sqtl_coloc`), not this one.",
        "- Scope: concordance is evaluated on genes whose introns map to the GENCODE "
        "cache (~80% of GTEx brain sQTL introns) and that carry >= {} interpretable "
        "switch transcripts.".format(_MIN_TX),
    ]
    (out_dir / "SQTL_CONCORDANCE.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="sQTL/isoform-switch direction concordance.")
    p.add_argument("--analysis", required=True, help="brainseq-sczd, brainseq-aging, gtex-aging")
    p.add_argument("--region", default=None)
    p.add_argument("--variant", default="standard")
    p.add_argument("--fdr", type=float, default=_QVAL)
    p.add_argument("--xqtl-dir", default=str(DEFAULT_XQTL_DIR))
    p.add_argument("--tissue", default=None, help="override GTEx tissue prefix")
    p.add_argument("--gtf-cache", default=str(DEFAULT_GTF_CACHE))
    args = p.parse_args()
    summary = run(args.analysis, args.region, args.variant, args.fdr,
                  Path(args.xqtl_dir), args.tissue, Path(args.gtf_cache))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
