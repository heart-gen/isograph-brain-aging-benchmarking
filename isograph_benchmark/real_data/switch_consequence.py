"""Do IsoGraph co-switch modules align with PROTEIN-CONSEQUENTIAL isoform differences?

The genetic-anchoring layer shows the switch modules are real, coordinated, and heritable.
This asks the functional-content question a reviewer will: coordinated changes in *what*?
For every switching transcript pair IsoGraph reports, we score the structural consequence of
the switch (first/last/internal exon change, CDS change, UTR change, biotype switch, coding-
status change) with the SAME pairwise comparator that built the module structural annotations
(`isograph.explain.structure.compare_pair`), and test whether the switch pairs are enriched
for coding-consequential differences relative to a WITHIN-GENE null of random transcript pairs.

Why a within-gene null: the null pairs are drawn from each switch gene's own transcript set, so
transcript count, gene length, and isoform-annotation depth — the things that trivially inflate
any structural difference — cancel exactly. The test isolates whether IsoGraph's switch axis
preferentially selects coding-consequential isoform pairs. Results are stratified by GO-invisible
vs GO-visible modules: if the enrichment holds (or strengthens) in GO-invisible modules, the
protein-remodeling biology is exactly what abundance/pathway enrichment cannot see.

One region+resolution per invocation (mirrors interpret_modules). Deterministic given --seed.

Output (<artifact_dir>/switch_consequence/): pair_consequence.parquet (per observed switch
pair: flags + module tags), consequence_enrichment.parquet (per class x GO-stratum: observed
rate, permutation-null mean, enrichment, empirical p), SWITCH_CONSEQUENCE.md.
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

from isograph.explain.structure import compare_pair, parse_gtf
from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.coloc_prep import load_switch_genes
from isograph_benchmark.real_data.interpret_modules import (
    DEFAULT_GTF_CACHE,
    DEFAULT_GTF_PATH,
)
from isograph_benchmark.real_data.qtl_anchoring import _bare

# Boolean consequence classes (from compare_pair) + composites we derive.
#   coding_consequence = CDS change or coding-status change.
#   nmd_switch         = one switch isoform is a predicted NMD target, the other is not
#                        (a switch INTO/OUT OF nonsense-mediated decay).
_CLASSES = ["first_exon_changed", "last_exon_changed", "internal_exon_difference",
            "cds_changed", "utr_changed", "biotype_switch", "coding_status_change",
            "coding_consequence", "nmd_switch"]


def _strip_version(tx_id: str) -> str:
    return tx_id.split(".", 1)[0]


def _ordered_exons(rec) -> list[tuple[int, int]]:
    """Exons in transcript 5'->3' order (ascending genomic on +, descending on -)."""
    ex = sorted(rec.exons)
    return ex if rec.strand == "+" else list(reversed(ex))


def _spliced_offset(pos: int, ordered: list[tuple[int, int]], strand: str) -> int | None:
    """Distance (nt) from the transcript 5' end to a genomic position along the spliced
    mRNA, or None if the position is not exonic."""
    off = 0
    for (s, e) in ordered:
        if s <= pos <= e:
            return off + ((pos - s) if strand == "+" else (e - pos))
        off += (e - s + 1)
    return None


def _is_nmd(rec) -> bool:
    """Predicted NMD target by the 50-nt rule: a coding transcript whose stop codon lies
    >50 nt upstream (5') of the last exon-exon junction along the spliced mRNA.

    Non-coding transcripts and single-exon transcripts are not NMD targets by this rule.
    """
    if not rec.cds:
        return False
    ordered = _ordered_exons(rec)
    if len(ordered) < 2:
        return False
    cds_sorted = sorted(rec.cds)
    # 3' end of the CDS (translation stop) in transcript orientation.
    stop_g = cds_sorted[-1][1] if rec.strand == "+" else cds_sorted[0][0]
    stop_spliced = _spliced_offset(stop_g, ordered, rec.strand)
    if stop_spliced is None:
        return False
    # last exon-exon junction = spliced coordinate of the final exon's 5' base.
    junction_spliced = sum(e - s + 1 for (s, e) in ordered[:-1])
    return (junction_spliced - stop_spliced) > 50


def _gene_transcripts(tx_db: dict) -> dict[str, list[str]]:
    """bare gene_id -> list of its transcript_ids present in the GTF."""
    by_gene: dict[str, list[str]] = {}
    for tx_id, rec in tx_db.items():
        by_gene.setdefault(_strip_version(rec.gene_id), []).append(tx_id)
    return by_gene


def _pair_flags(tx_db: dict, t1: str, t2: str) -> dict | None:
    """compare_pair -> the boolean classes + a coding_consequence composite, or None if a
    transcript is absent from the GTF."""
    if t1 not in tx_db or t2 not in tx_db:
        return None
    r = compare_pair(tx_db[t1], tx_db[t2])
    out = {c: bool(r[c]) for c in _CLASSES if c in r}
    out["coding_consequence"] = bool(r["cds_changed"]) or bool(r["coding_status_change"])
    # switch INTO/OUT OF NMD: exactly one isoform is a predicted NMD target.
    out["nmd_switch"] = _is_nmd(tx_db[t1]) != _is_nmd(tx_db[t2])
    return out


def _observed_pairs(artifact_dir: Path, switch_genes: pd.DataFrame) -> pd.DataFrame:
    """Observed IsoGraph switch pairs restricted to phenotype-significant switch genes,
    tagged with module + GO-invisible."""
    sp = pd.read_parquet(artifact_dir / "module_interpret" / "structure_switch_pairs.parquet")
    sp["gene"] = _bare(sp["gene_id"])
    # a gene can sit in >1 module; take the GO-invisible tag as True if ANY of its
    # phenotype-significant modules is GO-invisible (the analysis unit here is the gene).
    tag = (switch_genes.groupby("gene")["go_invisible"].max().reset_index())
    sp = sp.merge(tag, on="gene", how="inner")
    return sp[["gene", "transcript_id_1", "transcript_id_2", "go_invisible"]]


def run(region: str, artifact_dir: Path, fdr: float, n_perm: int, seed: int) -> pd.DataFrame:
    iso_dir = artifact_dir
    switch_genes = load_switch_genes(iso_dir, fdr)
    if switch_genes.empty:
        print(f"{region}: no phenotype-significant switch genes; skipping.")
        return pd.DataFrame()
    obs = _observed_pairs(artifact_dir, switch_genes)
    if obs.empty:
        print(f"{region}: no switch pairs for phenotype-significant genes; skipping.")
        return pd.DataFrame()

    tx_db = parse_gtf(DEFAULT_GTF_PATH, cache=DEFAULT_GTF_CACHE)
    gene_tx = _gene_transcripts(tx_db)

    # Observed consequence per switch pair.
    rows = []
    for r in obs.itertuples():
        f = _pair_flags(tx_db, r.transcript_id_1, r.transcript_id_2)
        if f is None:
            continue
        f.update(gene=r.gene, go_invisible=bool(r.go_invisible),
                 transcript_id_1=r.transcript_id_1, transcript_id_2=r.transcript_id_2)
        rows.append(f)
    pair_df = pd.DataFrame(rows)
    if pair_df.empty:
        print(f"{region}: no switch pairs resolvable in the GTF; skipping.")
        return pd.DataFrame()

    # Per-gene pool of ALL annotated transcript-pair consequences (cache once), for the
    # within-gene null. Genes with <2 GTF transcripts contribute no null pairs.
    rng = np.random.default_rng(seed)
    pool: dict[str, list[dict]] = {}
    for gene in pair_df["gene"].unique():
        txs = gene_tx.get(gene, [])
        pairs = [p for p in ((_pair_flags(tx_db, a, b))
                             for a, b in itertools.combinations(txs, 2))
                 if p is not None]
        if pairs:
            pool[gene] = pairs

    # observed pair count per gene (to match the null draw size per gene)
    n_obs_per_gene = pair_df.groupby("gene").size().to_dict()

    enrich = []
    for stratum, sub in _strata(pair_df):
        genes = sub["gene"].unique()
        n_pairs = len(sub)
        if n_pairs == 0:
            continue
        obs_rate = {c: float(sub[c].mean()) for c in _CLASSES}
        # permutation null: per gene, draw n_obs_per_gene[gene] random pairs from its pool
        null_rates = {c: np.empty(n_perm) for c in _CLASSES}
        draw_genes = [g for g in genes if g in pool]
        for i in range(n_perm):
            sampled: list[dict] = []
            for g in draw_genes:
                p = pool[g]
                k = n_obs_per_gene.get(g, 0)
                idx = rng.integers(0, len(p), size=k)
                sampled.extend(p[j] for j in idx)
            if not sampled:
                for c in _CLASSES:
                    null_rates[c][i] = np.nan
                continue
            s = pd.DataFrame(sampled)
            for c in _CLASSES:
                null_rates[c][i] = float(s[c].mean())
        for c in _CLASSES:
            nn = null_rates[c][~np.isnan(null_rates[c])]
            null_mean = float(nn.mean()) if nn.size else np.nan
            # one-sided empirical p: P(null >= observed)
            p_emp = (float((nn >= obs_rate[c]).sum() + 1) / (nn.size + 1)
                     if nn.size else np.nan)
            enrich.append({
                "region": region, "stratum": stratum, "consequence": c,
                "n_obs_pairs": n_pairs, "n_genes": len(genes),
                "obs_rate": obs_rate[c], "null_mean": null_mean,
                "enrichment": (obs_rate[c] / null_mean
                               if null_mean and null_mean > 0 else np.nan),
                "p_emp": p_emp,
            })

    enrich_df = pd.DataFrame(enrich)
    out_dir = ensure_dir(artifact_dir / "switch_consequence")
    pair_df.to_parquet(out_dir / "pair_consequence.parquet", index=False)
    enrich_df.to_parquet(out_dir / "consequence_enrichment.parquet", index=False)
    _write_report(out_dir, region, enrich_df, pair_df)
    sig = enrich_df[(enrich_df["consequence"] == "coding_consequence")]
    msg = "; ".join(f"{r.stratum} enrich {r.enrichment:.2f} p={r.p_emp:.3g}"
                    for r in sig.itertuples())
    print(f"{region}: {len(pair_df)} switch pairs; coding_consequence {msg}")
    return enrich_df


def _strata(pair_df: pd.DataFrame):
    """Yield (label, subframe): all switch pairs, then GO-invisible / GO-visible."""
    yield "all", pair_df
    yield "go_invisible", pair_df[pair_df["go_invisible"]]
    yield "go_visible", pair_df[~pair_df["go_invisible"]]


def _write_report(out_dir: Path, region: str, enrich: pd.DataFrame,
                  pair_df: pd.DataFrame) -> None:
    lines = [
        f"# Switch coding-consequence enrichment — {region}", "",
        "Per consequence class: fraction of IsoGraph switch pairs with the change "
        "(`obs_rate`) vs a within-gene permutation null of random transcript pairs "
        "(`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. "
        "`coding_consequence` = CDS change or coding-status change. Stratified by whether the "
        "switch gene's phenotype-significant module is GO-invisible.", "",
        f"- switch pairs scored: **{len(pair_df)}** "
        f"(GO-invisible {int(pair_df['go_invisible'].sum())}, "
        f"GO-visible {int((~pair_df['go_invisible']).sum())})", "",
        "| stratum | consequence | obs rate | null mean | enrichment | p |",
        "|---------|-------------|----------|-----------|------------|---|",
    ]
    for r in enrich.sort_values(["stratum", "consequence"]).itertuples():
        lines.append(f"| {r.stratum} | {r.consequence} | {r.obs_rate:.3f} | "
                     f"{r.null_mean:.3f} | {r.enrichment:.2f} | {r.p_emp:.3g} |")
    (out_dir / "SWITCH_CONSEQUENCE.md").write_text("\n".join(lines) + "\n")


def _default_artifact_dir(region: str, resolution: str) -> Path:
    sub = "brainseq" if region in {"caudate", "caudate_sczd", "hippocampus", "dlpfc"} else "gtex"
    return rel("real_data", sub, region, "_m", resolution)


def main() -> None:
    p = argparse.ArgumentParser(
        description="Within-gene permutation test of switch coding-consequence enrichment.")
    p.add_argument("--region", required=True,
                   help="brain region dir (e.g. caudate_sczd, caudate, cortex).")
    p.add_argument("--resolution", default="isograph_vae",
                   help="artifact subdir (isograph_vae or isograph_vae_res2).")
    p.add_argument("--artifact-dir", default=None,
                   help="override the artifact dir; else derived from region+resolution.")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--seed", type=int, default=13)
    args = p.parse_args()
    artifact_dir = (Path(args.artifact_dir) if args.artifact_dir
                    else _default_artifact_dir(args.region, args.resolution))
    if not (artifact_dir / "module_interpret" / "structure_switch_pairs.parquet").exists():
        raise SystemExit(f"no structure_switch_pairs under {artifact_dir}; "
                         f"run interpret_modules first.")
    run(args.region, artifact_dir, args.fdr, args.n_perm, args.seed)


if __name__ == "__main__":
    main()
