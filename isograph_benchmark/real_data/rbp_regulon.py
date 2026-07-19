"""Stage 2: are IsoGraph co-switch modules RBP regulons?

The coding-consequence analysis showed the switch axis is a UTR-remodeling layer (UTR change
enriched 1.27x, 10/10 regions). 3'UTR usage is heavily RBP-controlled, which motivates the
coordination hypothesis: a co-switch module is coordinated because its member genes share a
common trans-acting RBP — the one thing cis-sQTL anchoring does not explain.

Using the per-(transcript, RBP) motif counts from `rbp_scan.py`, for each switch gene we call
whether the switch GAINS or LOSES an RBP binding site (motif present in one switch-pair isoform,
absent in the other). Then, per module, we test whether a given RBP's site-switching is
over-represented among the module's genes relative to the pooled switch-gene background
(hypergeometric) — a significant RBP marks the module as a candidate regulon for it. The
within-pair "gained/lost" call compares two isoforms of the SAME gene, controlling transcript
length/composition. Results are stratified GO-invisible vs GO-visible.

Output (real_data/_m/rbp/): rbp_switch_calls.parquet (per gene x RBP), rbp_regulon.parquet
(per module x RBP enrichment), RBP_REGULON.md, and a cross-region meta.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from isograph_benchmark.paths import rel
from isograph_benchmark.real_data.coloc_prep import load_switch_genes
from isograph_benchmark.real_data.qtl_anchoring import _bare

_RBP_DIR = rel("real_data", "_m", "rbp")
_COUNTS = _RBP_DIR / "rbp_counts.parquet"
_TREE_OF = {**{r: "brainseq" for r in ("caudate_sczd", "caudate", "hippocampus", "dlpfc")}}
_REGIONS = [
    ("brainseq", "caudate_sczd"), ("brainseq", "caudate"), ("brainseq", "hippocampus"),
    ("brainseq", "dlpfc"),
    *[("gtex", r) for r in (
        "amygdala", "anterior_cingulate_cortex_ba24", "caudate_basal_ganglia",
        "cerebellar_hemisphere", "cerebellum", "cortex", "frontal_cortex_ba9",
        "hippocampus", "hypothalamus", "nucleus_accumbens_basal_ganglia",
        "putamen_basal_ganglia", "spinal_cord_cervical_c_1", "substantia_nigra")],
]


def _presence(counts: pd.DataFrame) -> dict[tuple[str, str], int]:
    """(transcript_id, rbp) -> hit count (missing = 0)."""
    return {(t, r): int(c) for t, r, c in
            counts[["transcript_id", "rbp", "count"]].itertuples(index=False)}


def _switch_calls(region_tree: str, region: str, pres: dict, rbps: list[str],
                  fdr: float) -> pd.DataFrame:
    """Per (gene, RBP): does the switch gain/lose the motif (present in one isoform only)?"""
    art = rel("real_data", region_tree, region, "_m", "isograph_vae")
    sp_path = art / "module_interpret" / "structure_switch_pairs.parquet"
    if not sp_path.exists():
        return pd.DataFrame()
    sg = load_switch_genes(art, fdr)
    if sg.empty:
        return pd.DataFrame()
    sp = pd.read_parquet(sp_path)
    sp["gene"] = _bare(sp["gene_id"])
    tag = sg.groupby("gene").agg(go_invisible=("go_invisible", "max"),
                                 module_id=("module_id", "first")).reset_index()
    sp = sp.merge(tag, on="gene", how="inner")
    if sp.empty:
        return pd.DataFrame()

    rows = []
    for gene, sub in sp.groupby("gene"):
        go_inv = bool(sub["go_invisible"].iloc[0])
        mod = sub["module_id"].iloc[0]
        for rbp in rbps:
            # gained/lost across ANY of the gene's switch pairs
            switched = False
            for r in sub.itertuples():
                c1 = pres.get((r.transcript_id_1, rbp), 0)
                c2 = pres.get((r.transcript_id_2, rbp), 0)
                if (c1 > 0) != (c2 > 0):
                    switched = True
                    break
            rows.append((region, gene, mod, go_inv, rbp, switched))
    return pd.DataFrame(rows, columns=["region", "gene", "module_id", "go_invisible",
                                       "rbp", "switched"])


def _regulon_enrich(calls: pd.DataFrame) -> pd.DataFrame:
    """Per (region, module, RBP): hypergeometric over-representation of RBP site-switching
    among the module's genes vs the region's switch-gene pool."""
    out = []
    for region, rc in calls.groupby("region"):
        genes = rc["gene"].unique()
        N = len(genes)
        # pool: number of genes with this RBP switched (region-wide)
        pool = rc.groupby("rbp")["switched"].sum().to_dict()
        gene_mod = rc.drop_duplicates("gene").set_index("gene")["module_id"]
        gene_inv = rc.drop_duplicates("gene").set_index("gene")["go_invisible"]
        for module, mg in rc.groupby("module_id"):
            module_genes = mg["gene"].unique()
            n = len(module_genes)
            if n < 3:
                continue
            go_inv = bool(gene_inv.reindex(module_genes).mode().iloc[0]) \
                if len(module_genes) else False
            per_rbp = mg.groupby("rbp")["switched"].sum()
            for rbp, k in per_rbp.items():
                K = int(pool.get(rbp, 0))
                if K == 0 or k == 0:
                    continue
                p = float(hypergeom.sf(int(k) - 1, N, K, n))
                exp = n * K / N
                out.append({
                    "region": region, "module_id": module, "go_invisible": go_inv,
                    "rbp": rbp, "module_size": n, "n_switched": int(k),
                    "pool_switched": K, "pool_size": N,
                    "expected": exp, "enrichment": (k / exp if exp > 0 else np.nan),
                    "p": p})
    res = pd.DataFrame(out)
    if not res.empty:
        res["q"] = multipletests(res["p"], method="fdr_bh")[1]
    return res


def run(fdr: float) -> None:
    if not _COUNTS.exists():
        raise SystemExit(f"{_COUNTS} missing; run rbp_scan.py (motif env) first.")
    counts = pd.read_parquet(_COUNTS)
    pres = _presence(counts)
    rbps = sorted(counts["rbp"].unique())
    print(f"loaded motif counts: {len(rbps)} RBPs, {counts['transcript_id'].nunique():,} transcripts")

    calls = []
    for tree, region in _REGIONS:
        c = _switch_calls(tree, region, pres, rbps, fdr)
        if not c.empty:
            calls.append(c)
            print(f"  {region}: {c['gene'].nunique()} switch genes")
    calls = pd.concat(calls, ignore_index=True)
    _RBP_DIR.mkdir(parents=True, exist_ok=True)
    calls.to_parquet(_RBP_DIR / "rbp_switch_calls.parquet", index=False)

    reg = _regulon_enrich(calls)
    reg.to_parquet(_RBP_DIR / "rbp_regulon.parquet", index=False)
    _write_report(reg, calls)
    n_sig = int((reg["q"] < 0.05).sum()) if not reg.empty else 0
    print(f"candidate RBP regulons (q<0.05): {n_sig} module-RBP pairs across "
          f"{reg['region'].nunique() if not reg.empty else 0} regions")


def _write_report(reg: pd.DataFrame, calls: pd.DataFrame) -> None:
    sig = reg[reg["q"] < 0.05].sort_values("q") if not reg.empty else reg
    lines = [
        "# Candidate RBP regulons among IsoGraph co-switch modules", "",
        "For each module and RBP, whether RBP binding-site switching (motif gained/lost "
        "between the switch-pair isoforms) is over-represented among the module's genes vs "
        "the region's switch-gene pool (hypergeometric, BH across module x RBP tests). A "
        "significant RBP marks the module as a candidate regulon for it.", "",
        f"- switch genes tested: **{calls['gene'].nunique()}** over "
        f"{calls['region'].nunique()} regions",
        f"- module x RBP tests: **{len(reg)}**; candidate regulons at q<0.05: "
        f"**{int((reg['q'] < 0.05).sum()) if not reg.empty else 0}** "
        f"(GO-invisible {int(((reg['q'] < 0.05) & reg['go_invisible']).sum()) if not reg.empty else 0})",
        "",
        "| region | module | RBP | module size | switched | enrichment | q | GO-inv |",
        "|--------|--------|-----|-------------|----------|------------|---|--------|",
    ]
    for r in sig.head(30).itertuples():
        lines.append(f"| {r.region} | {r.module_id} | {r.rbp} | {r.module_size} | "
                     f"{r.n_switched} | {r.enrichment:.2f} | {r.q:.2e} | "
                     f"{'yes' if r.go_invisible else 'no'} |")
    (_RBP_DIR / "RBP_REGULON.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="Per-module RBP-regulon enrichment (stage 2).")
    p.add_argument("--fdr", type=float, default=0.05)
    args = p.parse_args()
    run(args.fdr)


if __name__ == "__main__":
    main()
