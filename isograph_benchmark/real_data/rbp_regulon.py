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
from isograph_benchmark.real_data.qtl_anchoring import _bare, build_gene_sets

_RBP_DIR = rel("real_data", "_m", "rbp")
_COUNTS = _RBP_DIR / "rbp_counts.parquet"
_FLANK_NOTE = 100        # intronic flank window (nt); mirrors rbp_scan_intronic._FLANK
# per-scope Stage-1 count tables; "combined" unions the mature + intronic presence
_SCOPE_COUNTS = {
    "mature": [_RBP_DIR / "rbp_counts.parquet"],
    "intronic": [_RBP_DIR / "rbp_counts_intronic.parquet"],
    "combined": [_RBP_DIR / "rbp_counts.parquet", _RBP_DIR / "rbp_counts_intronic.parquet"],
}
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


def _gene_tags(art, fdr: float) -> tuple[pd.DataFrame, str]:
    """(gene, module_id, go_invisible) tags + provenance.

    Primary pool = phenotype-associated switch genes (`load_switch_genes`). When that is empty
    (regions with no FDR-significant phenotype/age association, e.g. several GTEx tissues), fall
    back to the region's FULL module-gene pool so the RBP-regulon question — "is this co-switch
    module coordinated by an RBP?" — is still asked; the GO-invisible tag is preserved from
    `build_gene_sets`, and a `pool_source` flag records which universe was used.
    """
    sg = load_switch_genes(art, fdr)
    if not sg.empty:
        return sg[["gene", "module_id", "go_invisible"]].drop_duplicates(), "switch_genes"
    enrich_path = art.parent / "module_enrichment" / "isograph_modules.parquet"
    inv = build_gene_sets(art, enrich_path, fdr).get("go_invisible_modules", set())
    mods = pd.read_parquet(art / "modules.parquet")
    mods["gene"] = _bare(mods["gene_id"])
    mods["module_id"] = mods["module_id"].astype(str)
    mods["go_invisible"] = mods["gene"].isin(inv)
    return mods[["gene", "module_id", "go_invisible"]].drop_duplicates(), "module_genes"


def _switch_calls(region_tree: str, region: str, pres: dict, rbps: list[str],
                  fdr: float) -> pd.DataFrame:
    """Per (gene, RBP): does the switch gain/lose the motif (present in one isoform only)?"""
    art = rel("real_data", region_tree, region, "_m", "isograph_vae")
    sp_path = art / "module_interpret" / "structure_switch_pairs.parquet"
    if not sp_path.exists():
        return pd.DataFrame()
    tag, pool_source = _gene_tags(art, fdr)
    if tag.empty:
        return pd.DataFrame()
    sp = pd.read_parquet(sp_path)
    sp["gene"] = _bare(sp["gene_id"])
    tag = tag.groupby("gene").agg(go_invisible=("go_invisible", "max"),
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
            rows.append((region, gene, mod, go_inv, rbp, switched, pool_source))
    return pd.DataFrame(rows, columns=["region", "gene", "module_id", "go_invisible",
                                       "rbp", "switched", "pool_source"])


def _regulon_enrich(calls: pd.DataFrame) -> pd.DataFrame:
    """Per (region, module, RBP): hypergeometric over-representation of RBP site-switching
    among the module's genes vs the region's switch-gene pool."""
    out = []
    for region, rc in calls.groupby("region"):
        genes = rc["gene"].unique()
        N = len(genes)
        pool_source = rc["pool_source"].iloc[0] if "pool_source" in rc.columns else "switch_genes"
        # pool: number of genes with this RBP switched (region-wide)
        pool = rc.groupby("rbp")["switched"].sum().to_dict()
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
                    "p": p, "pool_source": pool_source})
    res = pd.DataFrame(out)
    if not res.empty:
        res["q"] = multipletests(res["p"], method="fdr_bh")[1]
    return res


def run(fdr: float, scope: str = "mature") -> None:
    paths = _SCOPE_COUNTS[scope]
    missing = [p for p in paths if not p.exists()]
    if missing:
        hint = "rbp_scan.py" if scope == "mature" else "rbp_scan_intronic.py"
        raise SystemExit(f"{missing[0]} missing; run {hint} (motif env) first "
                         f"(scope={scope} needs {[p.name for p in paths]}).")
    counts = pd.concat([pd.read_parquet(p) for p in paths], ignore_index=True)
    if len(paths) > 1:                     # combined: union presence across scopes
        counts = counts.groupby(["transcript_id", "rbp"], as_index=False)["count"].sum()
    pres = _presence(counts)
    rbps = sorted(counts["rbp"].unique())
    print(f"[{scope}] loaded motif counts: {len(rbps)} RBPs, "
          f"{counts['transcript_id'].nunique():,} transcripts")

    calls = []
    for tree, region in _REGIONS:
        c = _switch_calls(tree, region, pres, rbps, fdr)
        if not c.empty:
            calls.append(c)
            print(f"  {region}: {c['gene'].nunique()} switch genes")
    calls = pd.concat(calls, ignore_index=True)
    _RBP_DIR.mkdir(parents=True, exist_ok=True)
    suffix = "" if scope == "mature" else f"_{scope}"
    calls.to_parquet(_RBP_DIR / f"rbp_switch_calls{suffix}.parquet", index=False)

    reg = _regulon_enrich(calls)
    reg.to_parquet(_RBP_DIR / f"rbp_regulon{suffix}.parquet", index=False)
    _write_report(reg, calls, _RBP_DIR / f"RBP_REGULON{suffix}.md", scope)
    n_sig = int((reg["q"] < 0.05).sum()) if not reg.empty else 0
    print(f"[{scope}] candidate RBP regulons (q<0.05): {n_sig} module-RBP pairs across "
          f"{reg['region'].nunique() if not reg.empty else 0} regions")


def _write_report(reg: pd.DataFrame, calls: pd.DataFrame,
                  out_path: Path | None = None, scope: str = "mature") -> None:
    out_path = out_path or (_RBP_DIR / "RBP_REGULON.md")
    scope_note = {
        "mature": "Motifs are scanned on the **mature transcript** (exonic + UTR) sequence.",
        "intronic": "Motifs are scanned on **intronic splice-site flanks** "
                    f"(pre-mRNA sense; up to {_FLANK_NOTE} nt into each intron), the binding "
                    "niche for splicing-regulatory RBPs invisible to the mature-transcript scan.",
        "combined": "Motif presence unions the **mature-transcript** and **intronic "
                    "splice-site flank** scans.",
    }[scope]
    sig = reg[reg["q"] < 0.05].sort_values("q") if not reg.empty else reg
    lines = [
        f"# Candidate RBP regulons among IsoGraph co-switch modules ({scope} scope)", "",
        scope_note, "",
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
    out_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="Per-module RBP-regulon enrichment (stage 2).")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--scope", choices=("mature", "intronic", "combined"), default="mature",
                   help="which Stage-1 count table(s) to consume (default: mature, canonical)")
    args = p.parse_args()
    run(args.fdr, args.scope)


if __name__ == "__main__":
    main()
