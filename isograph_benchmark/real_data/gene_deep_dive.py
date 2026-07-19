"""Per-gene mechanistic deep-dive for the colocalization-prioritized disease genes.

For each target gene this joins the genetic-anchoring and downstream-biology layers we have
already built into a single mechanistic vignette:

  1. genetic anchor      -- coloc_isoform_events_combined / coloc_direction_combined:
                            which trait, sQTL vs eQTL, tissue, CLPP, risk allele + rsID,
                            signed direction of the risk allele on the QTL.
  2. the switch          -- whether the colocalized junction maps onto an IsoGraph switch
                            pair (junction_in_switch_pair), the switch pair, GTEx concordance,
                            GO-invisibility, and independent BrainSeq replication.
  3. coding consequence  -- the resolved event's structural_consequence.
  4. regulatory logic    -- rbp_switch_calls / rbp_regulon: which RBP motifs are called
                            switched in the gene, and which of those are also enriched in the
                            gene's IsoGraph module (candidate splice regulators of the switch).
  5. constraint / clinic -- gnomAD LOEUF + missense o/e from the clinical-consequence layer.

Each gene is classified as splicing-led (an sQTL that colocalizes onto a concordant IsoGraph
switch pair -- the IsoGraph-unique, GO-invisible case), expression-led (eQTL gene-level only),
or splicing-unresolved (an sQTL that does not map onto the switch pair). Writes one markdown
vignette per gene plus a panel-wide parquet of one row per gene.

Deterministic; pure joins over existing parquets (no heavy compute). Usage:
  python -m isograph_benchmark.real_data.gene_deep_dive
  python -m isograph_benchmark.real_data.gene_deep_dive --genes SNCA,CTSH,PPP6R2
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel

# Originally hand-prioritized panel (genetic_anchoring_v2), kept for reference. The CLI
# defaults to EVERY colocalized gene so multi-locus / multi-trait cases are not dropped;
# main-figure framing is reserved for the resolved cross-disease headliners (SNCA, CTSH).
PANEL = [
    "SNCA", "TPP1", "SCFD1", "PGS1", "PPP6R2", "GGNBP2",
    "CTSH", "TPCN1", "MRPS10", "MYO18A", "PBX1", "RBFA", "MED15",
]
_Q_ENRICH = 0.05


def _all_coloc_genes() -> list[str]:
    ev = pd.read_parquet(rel("real_data", "coloc", "_m", "coloc_isoform_events_combined.parquet"))
    return sorted(ev["gene_name"].dropna().unique())


def _ens(series: pd.Series) -> pd.Series:
    return series.astype(str).str.split(".").str[0]


def _load() -> dict:
    ev = pd.read_parquet(rel("real_data", "coloc", "_m", "coloc_isoform_events_combined.parquet"))
    di = pd.read_parquet(rel("real_data", "coloc", "_m", "coloc_direction_combined.parquet"))
    calls = pd.read_parquet(rel("real_data", "_m", "rbp", "rbp_switch_calls.parquet"))
    regulon = pd.read_parquet(rel("real_data", "_m", "rbp", "rbp_regulon.parquet"))
    con = []
    for tree in ("brainseq", "gtex"):
        for f in rel("real_data", tree).glob(
                "*/_m/isograph_vae/clinical_consequence/gene_constraint.parquet"):
            con.append(pd.read_parquet(f))
    constraint = pd.concat(con, ignore_index=True) if con else pd.DataFrame(
        columns=["gene", "loeuf", "mis_oe", "region"])
    for df in (ev, di, calls, constraint):
        if "gene" in df.columns:
            df["ens"] = _ens(df["gene"])
    return {"ev": ev, "di": di, "calls": calls, "regulon": regulon, "constraint": constraint}


def _rbp_layer(gene_calls: pd.DataFrame, regulon: pd.DataFrame) -> dict:
    """Switched RBP motifs for the gene, flagged by whether they are also module-enriched."""
    switched = gene_calls[gene_calls["switched"]]
    if switched.empty:
        return {"switched_rbps": [], "module_enriched_rbps": []}
    # recurrence of each switched RBP across the gene's region/module memberships
    rec = switched.groupby("rbp")["region"].nunique().sort_values(ascending=False)
    # module-level enrichment: match the gene's (region, module_id) to the regulon
    mods = switched[["region", "module_id"]].drop_duplicates()
    reg = regulon.merge(mods, on=["region", "module_id"], how="inner")
    reg = reg[reg["q"] < _Q_ENRICH]
    # keep only RBPs that are BOTH switched in the gene AND enriched in its module
    reg = reg[reg["rbp"].isin(set(switched["rbp"]))].sort_values("q")
    enr = (reg.groupby("rbp")["q"].min().sort_values().index.tolist())
    return {
        "switched_rbps": [f"{r}({n})" for r, n in rec.items()],
        "module_enriched_rbps": enr,
    }


def _classify(g_ev: pd.DataFrame) -> tuple[str, bool]:
    """Return (verdict, go_invisible_resolved)."""
    sqtl = g_ev[g_ev["kind"] == "sQTL"]
    resolved = sqtl[(sqtl["junction_in_switch_pair"]) & (sqtl["concordant"] == True)]  # noqa: E712
    if not resolved.empty:
        go_inv = bool(resolved["go_invisible"].any())
        return "splicing-led (IsoGraph-resolved switch)", go_inv
    if (g_ev["kind"] == "eQTL").all():
        return "expression-led (eQTL gene-level)", bool(g_ev["go_invisible"].any())
    if not sqtl.empty:
        return "splicing (sQTL not resolved to switch pair)", bool(g_ev["go_invisible"].any())
    return "unclassified", False


def _gene_row(gene: str, d: dict) -> dict | None:
    ev = d["ev"][d["ev"]["gene_name"] == gene]
    if ev.empty:
        return None
    ens = ev["ens"].iloc[0]
    di = d["di"][d["di"]["ens"] == ens]
    calls = d["calls"][d["calls"]["ens"] == ens]
    con = d["constraint"][d["constraint"]["ens"] == ens]
    rbp = _rbp_layer(calls, d["regulon"])
    verdict, go_inv = _classify(ev)
    resolved = ev[(ev["kind"] == "sQTL") & (ev["junction_in_switch_pair"]) &
                  (ev["concordant"] == True)]  # noqa: E712
    top = ev.sort_values("clpp", ascending=False).iloc[0]
    return {
        "gene": gene,
        "ens": ens,
        "traits": ",".join(sorted(ev["trait"].str.upper().unique())),
        "kinds": ",".join(sorted(ev["kind"].unique())),
        "n_events": int(len(ev)),
        "max_clpp": float(ev["clpp"].max()),
        "top_tissue": str(top["tissue"]),
        "top_rsid": str(top["best_rsid"]),
        "top_risk_allele": str(top["risk_allele"]),
        "resolved_to_switch_pair": bool(not resolved.empty),
        "n_resolved_events": int(len(resolved)),
        "multi_locus": bool(len(resolved) >= 2),
        "concordant_traits": ",".join(sorted(resolved["trait"].str.upper().unique()))
                             if not resolved.empty else "",
        "brainseq_replicates": bool(ev["brainseq_replicates_switch"].eq(True).any()),
        "go_invisible": bool(go_inv),
        "loeuf": float(con["loeuf"].min()) if not con.empty else np.nan,
        "mis_oe": float(con["mis_oe"].min()) if not con.empty else np.nan,
        "n_switched_rbps": len(rbp["switched_rbps"]),
        "module_enriched_rbps": ",".join(rbp["module_enriched_rbps"][:10]),
        "verdict": verdict,
        "_ev": ev, "_di": di, "_rbp": rbp, "_con": con, "_resolved": resolved,
    }


def _vignette(row: dict) -> str:
    ev, di, rbp = row["_ev"], row["_di"], row["_rbp"]
    L = [f"# {row['gene']} — mechanistic deep-dive", ""]
    L.append(f"**Verdict:** {row['verdict']}"
             f"{' · GO-invisible' if row['go_invisible'] else ''}"
             f"{' · replicates in BrainSeq' if row['brainseq_replicates'] else ''}")
    L.append("")
    L.append(f"- **Traits:** {row['traits']}  ·  **QTL kinds:** {row['kinds']}  ·  "
             f"**max CLPP:** {row['max_clpp']:.2f} ({row['top_tissue']})")
    loeuf = "n/a" if np.isnan(row["loeuf"]) else f"{row['loeuf']:.3f}"
    L.append(f"- **Constraint:** LOEUF {loeuf}  ·  missense o/e "
             f"{'n/a' if np.isnan(row['mis_oe']) else f'{row['mis_oe']:.2f}'}")
    L.append("")
    L.append("## 1–3. Genetic anchor → switch → coding consequence")
    L.append("| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |")
    L.append("|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|")
    for r in ev.sort_values(["kind", "clpp"], ascending=[True, False]).itertuples():
        sc = (str(r.structural_consequence) or "").strip() or "—"
        L.append(f"| {r.trait} | {r.kind} | {r.tissue} | {r.clpp:.2f} | "
                 f"{'yes' if r.junction_in_switch_pair else 'no'} | "
                 f"{'yes' if r.concordant == True else ('no' if r.concordant == False else '—')} | "  # noqa: E712
                 f"{'yes' if r.go_invisible else 'no'} | {sc[:48]} | {str(r.resolved_event)[:90]} |")
    if not di.empty:
        L.append("")
        L.append("### Signed risk-allele direction (colocalized loci with allele matching)")
        L.append("| trait | tissue | rsID | risk allele | risk QTL effect | direction |")
        L.append("|-------|--------|------|-------------|-----------------|-----------|")
        for r in di.itertuples():
            L.append(f"| {r.trait} | {r.tissue} | {r.best_rsid} | {r.risk_allele} | "
                     f"{r.risk_qtl_effect} | {r.direction} |")
    L.append("")
    L.append("## 4. Regulatory logic (RBP motifs in switched exons)")
    if rbp["module_enriched_rbps"]:
        L.append(f"- **Switched *and* module-enriched (q<{_Q_ENRICH}) RBPs:** "
                 f"{', '.join(rbp['module_enriched_rbps'][:15])}")
    L.append(f"- **All switched-motif RBPs (recurrence across regions):** "
             f"{', '.join(rbp['switched_rbps'][:20]) or '—'}")
    L.append("")
    L.append("## 5. Interpretation")
    L.append(_interpretation(row))
    L.append("")
    return "\n".join(L)


def _interpretation(row: dict) -> str:
    if row["resolved_to_switch_pair"]:
        s = (f"A splicing QTL colocalizes onto an IsoGraph switch pair for {row['gene']} "
             f"in {row['concordant_traits']}")
        if "," in row["concordant_traits"]:
            s += " — the same switch is genetically anchored across more than one trait"
        s += (". This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts "
              "through isoform choice, not gene dosage")
        if row["go_invisible"]:
            s += ", in a GO-invisible module a pathway-enrichment scan would miss"
        if row["brainseq_replicates"]:
            s += "; the switch also replicates in an independent BrainSeq cohort"
        return s + "."
    if "expression-led" in row["verdict"]:
        return (f"{row['gene']} colocalizes as an eQTL (gene-level expression), with no splicing "
                "event resolving to an IsoGraph switch pair — an honest expression-confounded case "
                "that abundance networks would also capture.")
    return (f"{row['gene']} has a colocalizing sQTL, but the junction does not map onto the "
            "IsoGraph switch pair for the tissue — splicing-associated but not resolved to a "
            "switch; a candidate for deeper transcript-level follow-up.")


def run(genes: list[str]) -> pd.DataFrame:
    d = _load()
    out_dir = ensure_dir(rel("real_data", "_m", "deep_dive"))
    rows = []
    for g in genes:
        row = _gene_row(g, d)
        if row is None:
            print(f"  {g}: no colocalized isoform events; skipping.")
            continue
        (out_dir / f"{g}.md").write_text(_vignette(row))
        rows.append({k: v for k, v in row.items() if not k.startswith("_")})
    panel = pd.DataFrame(rows)
    # order: resolved splicing-led first, multi-locus above single, then by CLPP
    panel = panel.sort_values(["resolved_to_switch_pair", "n_resolved_events", "max_clpp"],
                              ascending=[False, False, False]).reset_index(drop=True)
    panel.to_parquet(out_dir / "deep_dive_panel.parquet", index=False)
    _write_panel_md(panel, out_dir)
    _write_supp_tables(d, set(panel["ens"]), out_dir)
    print(f"deep-dive over {len(panel)} genes -> {out_dir}")
    return panel


def _emit(df: pd.DataFrame, out_dir, name: str) -> None:
    """Write a supplementary table as both parquet (analysis) and tsv (reader-facing)."""
    df.to_parquet(out_dir / f"{name}.parquet", index=False)
    df.to_csv(out_dir / f"{name}.tsv", sep="\t", index=False)
    print(f"  supp table {name}: {len(df)} rows")


def _write_supp_tables(d: dict, ens_set: set, out_dir) -> None:
    """Machine-readable per-gene tables so readers can reconstruct any gene's deep-dive.

    (1) events  - one row per colocalized isoform event (anchor -> switch -> consequence),
                  merged with the signed risk-allele direction where allele matching succeeded;
    (2) rbp     - per gene, RBP motifs both switched in the gene and enriched in its module;
    (3) exons   - per gene/region/exon: switched vs constitutive, CDS overlap, ClinVar P/LP.
    """
    # (1) per-event table + signed direction
    ev = d["ev"][d["ev"]["ens"].isin(ens_set)].copy()
    # one direction row per locus (the table is keyed per-junction; dedup avoids fan-out)
    dkeep = (d["di"][["ens", "trait", "tissue", "best_rsid", "variant_id", "ref", "alt",
                      "risk_beta", "slope", "direction"]]
             .drop_duplicates(["ens", "trait", "tissue", "best_rsid"]))
    events = ev.merge(dkeep, on=["ens", "trait", "tissue", "best_rsid"], how="left")
    cols = ["gene_name", "ens", "trait", "case", "kind", "tissue", "best_rsid", "risk_allele",
            "risk_qtl_effect", "direction", "variant_id", "ref", "alt", "junction", "clpp",
            "go_invisible", "switch_pair", "junction_in_switch_pair", "structural_consequence",
            "concordant", "brainseq_region", "brainseq_replicates_switch", "resolved_event"]
    events = events[[c for c in cols if c in events.columns]].sort_values(
        ["gene_name", "trait", "clpp"], ascending=[True, True, False])
    _emit(events, out_dir, "deep_dive_events")

    # (2) per-gene switched + module-enriched RBP regulators
    ens2sym = dict(zip(d["ev"]["ens"], d["ev"]["gene_name"]))
    calls = d["calls"][(d["calls"]["ens"].isin(ens_set)) & (d["calls"]["switched"])]
    reg = d["regulon"][d["regulon"]["q"] < _Q_ENRICH][
        ["region", "module_id", "rbp", "enrichment", "q"]]
    rbp = calls.merge(reg, on=["region", "module_id", "rbp"], how="inner")
    rbp["gene_name"] = rbp["ens"].map(ens2sym)
    rbp = rbp[["gene_name", "ens", "region", "module_id", "go_invisible", "rbp",
               "enrichment", "q"]].sort_values(["gene_name", "q"])
    _emit(rbp, out_dir, "deep_dive_rbp")

    # (3) per-gene/exon clinical annotation (SNCA-style read for every gene)
    frames = []
    for tree in ("brainseq", "gtex"):
        for f in rel("real_data", tree).glob(
                "*/_m/isograph_vae/clinical_consequence/exon_clinvar.parquet"):
            x = pd.read_parquet(f)
            x["ens"] = _ens(x["gene"])
            x = x[x["ens"].isin(ens_set)]
            if not x.empty:
                x["region"] = f.parents[2].name if "region" not in x.columns else x["region"]
                frames.append(x)
    if frames:
        exons = pd.concat(frames, ignore_index=True)
        exons["gene_name"] = exons["ens"].map(ens2sym)
        keep = ["gene_name", "ens", "region", "chrom", "start", "end", "length", "switched",
                "cds_overlap", "n_plp", "n_clinvar", "go_invisible"]
        exons = exons[[c for c in keep if c in exons.columns]].sort_values(
            ["gene_name", "region", "start"])
        _emit(exons, out_dir, "deep_dive_exon_clinical")


def _write_panel_md(panel: pd.DataFrame, out_dir) -> None:
    L = ["# Per-gene mechanistic deep-dive — panel summary", "",
         "Every colocalized disease gene, one row each (ranked splicing-led first, multi-locus "
         "above single). `resolved events` = colocalizing sQTLs that map onto a concordant "
         "IsoGraph switch pair (the splicing-led, DTU-without-DGE class); `multi-locus` flags "
         "genes with >=2 such events (multiple significant colocalizations, incl. cross-trait). "
         "Main-figure framing is reserved for the resolved cross-disease headliners (SNCA, "
         "CTSH); the remainder are supporting vignettes. See `<GENE>.md` for each.", "",
         "| gene | traits | kinds | max CLPP | LOEUF | resolved events | multi-locus | concordant traits | BrainSeq rep | GO-inv | verdict |",
         "|------|--------|-------|----------|-------|-----------------|-------------|-------------------|--------------|--------|---------|"]
    for r in panel.itertuples():
        loeuf = "n/a" if pd.isna(r.loeuf) else f"{r.loeuf:.2f}"
        L.append(f"| {r.gene} | {r.traits} | {r.kinds} | {r.max_clpp:.2f} | {loeuf} | "
                 f"{r.n_resolved_events} | {'yes' if r.multi_locus else '—'} | "
                 f"{r.concordant_traits or '—'} | "
                 f"{'yes' if r.brainseq_replicates else 'no'} | "
                 f"{'yes' if r.go_invisible else 'no'} | {r.verdict} |")
    (out_dir / "DEEP_DIVE_PANEL.md").write_text("\n".join(L) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description="Per-gene mechanistic deep-dive of coloc genes.")
    ap.add_argument("--genes", default="",
                    help="comma-separated gene symbols (default: every colocalized gene).")
    args = ap.parse_args()
    genes = [g.strip() for g in args.genes.split(",") if g.strip()] or _all_coloc_genes()
    run(genes)


if __name__ == "__main__":
    main()
