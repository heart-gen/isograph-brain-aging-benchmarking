"""Nominate a wet-lab perturbation panel for a candidate splicing regulator (default KHDRBS1).

The RBP-regulon stage (`rbp_regulon.py`) answers a module-level question -- *which IsoGraph
co-switch modules are over-represented for site-switching of a given RBP*. A knockdown /
knockout experiment needs the gene-level question instead: *which individual genes actually
gain or lose that RBP's motif between their IsoGraph switch-pair isoforms*, inside a module
that was independently nominated as that RBP's regulon. This CLI does exactly that join and
freezes the resulting target list, so the panel sent to a wet-lab collaborator is derived from
the pipeline rather than picked by hand from the manuscript text.

Eligibility (the candidate universe) is the conjunction of the two layers:

  1. `rbp_switch_calls`  -- gene x RBP, `switched == True`: the gene's switch pair changes
                            motif presence (present in one isoform, absent in the other).
  2. `rbp_regulon`       -- (region, module) x RBP, hypergeometric `q <= --q`: the gene's
                            module is itself an enriched regulon for the RBP.

Genes are then ranked deterministically. The default `--rank-by adjusted` leads on the
covariate-adjusted regulon GLM (`rbp_regulon.odds_ratio`, which adjusts for length, GC,
UTR-vs-CDS composition and transcript count): genes whose best module clears
`--min-adjusted-or` come first, ordered by the lower bound of that module's adjusted 95% CI,
so a module whose adjusted interval covers 1 cannot outrank one whose interval excludes it.
Recurrence across brain regions, eCLIP-binding support (`rbp_binding_regulon`), unadjusted
module q, the gene's own switch-driver strength (max |switch_r| from `module_gene_roles`) and
finally symbol break the remaining ties. `--rank-by hypergeometric` reproduces the
pre-adjustment ordering (recurrence first, then eCLIP, then unadjusted q and enrichment) so
the two panels can be diffed.

Eligibility deliberately stays on the hypergeometric q even under `--rank-by adjusted`. For
KHDRBS1 no module reaches adjusted q <= 0.05, so gating on the adjusted arm would empty the
panel; the adjusted arm is strong enough to order candidates and to demote modules whose
enrichment is an opportunity artefact, not to gate them. The report states this, prints the
best adjusted q, and flags every sampled module the adjustment contradicts. Genes are also
required to sit in a phenotype-associated module in at least `--min-pheno-regions` region(s)
(`pool_source == switch_genes`).

Recurrence alone concentrates the panel on whichever single module recurs most widely -- for
KHDRBS1 that is one inflammatory/immune module, whose top genes would all test the same
biology (and are largely microglial, not neuronal). `--max-per-module` therefore caps how many
panel genes may come from the same (region, module), so the panel samples the regulator's
program rather than its largest module (the cap applies to the gene's best -- most enriched --
module, the one it is reported under). Each selected gene carries the GO terms of that module,
so the biology being sampled is visible in the report.

The panel is tiered:

  Tier 0  regulator            -- the RBP itself (the perturbation target).
  Tier 1  core recurrent       -- top-ranked direct targets; the primary readout that the RBP
                                 drives a coordinated switching program.
  Tier 2  disease anchor       -- eligible genes that also carry a colocalized disease locus
                                 (`deep_dive_panel`), filled splicing-led (the sQTL resolves
                                 onto the IsoGraph switch pair) first, then sQTL-anchored but
                                 unresolved. Expression-led coloc is eQTL-driven and would
                                 confound a splicing readout, so it is excluded unless
                                 `--include-expression-led` (always flagged by `anchor_class`).
  Tier 3  GO-invisible program -- recurrent targets in GO-invisible modules with no coloc
                                 anchor: the program-level control, so the experiment is not
                                 reduced to validating GWAS loci.

Outputs (08_integration/_m/rbp_target_panel/<RBP>/):
  rbp_target_evidence.parquet/tsv     one row per eligible (gene, region)
  rbp_target_candidates.parquet/tsv   one row per eligible gene, ranked (the full universe)
  rbp_target_panel.parquet/tsv        the tiered panel
  rbp_target_switch_pairs.parquet/tsv the motif-differential switch pairs of the panel genes
  expression_check_list.tsv           symbol / Ensembl / tier -- the list to send the collaborator
  RBP_TARGET_PANEL.md                 report with the ranking rule, the panel, and the caveats

Deterministic; pure joins over existing parquets. Usage:
  python -m isograph_benchmark.real_data.rbp_target_panel
  python -m isograph_benchmark.real_data.rbp_target_panel --rbp QKI --n-core 10
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.real_data.partition_provenance import load_region_enrichment
from isograph_benchmark.paths import ensure_dir, region_store, stage_out
from isograph_benchmark.real_data.go_invisible_gate import gene_symbol_map
from isograph_benchmark.real_data.interpret_modules import DEFAULT_GTF_PATH
from isograph_benchmark.real_data.qtl_anchoring import _bare
from isograph_benchmark.real_data.rbp_regulon import (
    _REGIONS, _UNIT_COUNTS, _load_counts, _presence,
)

_RBP_DIR = stage_out("regulation", "rbp")
_OUT_ROOT = stage_out("integration", "rbp_target_panel")
_DEEP_DIVE = stage_out("integration", "deep_dive")
_TREE_OF = {region: tree for tree, region in _REGIONS}


_DEFAULT_RANK_BY = "adjusted"
_RANK_BY_CHOICES = ("adjusted", "hypergeometric")


def _scoped(base: str, scope: str, bg_mode: str) -> Path:
    """Stage-2 filename convention: `<base>[_<scope>][_flatbg].parquet`."""
    suffix = "" if scope == "mature" else f"_{scope}"
    suffix += "" if bg_mode == "composition" else "_flatbg"
    return _RBP_DIR / f"{base}{suffix}.parquet"


def _read(path: Path, what: str) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"{path} missing; run rbp_regulon.py first ({what}).")
    return pd.read_parquet(path)


# --------------------------------------------------------------------------- evidence layers
def _module_trust() -> pd.DataFrame:
    """(region, module_id) -> phenotype FDR / GO-term count / module size.

    Same table `load_switch_genes` / `build_gene_sets` use to define phenotype-associated and
    GO-invisible modules, so the trust columns here match the definitions used upstream.
    """
    rows = []
    for tree, region in _REGIONS:
        d = load_region_enrichment(
            tree, region, "isograph", context=f"rbp_target_panel {tree}/{region}"
        )
        if d is None:
            continue
        d = d[d["method"] == "isograph"] if "method" in d.columns else d
        keep = ["module_id", "n_genes", "n_go_terms", "pheno_fdr"]
        if "top_go_terms" in d.columns:
            keep.append("top_go_terms")
        d = d[keep].copy()
        d["module_id"] = d["module_id"].astype(str)
        d["region"] = region
        d["module_go"] = ([", ".join(list(t)[:3]) for t in d["top_go_terms"]]
                          if "top_go_terms" in d.columns else "")
        rows.append(d.drop(columns=[c for c in ("top_go_terms",) if c in d.columns]))
    if not rows:
        return pd.DataFrame(columns=["region", "module_id", "n_genes", "n_go_terms",
                                     "pheno_fdr", "module_go"])
    return pd.concat(rows, ignore_index=True).rename(columns={"n_genes": "module_n_genes"})


def _gene_roles() -> pd.DataFrame:
    """(region, gene) -> IsoGraph module role. switch_only genes are the cleanest readout:
    the module places them there on switching alone, with no abundance channel."""
    rows = []
    for tree, region in _REGIONS:
        p = region_store(tree, region, "isograph_vae", "module_gene_roles.parquet")
        if not p.exists():
            continue
        d = pd.read_parquet(p)
        d["gene"] = _bare(d["gene_id"])
        d["region"] = region
        rows.append(d[["region", "gene", "module_role", "switch_r", "switch_active"]])
    if not rows:
        return pd.DataFrame(columns=["region", "gene", "module_role", "switch_r", "switch_active"])
    return pd.concat(rows, ignore_index=True).drop_duplicates(["region", "gene"])


def _binding_support(rbp: str) -> pd.DataFrame:
    """(region, module_id) -> whether the regulon call is supported by ENCODE eCLIP binding
    that is preferential for the switched isoform (`rbp_binding.py`). Mature scope only."""
    p = _RBP_DIR / "rbp_binding_regulon.parquet"
    if not p.exists():
        return pd.DataFrame(columns=["region", "module_id", "binding_supported", "preferential"])
    d = pd.read_parquet(p)
    d = d[d["rbp"] == rbp]
    keep = ["region", "module_id", "binding_supported", "preferential", "mcnemar_p"]
    d = d[[c for c in keep if c in d.columns]].copy()
    d["module_id"] = d["module_id"].astype(str)
    return d


def _coloc_anchors() -> pd.DataFrame:
    """Bare Ensembl -> colocalized-disease-locus summary from the per-gene deep dive."""
    p = _DEEP_DIVE / "deep_dive_panel.parquet"
    if not p.exists():
        return pd.DataFrame(columns=["gene", "anchor_traits", "anchor_kinds", "anchor_max_clpp",
                                     "anchor_verdict", "anchor_splicing_led"])
    d = pd.read_parquet(p)
    out = pd.DataFrame({
        "gene": _bare(d["ens"]),
        "anchor_traits": d["traits"],
        "anchor_kinds": d["kinds"],
        "anchor_max_clpp": pd.to_numeric(d["max_clpp"], errors="coerce"),
        "anchor_verdict": d["verdict"],
    })
    v = out["anchor_verdict"].astype(str)
    # three coloc classes, in decreasing usefulness for a splicing readout: the sQTL resolves
    # onto the IsoGraph switch pair; an sQTL that does not resolve onto it; gene-level eQTL.
    out["anchor_class"] = np.where(
        v.str.startswith("splicing-led"), "splicing-led",
        np.where(v.str.startswith("splicing"), "splicing-unresolved",
                 np.where(v.str.startswith("expression-led"), "expression-led", "other")))
    out["anchor_splicing_led"] = out["anchor_class"] == "splicing-led"
    return out.drop_duplicates("gene")


# --------------------------------------------------------------------------- candidate build
def build_evidence(rbp: str, q: float, scope: str, bg_mode: str) -> pd.DataFrame:
    """One row per eligible (gene, region): motif switched AND module is an enriched regulon."""
    calls = _read(_scoped("rbp_switch_calls", scope, bg_mode), "stage-2 switch calls")
    regulon = _read(_scoped("rbp_regulon", scope, bg_mode), "stage-2 regulon enrichment")
    for name, d in (("switch calls", calls), ("regulon", regulon)):
        if rbp not in set(d["rbp"]):
            raise SystemExit(f"{rbp!r} not present in the {name} table; check --rbp/--scope.")

    mods = regulon[(regulon["rbp"] == rbp) & (regulon["q"] <= q)].copy()
    if mods.empty:
        raise SystemExit(f"no module reaches q <= {q} for {rbp} (scope={scope}, bg={bg_mode}).")
    mods["module_id"] = mods["module_id"].astype(str)
    mod_cols = ["region", "module_id", "module_size", "n_switched", "enrichment", "q"]
    # The covariate-adjusted GLM arm: effect size and interval, not just its q-value.  Ranking
    # on the adjusted interval is what separates a module that is enriched for this RBP from
    # one that is merely long, GC-rich or UTR-heavy (see --rank-by).
    mod_cols += [c for c in ("qvalue_glm", "odds_ratio", "ci_low", "ci_high", "model_status")
                 if c in mods.columns]
    mods = mods[mod_cols].rename(columns={"q": "module_q", "enrichment": "module_enrichment",
                                          "qvalue_glm": "module_q_glm",
                                          "odds_ratio": "module_or_glm",
                                          "ci_low": "module_or_ci_low",
                                          "ci_high": "module_or_ci_high",
                                          "model_status": "module_glm_status"})

    tgt = calls[(calls["rbp"] == rbp) & (calls["switched"])].copy()
    tgt["module_id"] = tgt["module_id"].astype(str)
    tgt = tgt[["region", "gene", "module_id", "go_invisible", "pool_source"]]
    ev = tgt.merge(mods, on=["region", "module_id"], how="inner")
    if ev.empty:
        raise SystemExit(f"no gene is both switched and in an enriched {rbp} module.")

    ev = ev.merge(_module_trust(), on=["region", "module_id"], how="left")
    ev = ev.merge(_gene_roles(), on=["region", "gene"], how="left")
    ev = ev.merge(_binding_support(rbp), on=["region", "module_id"], how="left")
    for c in ("binding_supported", "preferential"):
        if c in ev.columns:
            ev[c] = ev[c].astype("boolean").fillna(False).astype(bool)
    ev["pheno_module"] = ev["pool_source"].eq("switch_genes")
    return ev.sort_values(["gene", "region"]).reset_index(drop=True)


# Columns copied verbatim from the single evidence row `_best_module` selects. They have to
# travel together: a point estimate, a lower bound and an upper bound aggregated independently
# across a gene's other modules describe no model that was ever fitted.
_BEST_MODULE_CARRY = {"module_or_glm": "best_module_or",
                      "module_or_ci_low": "best_module_or_ci_low",
                      "module_or_ci_high": "best_module_or_ci_high",
                      "module_enrichment": "best_module_enrichment",
                      "module_q": "best_module_hyper_q",
                      "module_q_glm": "best_module_q_glm"}


def _best_module(ev: pd.DataFrame, rank_by: str = _DEFAULT_RANK_BY) -> pd.DataFrame:
    """Per gene, the single strongest (region, module) it belongs to, with that module's GO
    annotation -- so the biology each panel gene samples is explicit.

    Under ``rank_by="adjusted"`` "strongest" means the largest lower bound on the adjusted
    odds ratio, so a gene is reported under the module that best survives opportunity
    adjustment rather than the one with the smallest unadjusted q.
    """
    if rank_by == "adjusted" and "module_or_ci_low" in ev.columns:
        b = ev.sort_values(["module_or_ci_low", "module_q"], ascending=[False, True],
                           na_position="last")
    else:
        b = ev.sort_values(["module_q", "module_enrichment"], ascending=[True, False])
    b = b.drop_duplicates("gene")
    cols = {"region": "best_region", "module_id": "best_module"}
    if "module_go" in b.columns:
        cols["module_go"] = "best_module_go"
    cols.update({src: dst for src, dst in _BEST_MODULE_CARRY.items() if src in b.columns})
    out = b[["gene", *cols]].rename(columns=cols)
    # A table with no GLM arm leaves these unset rather than substituting the unadjusted
    # enrichment, which would rank a raw count ratio as though it were an adjusted odds ratio.
    for dst in _BEST_MODULE_CARRY.values():
        if dst not in out.columns:
            out[dst] = np.nan
    return out


def rank_candidates(ev: pd.DataFrame, min_pheno_regions: int = 1,
                    rank_by: str = _DEFAULT_RANK_BY,
                    min_adjusted_or: float = 1.0) -> pd.DataFrame:
    """Collapse the (gene, region) evidence to one ranked row per gene.

    ``rank_by``:

    ``"adjusted"`` (default) ranks on the covariate-adjusted GLM: first on whether the gene's
    best module clears ``min_adjusted_or``, then on the lower bound of that module's adjusted
    odds ratio.  The lower bound rather than the point estimate, so the ordering penalises
    both weak and imprecisely estimated effects, and a module whose adjusted interval covers 1
    cannot outrank one whose interval excludes it.  Recurrence and eCLIP support stay as the
    next keys.

    ``"hypergeometric"`` reproduces the pre-adjustment ordering (recurrence, then eCLIP, then
    unadjusted q and enrichment) and exists so the two panels can be diffed.

    Eligibility is unchanged in both modes -- it stays on the hypergeometric q.  Requiring
    adjusted significance instead would empty the panel for regulators like KHDRBS1, where no
    module reaches adjusted q <= 0.05; the adjusted arm is strong enough to *order* the
    candidates and to demote modules it contradicts, not to gate them.
    """
    ev = ev.copy()
    ev["abs_switch_r"] = pd.to_numeric(ev["switch_r"], errors="coerce").abs()
    agg = ev.groupby("gene", as_index=False).agg(
        n_regions=("region", "nunique"),
        regions=("region", lambda s: ",".join(sorted(set(s)))),
        modules=("module_id", lambda s: ",".join(sorted(set(s)))),
        best_module_q=("module_q", "min"),
        max_enrichment=("module_enrichment", "max"),
        mean_enrichment=("module_enrichment", "mean"),
        n_regions_binding=("binding_supported", "sum"),
        n_regions_pheno_module=("pheno_module", "sum"),
        n_regions_go_invisible=("go_invisible", "sum"),
        n_regions_switch_active=("switch_active", lambda s: int(pd.Series(s).fillna(False).sum())),
        module_roles=("module_role", lambda s: ",".join(sorted({str(x) for x in s.dropna()}))),
        best_pheno_fdr=("pheno_fdr", "min"),
        max_abs_switch_r=("abs_switch_r", "max"),
    )
    for c in ("n_regions_binding", "n_regions_pheno_module", "n_regions_go_invisible"):
        agg[c] = agg[c].astype(int)
    # Merge before the flags: both are statements about the reported best module, so they must
    # read that module's own estimates rather than a gene-level maximum.
    agg = agg.merge(_best_module(ev, rank_by), on="gene", how="left")
    # Does the adjusted model support this gene's best module at all?
    agg["adjusted_or_ok"] = pd.to_numeric(
        agg["best_module_or"], errors="coerce").fillna(0.0) > min_adjusted_or
    agg["adjusted_ci_excludes_null"] = pd.to_numeric(
        agg["best_module_or_ci_low"], errors="coerce").fillna(0.0) > 1.0
    agg = agg.merge(_coloc_anchors(), on="gene", how="left")
    agg["anchor_splicing_led"] = agg["anchor_splicing_led"].astype("boolean").fillna(
        False).astype(bool)
    agg["has_anchor"] = agg["anchor_verdict"].notna()
    agg["eligible"] = agg["n_regions_pheno_module"] >= min_pheno_regions

    if rank_by == "adjusted":
        keys = ["eligible", "adjusted_or_ok", "best_module_or_ci_low", "n_regions",
                "n_regions_binding", "best_module_q", "max_abs_switch_r", "gene"]
        asc = [False, False, False, False, False, True, False, True]
    else:
        keys = ["eligible", "n_regions", "n_regions_binding", "best_module_q",
                "max_enrichment", "max_abs_switch_r", "gene"]
        asc = [False, False, False, True, False, False, True]
    agg = agg.sort_values(keys, ascending=asc, na_position="last").reset_index(drop=True)
    agg.insert(0, "rank", np.arange(1, len(agg) + 1))
    return agg


def assign_tiers(cand: pd.DataFrame, rbp: str, rbp_ens: str | None,
                 n_core: int, n_anchor: int, n_goinv: int, include_expression_led: bool,
                 max_per_module: int = 3) -> pd.DataFrame:
    """Tier the ranked universe. Each gene lands in exactly one tier; tiers are filled in
    anchor -> core -> GO-invisible order so the scarce evidence class is never crowded out.

    `max_per_module` caps how many panel genes may be charged to any one *best* (region,
    module) regulon -- the module the gene is reported under. Without it the ranking fills up
    from whichever module recurs across the most regions and the panel tests one module
    instead of the regulator's program. Capping the reported module (rather than every module
    the gene touches) keeps the cap and the `Modules sampled` table describing the same thing.
    """
    best = {g: (br, bm) for g, br, bm in
            zip(cand["gene"], cand["best_region"], cand["best_module"])}
    used: dict[tuple[str, str], int] = {}

    def admit(gene: str) -> bool:
        k = best.get(gene)
        if k is None:
            return True
        if used.get(k, 0) >= max_per_module:
            return False
        used[k] = used.get(k, 0) + 1
        return True

    taken: set[str] = set()
    rows = [{"tier": 0, "tier_label": "regulator", "gene": rbp_ens or "", "rank": np.nan,
             "selection_reason": f"{rbp} itself -- the perturbation target"}]

    def fill(pool: pd.DataFrame, n: int, tier: int, label: str, reason) -> None:
        for r in pool.itertuples():
            if len([x for x in rows if x["tier"] == tier]) >= n:
                break
            if r.gene in taken or not admit(r.gene):
                continue
            taken.add(r.gene)
            rows.append({"tier": tier, "tier_label": label, "gene": r.gene, "rank": r.rank,
                         "selection_reason": reason(r)})

    order = {"splicing-led": 0, "splicing-unresolved": 1, "expression-led": 2}
    if not include_expression_led:
        order.pop("expression-led")
    anchors = cand[cand["has_anchor"] & cand["eligible"]
                   & cand["anchor_class"].isin(order)].copy()
    anchors["anchor_priority"] = anchors["anchor_class"].map(order)
    anchors = anchors.sort_values(["anchor_priority", "rank"])
    # Disease anchors are the scarcest class and are meant to span loci, so they are exempt
    # from the module cap (they are selected on coloc evidence, not on recurrence).
    for r in anchors.head(n_anchor).itertuples():
        taken.add(r.gene)
        k = best.get(r.gene)
        if k is not None:
            used[k] = used.get(k, 0) + 1
        rows.append({"tier": 2, "tier_label": "disease anchor", "gene": r.gene, "rank": r.rank,
                     "selection_reason": f"colocalized {r.anchor_traits} ({r.anchor_kinds}); "
                                         f"{r.anchor_verdict}"})

    fill(cand[cand["eligible"]], n_core, 1, "core recurrent",
         lambda r: f"switched + enriched module in {r.n_regions} region(s); "
                   f"{r.n_regions_binding} eCLIP-supported")
    fill(cand[cand["eligible"] & (cand["n_regions_go_invisible"] > 0) & (~cand["has_anchor"])],
         n_goinv, 3, "GO-invisible program",
         lambda r: f"GO-invisible module in {r.n_regions_go_invisible} region(s), "
                   f"no coloc anchor")

    panel = pd.DataFrame(rows)
    panel = panel.merge(cand.drop(columns=["rank"]), on="gene", how="left")
    return panel.sort_values(["tier", "rank"], na_position="first").reset_index(drop=True)


def motif_switch_pairs(genes: set[str], rbp: str, scope: str, bg_mode: str) -> pd.DataFrame:
    """The actual isoform pairs whose motif presence differs -- the assay readout per target.

    Re-derives the gain/loss call from the Stage-1 motif counts exactly as `rbp_regulon` does,
    but keeps the transcript pair instead of collapsing it to a boolean.
    """
    paths = _UNIT_COUNTS["rbp"][scope]
    if any(not p.exists() for p in paths):
        return pd.DataFrame()
    counts = _load_counts(paths, "rbp", "composition" if bg_mode == "composition" else "flat",
                          None)
    pres = _presence(counts[counts["rbp"] == rbp], "rbp")

    frames = []
    for tree, region in _REGIONS:
        p = region_store(tree, region, "isograph_vae", "module_interpret",
                "structure_switch_pairs.parquet")
        if not p.exists():
            continue
        d = pd.read_parquet(p)
        d["gene"] = _bare(d["gene_id"])
        d = d[d["gene"].isin(genes)]
        if not d.empty:
            d["region"] = region
            frames.append(d)
    if not frames:
        return pd.DataFrame()

    sp = pd.concat(frames, ignore_index=True)
    sp["motif_1"] = [pres.get((t, rbp), 0) for t in sp["transcript_id_1"]]
    sp["motif_2"] = [pres.get((t, rbp), 0) for t in sp["transcript_id_2"]]
    sp = sp[(sp["motif_1"] > 0) != (sp["motif_2"] > 0)].copy()
    if sp.empty:
        return sp
    sp["motif_carrier"] = np.where(sp["motif_1"] > 0, sp["transcript_id_1"],
                                   sp["transcript_id_2"])
    out = (sp.groupby(["gene", "transcript_id_1", "transcript_id_2", "motif_1", "motif_2",
                       "motif_carrier"], as_index=False)
           .agg(regions=("region", lambda s: ",".join(sorted(set(s)))),
                n_regions=("region", "nunique")))
    return out.sort_values(["gene", "n_regions"], ascending=[True, False]).reset_index(drop=True)


# --------------------------------------------------------------------------- report / emit
def _emit(df: pd.DataFrame, out_dir: Path, stem: str) -> None:
    df.to_parquet(out_dir / f"{stem}.parquet", index=False, compression="zstd")
    df.to_csv(out_dir / f"{stem}.tsv", sep="\t", index=False)


def _adjustment_effect_section(panel: pd.DataFrame, rbp: str) -> list[str]:
    """Which sampled modules the opportunity adjustment supports, and which it contradicts.

    A module can be strongly enriched on the raw hypergeometric and *depleted* once length, GC
    and UTR-vs-CDS composition are adjusted for -- that is the failure mode the reviewer asked
    about, and it has to be visible next to the panel rather than buried in a column.
    """
    p = panel[panel["tier"] > 0].copy()
    if p.empty or "best_module_or" not in p.columns:
        return []
    # Every column here is a property of the (region, module) that names the group, so it is
    # constant within it and "first" reads it off. Reducing with max/min instead would rebuild
    # the same splice this table exists to expose.
    mods = (p.groupby(["best_region", "best_module"], as_index=False)
             .agg(n=("gene", "nunique"),
                  enr=("best_module_enrichment", "first"),
                  hyper_q=("best_module_hyper_q", "first"),
                  odds=("best_module_or", "first"),
                  lo=("best_module_or_ci_low", "first"),
                  hi=("best_module_or_ci_high", "first"),
                  genes=("gene_name", lambda s: ", ".join(sorted(s)))))
    mods["odds"] = pd.to_numeric(mods["odds"], errors="coerce")
    supported = int((pd.to_numeric(mods["lo"], errors="coerce") > 1).sum())
    contradicted = int((pd.to_numeric(mods["hi"], errors="coerce") < 1).sum())
    if contradicted:
        gloss = (f" The {contradicted} contradicted module(s) are enriched on raw counts but "
                 "*depleted* once length, GC and UTR composition are adjusted for, so their "
                 "genes are carried on the unadjusted evidence alone.")
    else:
        gloss = (" No sampled module is contradicted by the adjustment. A module can be: "
                 "enriched on raw counts yet *depleted* after adjusting for length, GC and "
                 "UTR composition, and the ranking sinks any such module below the panel cut "
                 "rather than excluding it, so an absence here is a result and not a filter.")
    out = ["", "### What the opportunity adjustment changes", "",
           f"Each sampled module's unadjusted enrichment next to its adjusted odds ratio. Of "
           f"{len(mods)} sampled modules, **{supported}** have an adjusted CI entirely above 1 "
           f"and **{contradicted}** entirely below it.{gloss}", "",
           "| region / module | n | unadj. enr. | unadj. q | adj. OR (95% CI) | verdict | genes |",
           "|---|---|---|---|---|---|---|"]
    for m in mods.sort_values("odds", ascending=False, na_position="last").itertuples():
        lo = pd.to_numeric(m.lo, errors="coerce"); hi = pd.to_numeric(m.hi, errors="coerce")
        if np.isfinite(lo) and lo > 1:
            verdict = "supported"
        elif np.isfinite(hi) and hi < 1:
            verdict = "**contradicted**"
        else:
            verdict = "not resolved"
        ci = (f"{m.odds:.2f} ({lo:.2f}–{hi:.2f})"
              if np.isfinite(m.odds) and np.isfinite(lo) and np.isfinite(hi) else "—")
        out.append(f"| {m.best_region}/{m.best_module} | {m.n} | {m.enr:.2f} | "
                   f"{m.hyper_q:.1e} | {ci} | {verdict} | {m.genes} |")
    return out


def _or_cell(r) -> str:
    """Adjusted OR with its CI; bold when the interval excludes 1, em-dash when unestimated."""
    o = pd.to_numeric(getattr(r, "best_module_or", np.nan), errors="coerce")
    lo = pd.to_numeric(getattr(r, "best_module_or_ci_low", np.nan), errors="coerce")
    hi = pd.to_numeric(getattr(r, "best_module_or_ci_high", np.nan), errors="coerce")
    if not np.isfinite(o):
        return "—"
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return f"{o:.2f}"
    txt = f"{o:.2f} ({lo:.2f}–{hi:.2f})"
    return f"**{txt}**" if lo > 1.0 else txt


def _write_report(panel: pd.DataFrame, cand: pd.DataFrame, pairs: pd.DataFrame, rbp: str,
                  q: float, scope: str, bg_mode: str, max_per_module: int,
                  out_dir: Path, rank_by: str = _DEFAULT_RANK_BY,
                  min_adjusted_or: float = 1.0) -> None:
    elig = cand[cand["eligible"]]
    qg = pd.to_numeric(elig["best_module_q_glm"], errors="coerce")
    n_glm = int((qg <= 0.05).sum())
    best_qg = float(qg.min()) if qg.notna().any() else float("nan")
    n_or = int(elig["adjusted_or_ok"].sum()) if "adjusted_or_ok" in elig.columns else 0
    n_ci = (int(elig["adjusted_ci_excludes_null"].sum())
            if "adjusted_ci_excludes_null" in elig.columns else 0)
    rank_rule = (
        "ranked by the covariate-adjusted GLM: genes whose best module has an adjusted odds "
        f"ratio above {min_adjusted_or:g} first, then by the lower bound of that module's "
        "adjusted 95% CI, then by recurrence across brain regions, then by the number of "
        "regions where the module's regulon call is eCLIP-binding-supported, then by "
        "unadjusted module q, then by the gene's own switch-driver strength (max |switch_r|), "
        "then alphabetically"
        if rank_by == "adjusted" else
        "ranked by recurrence across brain regions, then by the number of regions where the "
        "module's regulon call is eCLIP-binding-supported, then by module enrichment "
        "strength, then by the gene's own switch-driver strength (max |switch_r|), then "
        "alphabetically")
    L = [f"# {rbp} perturbation target panel", "",
         f"Generated by `isograph_benchmark.real_data.rbp_target_panel` "
         f"(scope={scope}, background={bg_mode}, module q <= {q}, "
         f"max {max_per_module} genes per module, rank_by={rank_by}).", "",
         "## How the list was derived", "",
         f"A gene is **eligible** only if both layers agree: its IsoGraph switch pair changes "
         f"{rbp} motif presence (`rbp_switch_calls.switched`), *and* its module is an "
         f"independently enriched {rbp} regulon (`rbp_regulon.q <= {q}`) that is itself "
         f"phenotype-associated. Eligible genes are {rank_rule}.", "",
         f"Eligibility stays on the **hypergeometric** q because no {rbp} module reaches "
         f"adjusted q <= 0.05 (best adjusted q = {best_qg:.3f}); gating on the adjusted arm "
         "would leave no panel at all. The adjusted arm is used to *order* the candidates and "
         "to demote modules it contradicts.", "",
         f"Candidate universe: **{len(elig)} eligible genes** "
         f"({len(cand) - len(elig)} more are switched-and-enriched but never in a "
         "phenotype-associated module); "
         f"{int((elig['n_regions'] > 1).sum())} recur in >1 region; "
         f"{int((elig['n_regions_binding'] > 0).sum())} have eCLIP support in >=1 region.", "",
         "## Panel", "",
         "`adj. OR` is the covariate-adjusted odds ratio of the gene's best module with its "
         "95% CI; **bold** marks an interval that excludes 1. `module q` remains the "
         "unadjusted hypergeometric value the eligibility filter uses.", "",
         "| tier | gene | Ensembl | regions | eCLIP | module q | enr. | adj. OR (95% CI) | "
         "GO-inv | best module | module biology | why |",
         "|------|------|---------|---------|-------|----------|------|------------------|"
         "--------|-------------|----------------|-----|"]
    for r in panel.itertuples():
        if r.tier == 0:
            L.append(f"| 0 | **{r.gene_name}** | {r.gene or '—'} | — | — | — | — | — | — | — "
                     f"| — | {r.selection_reason} |")
            continue
        go = getattr(r, "best_module_go", "") or ("GO-invisible" if r.n_regions_go_invisible
                                                  else "—")
        or_txt = _or_cell(r)
        L.append(
            f"| {r.tier} | {r.gene_name} | {r.gene} | {int(r.n_regions)} | "
            f"{int(r.n_regions_binding)} | {r.best_module_q:.1e} | {r.max_enrichment:.2f} | "
            f"{or_txt} | "
            f"{int(r.n_regions_go_invisible)} | {r.best_region}/{r.best_module} | {go} | "
            f"{r.selection_reason} |")

    mods = (panel[panel["tier"] > 0].groupby(["best_region", "best_module"], as_index=False)
            .agg(n=("gene", "nunique"),
                 genes=("gene_name", lambda s: ", ".join(sorted(s)))))
    L += ["", "### Modules sampled", "",
          f"The panel draws from **{len(mods)} distinct (region, module) regulons** "
          f"(cap: {max_per_module} genes each).", "",
          "| region / module | n | genes |", "|---|---|---|"]
    for m in mods.sort_values("n", ascending=False).itertuples():
        L.append(f"| {m.best_region}/{m.best_module} | {m.n} | {m.genes} |")

    if rank_by == "adjusted" and "best_module_or" in cand.columns:
        L += _adjustment_effect_section(panel, rbp)

    L += ["", "## Caveats", "",
          f"- **No {rbp} module survives opportunity adjustment at FDR.** The length/GC/UTR-"
          f"adjusted GLM (`qvalue_glm`) reaches q <= 0.05 for {n_glm}/{len(elig)} of these "
          f"genes' best modules (best adjusted q = {best_qg:.3f}); {n_or} best modules have "
          f"an adjusted odds ratio above {min_adjusted_or:g} and {n_ci} have an adjusted CI "
          "that excludes 1 before multiplicity. The module-level claim is therefore 'enriched "
          "relative to the switch-gene pool', not 'enriched after full opportunity "
          "adjustment'. This is a prioritized hypothesis set for a perturbation experiment, "
          "and must not be described as an FDR-significant regulon.",
          "- Motif gain/loss is sequence-predicted. eCLIP support (`eCLIP` column) is "
          "module-level and isoform-preferential, not per-target binding evidence.",
          "- Genes whose region contributed via the fallback full module-gene pool "
          "(`pheno_module == False`) come from a module with no FDR-significant phenotype "
          "association in that region; `n_regions_pheno_module` records how many regions "
          "carry the stronger pool.",
          f"- Recurrence is module-driven: the most widely recurrent {rbp} regulon supplies "
          "far more high-ranked genes than the module cap admits. Check the *Modules sampled* "
          "table before assuming the panel spans distinct biology, and check that the sampled "
          "modules' cell-type biology matches the model system the panel will be tested in "
          "(an immune/inflammatory module is not testable in a purely neuronal culture).",
          "- Freeze this ranking **before** any collaborator expression (TPM) screen. Use the "
          "screen only to decide which prespecified targets are testable in their cell model; "
          "selecting targets because they happen to be well expressed is a post-hoc step.", ""]
    if not pairs.empty:
        L += ["## Motif-differential switch pairs (assay readout)", "",
              f"{len(pairs)} isoform pairs across {pairs['gene'].nunique()} panel genes gain or "
              f"lose the {rbp} motif; the `motif_carrier` isoform is the one predicted to "
              "respond to the perturbation. Full table: `rbp_target_switch_pairs.tsv`.", ""]
    (out_dir / "RBP_TARGET_PANEL.md").write_text("\n".join(L) + "\n")


def run(rbp: str = "KHDRBS1", q: float = 0.05, scope: str = "mature",
        bg_mode: str = "composition", n_core: int = 14, n_anchor: int = 5,
        n_goinv: int = 5, include_expression_led: bool = False,
        max_per_module: int = 3, min_pheno_regions: int = 1,
        gtf_path: Path | None = None, rank_by: str = _DEFAULT_RANK_BY,
        min_adjusted_or: float = 1.0) -> pd.DataFrame:
    if rank_by not in _RANK_BY_CHOICES:
        raise SystemExit(f"--rank-by must be one of {_RANK_BY_CHOICES} (got {rank_by!r}).")
    ev = build_evidence(rbp, q, scope, bg_mode)
    cand = rank_candidates(ev, min_pheno_regions, rank_by, min_adjusted_or)

    sym = gene_symbol_map(gtf_path or DEFAULT_GTF_PATH)
    ens_of = {v: k for k, v in sym.items()}
    for d in (ev, cand):
        d["gene_name"] = d["gene"].map(sym).fillna(d["gene"])

    panel = assign_tiers(cand, rbp, ens_of.get(rbp), n_core, n_anchor, n_goinv,
                         include_expression_led, max_per_module)
    panel["gene_name"] = panel["gene"].map(sym).fillna(panel["gene"])
    panel.loc[panel["tier"] == 0, "gene_name"] = rbp
    # the tier-0 row has no evidence, so the count columns upcast to float; keep them integral
    for c in ("n_regions", "n_regions_binding", "n_regions_pheno_module",
              "n_regions_go_invisible", "n_regions_switch_active", "rank"):
        if c in panel.columns:
            panel[c] = panel[c].astype("Int64")
    pairs = motif_switch_pairs(set(panel["gene"]) - {""}, rbp, scope, bg_mode)
    if not pairs.empty:
        pairs["gene_name"] = pairs["gene"].map(sym).fillna(pairs["gene"])

    out_dir = ensure_dir(_OUT_ROOT / rbp)
    _emit(ev, out_dir, "rbp_target_evidence")
    _emit(cand, out_dir, "rbp_target_candidates")
    _emit(panel, out_dir, "rbp_target_panel")
    if not pairs.empty:
        _emit(pairs, out_dir, "rbp_target_switch_pairs")
    panel[["tier", "tier_label", "gene_name", "gene", "n_regions", "n_regions_binding",
           "selection_reason"]].rename(columns={"gene": "ensembl_gene_id"}).to_csv(
        out_dir / "expression_check_list.tsv", sep="\t", index=False)
    _write_report(panel, cand, pairs, rbp, q, scope, bg_mode, max_per_module, out_dir,
                  rank_by, min_adjusted_or)

    elig = cand[cand["eligible"]]
    print(f"[{rbp}] eligible genes: {len(elig)}/{len(cand)} "
          f"({int((elig['n_regions'] > 1).sum())} recurrent, "
          f"{int((elig['n_regions_binding'] > 0).sum())} eCLIP-supported)")
    n_or = int(elig["adjusted_or_ok"].sum()) if "adjusted_or_ok" in elig.columns else 0
    n_ci = (int(elig["adjusted_ci_excludes_null"].sum())
            if "adjusted_ci_excludes_null" in elig.columns else 0)
    print(f"[{rbp}] adjusted GLM support among eligible: {n_or} best modules with OR > "
          f"{min_adjusted_or}, {n_ci} whose adjusted CI excludes 1 (rank_by={rank_by})")
    print(f"[{rbp}] panel: {len(panel)} genes "
          f"({', '.join(f'tier {t}: {n}' for t, n in panel['tier'].value_counts().sort_index().items())})")
    print(f"[{rbp}] -> {out_dir}")
    return panel


def main() -> None:
    p = argparse.ArgumentParser(
        description="Nominate a perturbation target panel for a candidate splicing regulator.")
    p.add_argument("--rbp", default="KHDRBS1", help="regulator to build the panel for")
    p.add_argument("--q", type=float, default=0.05,
                   help="module-level regulon FDR threshold (default 0.05)")
    p.add_argument("--scope", choices=("mature", "intronic", "combined"), default="mature",
                   help="which Stage-2 tables to consume (default: mature, canonical)")
    p.add_argument("--bg-mode", choices=("composition", "flat"), default="composition",
                   help="composition-matched background (default) or the legacy flat 0.25")
    p.add_argument("--n-core", type=int, default=14, help="tier-1 size (default 14)")
    p.add_argument("--n-anchor", type=int, default=5, help="tier-2 size (default 5)")
    p.add_argument("--n-goinv", type=int, default=5, help="tier-3 size (default 5)")
    p.add_argument("--include-expression-led", action="store_true",
                   help="allow eQTL-led coloc genes into tier 2 (flagged; off by default "
                        "because expression-led loci confound a splicing readout)")
    p.add_argument("--max-per-module", type=int, default=3,
                   help="cap on panel genes drawn from one (region, module) regulon "
                        "(default 3; 0 disables the cap)")
    p.add_argument("--min-pheno-regions", type=int, default=1,
                   help="require the gene's module to be phenotype-associated in at least "
                        "this many regions (default 1)")
    p.add_argument("--gtf", default=None, help="GENCODE GTF for gene symbols (cached)")
    p.add_argument("--rank-by", choices=_RANK_BY_CHOICES, default=_DEFAULT_RANK_BY,
                   help="rank on the covariate-adjusted GLM odds ratio (default) or reproduce "
                        "the pre-adjustment hypergeometric ordering")
    p.add_argument("--min-adjusted-or", type=float, default=1.0,
                   help="under --rank-by adjusted, genes whose best module has an adjusted "
                        "odds ratio at or below this sort below those that clear it "
                        "(default 1.0, i.e. no enrichment after adjustment)")
    args = p.parse_args()
    run(args.rbp, args.q, args.scope, args.bg_mode, args.n_core, args.n_anchor, args.n_goinv,
        args.include_expression_led,
        args.max_per_module if args.max_per_module > 0 else 10 ** 9,
        args.min_pheno_regions, Path(args.gtf) if args.gtf else None,
        args.rank_by, args.min_adjusted_or)


if __name__ == "__main__":
    main()
