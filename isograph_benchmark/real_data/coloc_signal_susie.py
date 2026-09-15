"""Signal-level colocalization (`coloc.susie`) of switch genes against brain xQTLs.

WHY THIS EXISTS
---------------
The genetic-anchoring layer currently carries two colocalization estimators, and both
answer a coarser question than the paper asks.

`coloc.abf` (``coloc_modality_contrast``) assumes **at most one causal variant per trait
in the analysed window**. At loci that plainly carry more than one signal -- SNCA,
PICALM, and anything in a gene-dense or long-range-LD region -- that assumption is not a
technicality: two independent QTL signals, only one of which is shared with the GWAS,
push posterior mass into H3 and the locus reads as "distinct causal variants" when it is
really "one shared signal plus one private one". eCAVIAR CLPP
(``05_genetic_anchoring/_h/03b.coloc_clpp.R``) is signal-aware on the QTL side but
consumes GTEx's shipped credible sets, so its yield is small and its posteriors are
individually modest.

`coloc.susie` colocalizes **signal against signal**: it fine-maps both traits, then asks
for every pair of credible sets whether that pair shares a causal variant. That is the
estimator the claim actually needs, and it is why the hierarchy this module implements is

    signal-level SuSiE coloc  >  coloc.abf  >  CLPP

CLPP is **retained as orthogonal sensitivity evidence and is not relabelled a wrong
statistic**: eCAVIAR and coloc make different assumptions (eCAVIAR does not model a
shared-variant prior; coloc does), so agreement between them is informative and
disagreement is worth reporting rather than hiding.

WHY SUSIE HAS TO BE RE-FIT ON THE QTL SIDE
------------------------------------------
GTEx v11 ships ``*.SuSiE_summary.parquet`` -- ``phenotype_id, gene_id, gene_name,
biotype, variant_id, pip, af, cs_id, cs_size``. That is a PIP per credible-set variant.
`coloc.susie` needs the per-effect log Bayes factors (`lbf_variable`), which are only in
a full SuSiE object, and GTEx does not release one. So the QTL side is re-fit here with
`susie_rss` on the cis nominal statistics from the v11 all-pairs release.

THE LD PROBLEM, STATED PLAINLY
------------------------------
Re-fitting needs an LD matrix, and the only one available is the 1000G EUR Phase 3 panel
already built per locus for the CLPP layer. GTEx brain donors are **not** a pure EUR
sample, so the QTL side is fine-mapped under a reference LD that does not match the
sample LD that generated the z-scores. `susie_rss` is known to manufacture credible sets
under that mismatch. Three guards, all of which are reported rather than assumed:

  1. Both sides use the SAME LD matrix and the SAME SNP set, so the mismatch is at least
     common to the two traits rather than differential between them.
  2. **The GTEx-agreement filter.** A re-fit QTL credible set is flagged
     ``cs_matches_gtex`` only when it shares a variant with GTEx's OWN shipped credible
     set for that phenotype. GTEx fit theirs with in-sample genotypes, so agreement is an
     external check that the re-fit signal is real and not an LD artefact. Results are
     reported with and without the filter; the filtered set is primary.
  3. ``susieR::kriging_rss`` diagnostics are recorded per fit, so a locus whose z-scores
     are inconsistent with the reference LD can be identified instead of trusted.

Where SuSiE fails, or finds no credible set on either side, the cell **falls back to
`coloc.abf`** and is flagged ``estimator = "abf"``. The grid therefore stays complete and
a reader can see exactly which cells the signal-level estimator could and could not
speak to.

PRIOR SENSITIVITY
-----------------
`PP4 >= 0.8` is a convention, not a fact, and coloc's posteriors depend on the
association priors (`p1`, `p2`) and the colocalization prior (`p12`). For every cell
reaching PP4 >= 0.5 in either modality the R stage sweeps `p12` and records the range
over which the H4 call survives, so a hit can be reported as "colocalized for p12 in
[a, b]" rather than as a bare posterior.

STAGES
------
  --stage prep   Freeze the target grid. Deliberately reads the grid `coloc.abf` already
                 ran on, so the two estimators are compared on identical cells rather
                 than on two independently rebuilt gene lists. Also writes the GTEx
                 credible sets that back the agreement filter. Login-node safe.
  [R stage A]    05_genetic_anchoring/_h/03c.coloc_gwas_susie.R <analysis>
                 -- fits and caches the per-locus GWAS SuSiE, once, for every tissue task
                 downstream to reuse.
  [R stage B]    05_genetic_anchoring/_h/06b.coloc_signal_susie.R <analysis> <tissue>
                 -- QTL SuSiE + coloc.susie (+ abf fallback + p12 sweep) per cell.
  --stage meta   Assemble, apply the estimator hierarchy, run the paired modality tests,
                 write the report.

Outputs under 05_genetic_anchoring/_m/coloc_signal_susie/.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data import coloc_modality_contrast as cmc

# The gene pools are defined once, in the abf module; this module never re-derives them.
ARMS_GENE_POOL = cmc.ARMS_GENE_POOL
ANALYSES = cmc.ANALYSES
GTEX_BRAIN = cmc.GTEX_BRAIN

# Pre-specified, and identical to the abf layer so the estimators are comparable.
MIN_SHARED_SNPS = cmc.MIN_SHARED_SNPS
PP4_CALL = cmc.PP4_CALL
PP4_CALL_LOOSE = cmc.PP4_CALL_LOOSE
P12_PRIMARY = cmc.P12_PRIMARY
# Swept by the R stage for every cell that reaches PP4 >= 0.5 in either modality.
P12_SWEEP: tuple[float, ...] = (1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
# No seed constant: every step of this layer is deterministic given its inputs.
# `susie_rss` has no RNG, and coloc's posteriors are closed-form Bayes factors. A seed
# here would imply a stochastic step that does not exist.


# --------------------------------------------------------------------------- #
# Evidence-strength descriptors, carried on every hierarchy cell
# --------------------------------------------------------------------------- #
# Neither descriptor gates a nomination. They exist so a PP4 is never quoted without the
# two facts that most change what it is worth: how far down the prior the call survives,
# and -- for a coloc.abf cell -- why signal-level coloc did not score it.

# PP4 is monotone in p12 (p12 enters only the H4 prior), so the smallest p12 at which the
# call holds defines the range it survives. Both estimators reach 1e-6; the abf layer stops
# at the primary 1e-5, which cannot change the label of a cell that calls at the primary.
PRIOR_ROBUSTNESS = ("robust", "intermediate", "primary_prior", "permissive_prior", "none")

FALLBACK_REASONS = (
    "gwas_locus_over_max_snps",   # stage A never fit the locus: the compute guard
    "gwas_too_few_snps",          # too few usable SNPs after harmonization
    "gwas_no_credible_set",       # fit, but the GWAS signal did not fine-map
    "gwas_not_in_stage_a",        # no stage-A status row at all
    "no_qtl_credible_set",        # GWAS fine-mapped; no signal pair was emitted
    "qtl_cs_not_matching_gtex",   # a QTL credible set existed; GTEx's own disagreed
    "susie_pair_filtered",        # a matched pair existed and meta removed it (SNP floor)
)


def p12_survival(long: pd.DataFrame, keys: list[str],
                 call: float = PP4_CALL) -> pd.DataFrame:
    """Per cell: PP4 at the most conservative prior tested, and the smallest p12 that calls.

    `long` holds one row per (cell, p12) or per (cell, signal pair, p12). A cell's PP4 at
    each prior is its best pair there -- the same collapse as at the primary prior.
    """
    cols = [*keys, "p12_sweep_min", "PP4_at_p12_sweep_min", "p12_min_call"]
    if long.empty:
        return pd.DataFrame(columns=cols)
    per = long.groupby([*keys, "p12"], dropna=False, as_index=False)["PP4"].max()
    # idxmin, not sort + first(): first() skips NaN column-wise and would pair the lowest
    # p12 with a PP4 read at a different prior.
    lo = (per.loc[per.groupby(keys, dropna=False)["p12"].idxmin()]
             .rename(columns={"p12": "p12_sweep_min", "PP4": "PP4_at_p12_sweep_min"}))
    calls = (per[per["PP4"] >= call]
             .groupby(keys, dropna=False, as_index=False)["p12"].min()
             .rename(columns={"p12": "p12_min_call"}))
    return lo.merge(calls, on=keys, how="left")[cols]


def prior_robustness(p12_min_call) -> str:
    """Label the prior range a call survives. `robust` means it holds at p12 = 1e-6."""
    if p12_min_call is None or not np.isfinite(p12_min_call):
        return "none"
    lo = min(P12_SWEEP)
    if p12_min_call <= lo or np.isclose(p12_min_call, lo):
        return "robust"
    if np.isclose(p12_min_call, P12_PRIMARY):
        return "primary_prior"
    return "intermediate" if p12_min_call < P12_PRIMARY else "permissive_prior"


def load_gwas_status(root: Path) -> pd.DataFrame:
    """Every stage-A `_status.tsv` under one results root, tagged with its analysis."""
    cols = ["analysis", "LOCUS_ID", "fitted", "n_cs", "reason"]
    gw = root / "gwas_susie"
    parts = []
    for d in (sorted(p for p in gw.iterdir() if p.is_dir()) if gw.exists() else []):
        f = d / "_status.tsv"
        if f.exists():
            s = pd.read_csv(f, sep="\t")
            s["analysis"] = d.name
            parts.append(s)
    return pd.concat(parts, ignore_index=True)[cols] if parts else pd.DataFrame(columns=cols)


def fallback_reasons(cells: pd.DataFrame, status: pd.DataFrame, pairs: pd.DataFrame,
                     require_gtex_match: bool = True) -> pd.Series:
    """Why each coloc.abf cell was not scored by coloc.susie; None on susie cells.

    Walks the pipeline in order and reports the FIRST place a cell left it, because the
    reasons are different claims about the evidence. A locus the SNP guard kept from SuSiE
    has had no signal-level test at all. A GWAS that fit but fine-mapped to no credible set
    WAS tested and has no resolvable signal, which weakens any colocalization claimed at it.
    A QTL credible set that GTEx's in-sample fine-mapping does not share is the
    reference-LD artefact the agreement filter exists to catch.

    `no_qtl_credible_set` is inferred from absence: stage B emits a row only when both
    sides fine-map, so it also absorbs the rare QTL fit that errored or shared fewer than
    50 SNPs with the GWAS.
    """
    keys = ["analysis", "trait", "LOCUS_ID", "gene", "tissue", "modality"]
    out = pd.Series([None] * len(cells), index=cells.index, dtype=object)
    fb = cells["estimator"] == "abf"
    if not fb.any():
        return out

    st = status.drop_duplicates(["analysis", "LOCUS_ID"])
    st = {(str(a), str(l)): (str(r), n) for a, l, r, n in
          st[["analysis", "LOCUS_ID", "reason", "n_cs"]].itertuples(index=False)}

    def _keys(d: pd.DataFrame) -> set:
        return set(map(tuple, d[keys].astype(str).to_numpy())) if len(d) else set()

    emitted = _keys(pairs)
    matched = (_keys(pairs[pairs["cs_matches_gtex"].fillna(False).astype(bool)])
               if require_gtex_match and "cs_matches_gtex" in pairs.columns else emitted)

    for i, k in zip(cells.index[fb], map(tuple, cells.loc[fb, keys].astype(str).to_numpy())):
        s = st.get((k[0], k[2]))
        if s is None:
            why = "gwas_not_in_stage_a"
        elif "MAX_SNPS" in s[0]:
            why = "gwas_locus_over_max_snps"
        elif "usable SNPs" in s[0]:
            why = "gwas_too_few_snps"
        elif "credible set" in s[0] or (pd.notna(s[1]) and s[1] == 0):
            why = "gwas_no_credible_set"
        elif k not in emitted:
            why = "no_qtl_credible_set"
        elif k not in matched:
            why = "qtl_cs_not_matching_gtex"
        else:
            why = "susie_pair_filtered"
        out.at[i] = why
    return out


# The GWAS SuSiE SNP guard of the uniform primary grid (03c.coloc_gwas_susie.R). A run at
# any other value is a scoped SENSITIVITY arm with its own root, so it can never replace
# primary shards -- which is what the 2026-09-10 aging__ad recovery at 30,000 did before
# this split existed, overwriting the 12,000 AD cache and shards in place.
MAX_SNPS_PRIMARY = 12_000


def signal_root(arm: str = "switch", max_snps: int = MAX_SNPS_PRIMARY) -> Path:
    """Results root for one (gene-pool arm, GWAS SNP guard). Resolves; creates nothing.

    Mirrors what the R stages write: `<signal>/[sensitivity/max_snps_<N>/][arms/<arm>/]`.
    """
    base = stage_out("anchoring.coloc_signal")
    if int(max_snps) != MAX_SNPS_PRIMARY:
        base = base / "sensitivity" / f"max_snps_{int(max_snps)}"
    return base if arm == "switch" else base / "arms" / arm


def out_dir(arm: str = "switch") -> Path:
    """Output dir for a gene-pool arm, mirroring `coloc_modality_contrast.out_dir`."""
    return ensure_dir(signal_root(arm))


# --------------------------------------------------------------------------- #
# Stage: prep
# --------------------------------------------------------------------------- #
def gtex_credible_sets(genes: set[str]) -> pd.DataFrame:
    """GTEx's own shipped SuSiE credible sets, for the agreement filter.

    This is the reference the re-fit QTL credible sets are checked against. It is the
    same loader the CLPP layer's QTL side uses, so "agrees with GTEx" means the same
    thing in both places.
    """
    from isograph_benchmark.real_data.coloc_prep import (
        DEFAULT_XQTL_DIR, _GTEX_BRAIN, load_qtl_credible_sets,
    )
    cs = load_qtl_credible_sets(DEFAULT_XQTL_DIR, list(_GTEX_BRAIN), genes=genes)
    if cs.empty:
        raise SystemExit("no GTEx brain credible sets for the target genes")
    return cs[["gene", "phenotype_id", "variant_id", "pip", "cs_id", "cs_size",
               "tissue", "kind"]]


def run_prep(arm: str = "switch", dest: Path | None = None) -> Path:
    """Freeze the target grid, taken from the arm the abf layer already ran.

    Rebuilding the gene list independently would risk the two estimators diverging for a
    reason that has nothing to do with the estimator. So prep *reads* the abf arm's
    frozen targets and fails loudly if they are absent.
    """
    if arm not in ARMS_GENE_POOL:
        raise SystemExit(f"unknown arm {arm!r}; choose from {ARMS_GENE_POOL}")
    src = cmc.out_dir(arm)
    need = ["targets.parquet", "sqtl_representative.parquet", "gwas_meta.tsv"]
    missing = [n for n in need if not (src / n).exists()]
    if missing:
        raise SystemExit(
            f"missing {missing} under {src}; run "
            f"`coloc_modality_contrast --stage prep --arm {arm}` first so both "
            f"estimators share one target grid")
    dest = ensure_dir(dest or out_dir(arm))

    targets = pd.read_parquet(src / "targets.parquet")
    rep = pd.read_parquet(src / "sqtl_representative.parquet")
    gm = pd.read_csv(src / "gwas_meta.tsv", sep="\t")
    targets.to_parquet(dest / "targets.parquet", index=False)
    rep.to_parquet(dest / "sqtl_representative.parquet", index=False)
    gm.to_csv(dest / "gwas_meta.tsv", sep="\t", index=False)

    genes = set(targets["gene"])
    cs = gtex_credible_sets(genes)
    cs.to_parquet(dest / "gtex_credible_sets.parquet", index=False)

    # Work list for the SLURM array: one task per (analysis, tissue).
    work = pd.DataFrame(
        [(a, t) for a in sorted(targets["analysis"].unique()) for t in GTEX_BRAIN],
        columns=["analysis", "tissue"])
    work.to_csv(dest / "work_list.tsv", sep="\t", index=False)

    print(f"  target grid mirrored from {src}")
    print(f"  {len(targets):,} (analysis, locus, gene) targets, "
          f"{targets['gene'].nunique():,} genes, "
          f"{targets['LOCUS_ID'].nunique():,} loci")
    print(f"  GTEx credible sets for the agreement filter: {len(cs):,} rows, "
          f"{cs['gene'].nunique():,} genes, {cs['phenotype_id'].nunique():,} phenotypes")
    print(f"  work list: {len(work)} (analysis, tissue) array tasks -> {dest/'work_list.tsv'}")
    print(f"  p12 sweep: {P12_SWEEP}")
    return dest


# --------------------------------------------------------------------------- #
# Stage: meta
# --------------------------------------------------------------------------- #
def results_dir(src: Path, sqtl: str = "representative") -> Path:
    """Where one (arm, sqtl) run's shards and meta outputs live.

    The gene-pool arm already splits on directory; the sqtl phenotype choice needs the
    same, because both modes name their shards <analysis>__<tissue>.parquet. Without the
    split an all-introns run overwrites the representative-intron results it exists to be
    compared against. Prep artefacts (targets, representative map, GTEx credible sets,
    GWAS SuSiE cache) are shared and stay in the arm dir.
    """
    return src if sqtl == "representative" else src / "all_introns"


def _load_susie(src: Path) -> pd.DataFrame:
    d = src / "susie"
    if not d.exists():
        raise SystemExit(f"no susie/ results in {src}; run 06b.coloc_signal_susie.sh first")
    parts = [pd.read_parquet(p) for p in sorted(d.iterdir()) if p.suffix == ".parquet"]
    if not parts:
        raise SystemExit(f"no parquet files under {d}")
    return pd.concat(parts, ignore_index=True)


def best_signal_pair(pairs: pd.DataFrame, require_gtex_match: bool = True
                     ) -> pd.DataFrame:
    """Collapse signal-pair rows to one row per (cell, modality).

    `coloc.susie` returns a posterior for every (GWAS credible set x QTL credible set)
    pair. The cell-level statement the paper makes is "does this gene's splicing signal
    colocalize with the disease signal at all", so the cell keeps its **best** pair --
    but the number of pairs tested is carried alongside as `n_signal_pairs`, because a
    maximum over many pairs is a different object from a single test and should not be
    quoted as if it were one.

    `require_gtex_match` restricts to QTL signals that agree with GTEx's own credible
    set. That is the primary arm: it is what keeps a reference-LD fine-mapping artefact
    from being promoted to a locus nomination.
    """
    d = pairs
    if require_gtex_match and "cs_matches_gtex" in d.columns:
        d = d[d["cs_matches_gtex"].fillna(False).astype(bool)]
    # The shared-SNP floor the abf layer applies. Stage B emits pairs from 50 shared SNPs
    # and has always deferred the real floor to this stage, which until 2026-09-11 never
    # applied it. No signal-level cell sat below it (minimum 1,525 shared SNPs), so adding
    # it changed no result; it closes the gap rather than correcting one.
    if "n_shared" in d.columns:
        d = d[d["n_shared"] >= MIN_SHARED_SNPS]
    if d.empty:
        return d
    keys = ["analysis", "trait", "LOCUS_ID", "gene", "tissue", "modality"]

    # The R stage emits one row per (signal pair x p12) so the prior can be swept for
    # free. The PRIMARY result must be read at the pre-specified prior only: taking the
    # maximum over the sweep would silently report every locus at the most permissive
    # prior in the grid, which is a different and much weaker claim. (SNCA/LBD reads
    # 0.974 at p12=1e-5 and 0.997 at 1e-4 -- the sweep is the sensitivity, not the
    # answer.)
    prim = d[d["p12"] == P12_PRIMARY]
    if prim.empty:
        raise SystemExit(f"no rows at the primary prior p12={P12_PRIMARY}")
    best = (prim.sort_values("PP4", ascending=False)
                .groupby(keys, dropna=False, as_index=False).first())

    # Prior sensitivity, carried alongside: the range the cell's best PP4 spans across the
    # whole sweep, so a hit can be reported as surviving a RANGE of priors.
    sw = (d.groupby(keys + ["p12"], dropna=False, as_index=False)["PP4"].max()
            .groupby(keys, dropna=False, as_index=False)
            .agg(PP4_min_over_p12=("PP4", "min"), PP4_max_over_p12=("PP4", "max")))
    # Count DISTINCT signal pairs, not rows: the R stage emits each pair once per p12, so
    # a raw row count multiplies by the sweep length. Counting (idx1, idx2) is exact even
    # when coloc.susie fails for one prior and not another, which dividing by the sweep
    # length would silently get wrong.
    #
    # phenotype_id belongs in that identity. SuSiE numbers credible sets from 1 WITHIN
    # each fit, so in the all-introns arm the first credible set of intron A and the first
    # of intron B are both (idx1=1, idx2=1) and deduplicate onto each other -- collapsing
    # a gene's whole intron set to the pair count of a single intron. In the
    # representative arm phenotype_id is constant per cell, so adding it changes nothing.
    idx = [c for c in ("phenotype_id", "idx1", "idx2") if c in pairs.columns]
    if "idx1" in idx and "idx2" in idx:
        n = (pairs.drop_duplicates(keys + idx)
                  .groupby(keys, dropna=False, as_index=False)
                  .size().rename(columns={"size": "n_signal_pairs"}))
    else:
        n = (pairs.groupby(keys, dropna=False, as_index=False)
                  .size().rename(columns={"size": "n_signal_pairs"}))
    out = (best.merge(sw, on=keys, how="left").merge(n, on=keys, how="left")
               .merge(p12_survival(d, keys), on=keys, how="left"))

    # How many intron phenotypes the winning PP4 was a maximum over. The cell keeps its
    # best intron (phenotype_id above names it), but in the all-introns arm that best is
    # selected over ~12 introns per gene, and a maximum over 12 is not the same object as
    # a single test -- exactly the reason n_signal_pairs is carried. 1 everywhere in the
    # representative arm, which is what makes the two arms comparable at all.
    if "phenotype_id" in pairs.columns:
        nph = (pairs.drop_duplicates(keys + ["phenotype_id"])
                    .groupby(keys, dropna=False, as_index=False)
                    .size().rename(columns={"size": "n_phenotypes_tested"}))
        out = out.merge(nph, on=keys, how="left")
    return out


def apply_hierarchy(susie: pd.DataFrame, abf: pd.DataFrame) -> pd.DataFrame:
    """Signal-level result where SuSiE could speak; `coloc.abf` where it could not.

    A cell is scored by SuSiE only when BOTH traits yielded at least one credible set
    (and, in the primary arm, the QTL set agrees with GTEx). Otherwise the abf posterior
    stands in, flagged, so the grid stays complete and the fallback rate is visible
    rather than silently shrinking the denominator.
    """
    keys = ["analysis", "trait", "LOCUS_ID", "gene", "tissue", "modality"]
    s = susie.copy()
    s["estimator"] = "susie"
    a = abf.copy()
    a["estimator"] = "abf"
    have = set(map(tuple, s[keys].astype(str).to_numpy())) if len(s) else set()
    mask = [tuple(r) not in have for r in a[keys].astype(str).to_numpy()]
    fallback = a[pd.Series(mask, index=a.index)]
    # Concat only the non-empty frames: pandas warns (and will change dtype behaviour)
    # when an all-NA or empty frame takes part, and an empty SuSiE side is a normal
    # outcome here -- a GWAS that fine-maps nowhere leaves every cell on abf.
    parts = [d for d in (s, fallback) if len(d)]
    if not parts:
        return s.iloc[0:0]
    return pd.concat(parts, ignore_index=True)


def _abf_cells(arm: str) -> pd.DataFrame:
    """Primary-arm `coloc.abf` posteriors, in the long (cell, modality) shape."""
    src = cmc.out_dir(arm)
    abf = cmc._load_abf(src)
    d = abf[(abf["p12"] == P12_PRIMARY) & (abf["nsnps"] >= MIN_SHARED_SNPS)].copy()
    keep = ["analysis", "trait", "LOCUS_ID", "gene", "tissue", "modality",
            "PP3", "PP4", "nsnps"]
    for c in ("symbol", "module_id", "go_invisible"):
        if c in d.columns:
            keep.append(c)
    # The abf layer already ran p12 in {1e-6, 5e-6, 1e-5}, so fallback cells carry the same
    # prior-survival descriptor as SuSiE cells, read off those arms at the same SNP floor.
    keys = ["analysis", "trait", "LOCUS_ID", "gene", "tissue", "modality"]
    sw = p12_survival(abf[abf["nsnps"] >= MIN_SHARED_SNPS], keys)
    return d[keep].merge(sw, on=keys, how="left")


def scope_fallback(abf: pd.DataFrame, analyses_run: set[str],
                   max_snps: int = MAX_SNPS_PRIMARY) -> pd.DataFrame:
    """The abf cells a results root may fall back to.

    The primary grid keeps every analysis, including one whose GWAS fine-mapped nowhere:
    its cells all stand on abf, and that is a result (see `apply_hierarchy`). A sensitivity
    root holds only the analyses re-run at its guard. Merging in every OTHER analysis's abf
    cells would present them as scored in that arm -- all by fallback -- when they were
    never run there, and the arm's hierarchy counts would silently include them.
    """
    if int(max_snps) == MAX_SNPS_PRIMARY:
        return abf
    return abf[abf["analysis"].isin(analyses_run)]


def run_meta(arm: str = "switch", src: Path | None = None,
             require_gtex_match: bool = True,
             sqtl: str = "representative",
             max_snps: int = MAX_SNPS_PRIMARY) -> Path:
    src = results_dir(src or signal_root(arm, max_snps), sqtl)
    src.mkdir(parents=True, exist_ok=True)
    pairs = _load_susie(src)
    print(f"  loaded {len(pairs):,} coloc.susie signal-pair rows over "
          f"{pairs['analysis'].nunique()} analyses, {pairs['gene'].nunique():,} genes")
    pairs.to_parquet(src / "signal_pairs.parquet", index=False)

    # ---- cell-level SuSiE result, primary and sensitivity ------------------
    arms = {"gtex_matched": True, "all_signals": False}
    cells = {}
    for label, flag in arms.items():
        b = best_signal_pair(pairs, require_gtex_match=flag)
        b.to_parquet(src / f"cells_susie_{label}.parquet", index=False)
        cells[label] = b
        print(f"  {label:<13} {len(b):,} (cell, modality) rows, "
              f"{b['gene'].nunique() if len(b) else 0:,} genes")
    primary_label = "gtex_matched" if require_gtex_match else "all_signals"
    primary = cells[primary_label]
    primary.to_parquet(src / "cells_susie.parquet", index=False)

    # ---- estimator hierarchy: susie where it spoke, abf where it did not ---
    # A sensitivity root's analyses are exactly those stage A cached at its guard.
    gw = signal_root("switch", max_snps) / "gwas_susie"
    ran = {p.name for p in gw.iterdir() if p.is_dir()} if gw.exists() else set()
    abf = scope_fallback(_abf_cells(arm), ran, max_snps)
    merged = apply_hierarchy(primary, abf)
    # Evidence-strength descriptors on every cell (see FALLBACK_REASONS / PRIOR_ROBUSTNESS).
    merged["fallback_reason"] = fallback_reasons(
        merged, load_gwas_status(signal_root("switch", max_snps)), pairs,
        require_gtex_match=require_gtex_match)
    merged["prior_robustness"] = [prior_robustness(x) for x in merged["p12_min_call"]]
    fb = merged.loc[merged["estimator"] == "abf", "fallback_reason"].value_counts()
    print(f"  abf fallback reasons: {fb.to_dict()}")
    merged.to_parquet(src / "cells_hierarchy.parquet", index=False)
    n_susie = int((merged["estimator"] == "susie").sum())
    print(f"  hierarchy: {n_susie:,} / {len(merged):,} (cell, modality) rows scored by "
          f"coloc.susie; {len(merged) - n_susie:,} fell back to coloc.abf")

    # ---- paired modality contrast, on the hierarchy ------------------------
    # Reuses the abf layer's estimators verbatim so the signal-level number is directly
    # comparable to the published abf number rather than to a re-specified test.
    wide = merged.pivot_table(
        index=[c for c in ("analysis", "trait", "LOCUS_ID", "gene", "symbol",
                           "module_id", "go_invisible", "tissue") if c in merged.columns],
        columns="modality", values=["PP3", "PP4"], aggfunc="first")
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    need = ["PP4_sQTL", "PP4_eQTL", "PP3_sQTL", "PP3_eQTL"]
    for c in need:
        if c not in wide:
            wide[c] = np.nan
    wide = wide.dropna(subset=need)
    for k in ("sQTL", "eQTL"):
        den = wide[f"PP3_{k}"] + wide[f"PP4_{k}"]
        wide[f"cond_{k}"] = np.where(den > 0, wide[f"PP4_{k}"] / den, np.nan)
    wide["nsnps_sQTL"] = np.nan
    wide.to_parquet(src / "cells.parquet", index=False)

    genes = cmc.collapse_genes(wide, call=PP4_CALL)
    genes.to_parquet(src / "genes.parquet", index=False)

    rows = []
    for (analysis, trait), sub in genes.groupby(["analysis", "trait"]):
        rows.append({"analysis": analysis, "trait": trait,
                     **cmc.paired_tests(sub, call=PP4_CALL)})
    rows.append({"analysis": "POOLED", "trait": "ALL",
                 **cmc.paired_tests(genes, call=PP4_CALL)})
    contrast = pd.DataFrame(rows)
    contrast.to_parquet(src / "contrast.parquet", index=False)

    _write_report(src, arm, pairs, cells, merged, genes, contrast,
                  primary_label=primary_label, sqtl=sqtl)
    print(f"\n  wrote {src}")
    return src


def _fmt(v, nd=3):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "n/a"
    if isinstance(v, float):
        return f"{v:.{nd}g}" if (abs(v) < 1e-3 and v != 0) else f"{v:.{nd}f}"
    return str(v)


def _write_report(src: Path, arm: str, pairs: pd.DataFrame, cells: dict,
                  merged: pd.DataFrame, genes: pd.DataFrame,
                  contrast: pd.DataFrame, primary_label: str,
                  sqtl: str = "representative") -> None:
    L: list[str] = []
    A = L.append
    A("# Signal-level colocalization (coloc.susie)")
    A("")
    if "sensitivity" in src.parts:
        A(f"**SENSITIVITY ARM** (`{src}`): the GWAS SuSiE SNP guard was raised above the "
          f"primary {MAX_SNPS_PRIMARY:,} for the analyses in this directory only. Loci "
          "recovered this way were empirically enriched for GWAS-reference-LD "
          "inconsistency, and none may enter the biological narrative without a "
          "locus-specific LD audit (`08d.locus_ld_robustness.sh`). This is not the "
          "primary grid.")
        A("")
    A(f"Gene pool arm: `{arm}`. Primary signal filter: `{primary_label}`. "
      f"sQTL phenotypes: `{sqtl}`.")
    A("")
    if sqtl == "representative":
        A("Each gene contributes GTEx's single grouped-permutation representative intron, "
          "so a locus whose disease-relevant event is not that intron cannot colocalize "
          "here however real it is. The `all` arm is what tests a named event.")
    else:
        A("Every intron phenotype of each gene is tested, not only GTEx's grouped-"
          "permutation representative, so a named literature event is testable. The "
          "`coloc.abf` fallback rows still come from the representative-intron layer.")
        A("")
        A("A cell's PP4 here is a **maximum over that gene's introns**, and "
          "`n_phenotypes_tested` records how many it was taken over (1 in the "
          "representative arm, which is what makes the arms comparable). Read a PP4 "
          "against that count, not on its own; `phenotype_id` names the intron that won, "
          "which is what the event audit checks against the curated coordinates.")
    A("")
    A("`coloc.susie` fine-maps both traits and colocalizes credible set against "
      "credible set, so a locus carrying more than one causal signal is not forced "
      "into the single-causal-variant assumption `coloc.abf` makes. The QTL side is "
      "re-fit with `susie_rss` because GTEx v11 ships credible-set summaries, not "
      "SuSiE objects.")
    A("")
    A("## What was fit")
    A("")
    A(f"- signal-pair posteriors: {len(pairs):,}")
    for label, b in cells.items():
        A(f"- cells surviving `{label}`: {len(b):,} "
          f"({b['gene'].nunique() if len(b) else 0:,} genes)")
    n_s = int((merged['estimator'] == 'susie').sum())
    A(f"- estimator hierarchy: **{n_s:,} coloc.susie**, "
      f"{len(merged) - n_s:,} coloc.abf fallback")
    if "fallback_reason" in merged.columns:
        A("  - why coloc.abf, first place the cell left the signal-level pipeline:")
        fb = merged.loc[merged["estimator"] == "abf", "fallback_reason"].value_counts()
        for k, v in fb.items():
            A(f"    - `{k}`: {v:,}")
    if "prior_robustness" in merged.columns:
        called = merged[merged["PP4"] >= PP4_CALL]
        pr = called.groupby(["estimator", "prior_robustness"]).size()
        A(f"- prior robustness of the {len(called):,} cells calling at the primary prior: "
          + ", ".join(f"{e} `{lb}` {n:,}" for (e, lb), n in pr.items()))
    A("")
    A("## Reference-LD caveat")
    A("")
    A("The QTL side is fine-mapped under the 1000G EUR Phase 3 panel, because that is "
      "the only LD available for GTEx. GTEx brain donors are not a pure EUR sample, so "
      "a re-fit credible set can be a reference-LD artefact. The primary arm therefore "
      "keeps only QTL signals that share a variant with **GTEx's own** shipped credible "
      "set (`cs_matches_gtex`); the `all_signals` arm drops that requirement and is "
      "reported as a sensitivity, never as the headline. `kriging_rss` diagnostics are "
      "in `signal_pairs.parquet`.")
    A("")
    A("## Paired modality contrast")
    A("")
    A("| analysis | trait | genes | sQTL coloc | eQTL coloc | splicing-only | "
      "expression-only | McNemar P |")
    A("|---|---|---|---|---|---|---|---|")
    for r in contrast.itertuples(index=False):
        A(f"| {r.analysis} | {r.trait} | {r.n_genes} | {r.n_sqtl_coloc} | "
          f"{r.n_eqtl_coloc} | {r.splicing_only} | {r.expression_only} | "
          f"{_fmt(r.mcnemar_p)} |")
    A("")
    A("Conditioning on `PP4_sQTL >= 0.8` and then reading `PP4_eQTL` is a **selection**, "
      "so a splicing-preferential count is a set of locus nominations, not an unbiased "
      "splicing-specificity estimate. The unbiased test is the paired McNemar / Wilcoxon "
      "above.")
    A("")
    A("## Per-tissue consistency")
    A("")
    A("`genes.parquet` carries `n_tissue_sQTL_coloc`, `frac_tissue_sQTL_coloc`, "
      "`max_tissue_sQTL` and the full `tissue_pp4_sQTL` vector. The headline PP4 is a "
      "maximum over 13 correlated tissues; quote it with the consistency count beside "
      "it, never alone.")
    A("")
    top = genes.sort_values("PP4_sQTL", ascending=False).head(25)
    A("| gene | trait | PP4 sQTL | PP4 eQTL | tissues coloc | max tissue |")
    A("|---|---|---|---|---|---|")
    for r in top.itertuples(index=False):
        sym = getattr(r, "symbol", r.gene)
        A(f"| {sym} | {r.trait} | {_fmt(r.PP4_sQTL)} | {_fmt(r.PP4_eQTL)} | "
          f"{int(r.n_tissue_sQTL_coloc)}/{int(r.n_tissue)} | {r.max_tissue_sQTL} |")
    A("")
    (src / "COLOC_SIGNAL_SUSIE.md").write_text("\n".join(L) + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("prep", "meta"), required=True)
    ap.add_argument("--arm", choices=ARMS_GENE_POOL, default="switch",
                    help="gene pool; mirrors coloc_modality_contrast's arms")
    ap.add_argument("--all-signals", action="store_true",
                    help="do not require re-fit QTL credible sets to agree with GTEx's "
                         "own (reported as a sensitivity arm, never as primary)")
    ap.add_argument("--sqtl", choices=("representative", "all"), default="representative",
                    help="which sQTL phenotypes the shards were produced from; must match "
                         "COLOC_SIGNAL_SQTL, and selects the results subdirectory")
    ap.add_argument("--max-snps", type=int, default=MAX_SNPS_PRIMARY,
                    help="GWAS SuSiE SNP guard the shards were produced at; must match "
                         "COLOC_GWAS_MAX_SNPS. Anything but the primary "
                         f"{MAX_SNPS_PRIMARY} reads a sensitivity root")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)

    if args.stage == "prep":
        if args.max_snps != MAX_SNPS_PRIMARY:
            raise SystemExit("--max-snps does not apply to --stage prep: the target grid "
                             "is shared by every guard")
        run_prep(arm=args.arm, dest=args.out)
    else:
        run_meta(arm=args.arm, src=args.out,
                 require_gtex_match=not args.all_signals, sqtl=args.sqtl,
                 max_snps=args.max_snps)


if __name__ == "__main__":
    main()
