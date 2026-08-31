"""Per-gene sQTL-vs-eQTL colocalization contrast (coloc.abf on GTEx v11 all-pairs).

WHY THIS EXISTS
---------------
The splicing-specificity headline is currently a *set-level* statement: a ratio of two
odds ratios (sQTL-OR / eQTL-OR) computed over gene sets and meta-analysed across brain
xQTL analyses (`qtl_anchoring_meta`). That estimator has two properties a reviewer will
press on. Both of its component ORs are below 1 in every `all_modules` analysis, so the
contrast is a ratio of two depletions rather than an sQTL enrichment; and because the
sQTL and eQTL arms are computed over gene *sets*, nothing ties a splicing signal and an
expression signal to the same gene at the same locus.

This module puts the same claim on per-gene footing. For one switch gene, at one GWAS
locus, in one GTEx brain tissue, it asks whether the disease association colocalizes
with that gene's *splicing* QTL more than with that gene's *expression* QTL. Because
the comparison is made within a gene, whatever selected the gene into the analysis
cancels: the IsoGraph module membership that chose the gene is identical for both arms.

Only now possible because GTEx v11 cis all-pairs (nominal stats for every variant, not
just the significant ones) is on disk. eCAVIAR CLPP -- the estimator used by
10.coloc_clpp.R -- consumes credible sets, which is what GTEx used to ship, and it
requires fine-mapping resolution on BOTH sides, which is why its yield is 2-37 genes
per trait. coloc.abf needs the full cis window and has no such requirement.

WHAT MAKES THE PAIRING VALID
----------------------------
* Same donors. GTEx sQTL and eQTL for a tissue are called in the same samples, so the
  two arms share sample size, ancestry, covariates and genotypes exactly.
* Same GWAS. Both arms colocalize against the identical per-locus GWAS summary stats.
* Same variants. Both arms are restricted to the variants shared between the QTL cis
  window and the locus GWAS, so neither arm is scored on SNPs the other lacked.
* Same gene selection. Module membership enters both arms identically.

THREE ASYMMETRIES THAT WOULD OTHERWISE BIAS IT, AND HOW EACH IS HANDLED
----------------------------------------------------------------------
1. sQTL genes have MANY phenotypes (introns); eQTL genes have one. Taking the best
   intron by colocalization would hand the sQTL arm a free maximum over ~7 tests.
   Handled by using GTEx's OWN grouped-permutation representative phenotype per gene
   (the `phenotype_id` in `*.sGenes.txt.gz`, one row per gene, chosen by group-corrected
   `pval_beta`). That selection is made from the QTL data alone and never sees the GWAS,
   so it cannot be circular. `--sqtl-phenotypes all` runs the max-over-introns arm as a
   sensitivity so the size of the advantage is measured rather than assumed.
2. eQTLs are better powered than sQTLs in GTEx, so a raw PP4 comparison would favour
   expression for reasons that have nothing to do with disease. Handled by making the
   PRIMARY continuous statistic the conditional posterior PP4/(PP3+PP4) -- the
   probability that the two traits share a causal variant GIVEN that each has one in
   the window. That conditioning is the standard remedy for exactly this asymmetry.
   Unconditional PP4 is reported alongside it, never instead of it.
3. Aggregating over 13 brain tissues by maximum inflates PP4. It inflates it for both
   arms identically (same 13 tissues), and a cell enters the paired analysis only when
   BOTH modalities cleared the shared-SNP floor, so the maximum is taken over the same
   set of cells on both sides.

WHY THE BUILD MISMATCH DOES NOT MATTER HERE
-------------------------------------------
GTEx is GRCh38 and PGC3 SCZ is GRCh37, bridged by rsID through the GTEx v8 WGS lookup
(98.6% of v11 variants carry an rsID there, measured on chr22). Two reasons this is not
a threat. coloc.abf is LD-free and its approximate Bayes factors depend on z^2, so it is
invariant to allele orientation -- unlike the CLPP path, it needs no LD panel and
therefore no build-matched reference. And any variant the bridge drops is dropped from
the sQTL and eQTL arms of the same gene identically, so bridge incompleteness costs
power and cannot bias the contrast. (Allele harmonization is still performed, but only
to report the DIRECTION of effect, which is not part of the primary test.)

STAGES
------
  --stage bridge   one-off: GTEx v8 WGS lookup -> per-chromosome variant_id/rsID parquet
                   cache under inputs/raw/gtex_v11/variant_bridge/ (gitignored).
  --stage prep     per-analysis (gene, locus) targets + the GTEx representative sQTL
                   phenotype per (tissue, gene). Cheap; login-node safe.
  [R stage]        05_genetic_anchoring/_h/20.coloc_modality_abf.R <tissue> -- extracts
                   cis all-pairs, joins the locus GWAS, runs coloc::coloc.abf.
  --stage meta     the paired statistics, sensitivity arms and the report.

Outputs under 05_genetic_anchoring/_m/coloc_modality_contrast/.
"""
from __future__ import annotations

import argparse
import gzip
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.paths import ensure_dir, rel, stage_out

# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #
GTEX_V11 = Path("/ocean/projects/bio250020p/shared/resources/public-data/gtex_v11")
V8_LOOKUP = Path(
    "/ocean/projects/bio250020p/shared/resources/public-data/gtex-v8/"
    "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"
)
BRIDGE_DIR = rel("inputs", "raw", "gtex_v11", "variant_bridge")

# The 13 GTEx brain tissues, mirroring qtl_anchoring._GTEX_TISSUE. A switch gene may
# colocalize in any brain region, so all 13 are scanned for discovery power; the
# maximum over them is taken symmetrically for both modalities.
GTEX_BRAIN: tuple[str, ...] = (
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
)

# coloc analysis dir -> trait key. These are exactly the analyses the CLPP layer and
# module_coloc_convergence already report, so the per-gene contrast lands on the same
# cells as the existing set-level result.
# NOTE: `brainseq-sczd` is a legacy alias of `brainseq-sczd__scz` (identical
# clpp_results.tsv, 378 rows / 74 genes); only the trait-suffixed form is registered.
ANALYSES: dict[str, str] = {
    "aging__ad": "ad",
    "aging__pd": "pd",
    "aging__lbd": "lbd",
    "aging__als": "als",
    "aging__scz": "scz",
    "brainseq-sczd__scz": "scz",
}

# Pre-specified thresholds. Fixed before any contrast was computed.
MIN_SHARED_SNPS = 100     # a coloc window needs a window; below this PP is noise
PP4_CALL = 0.80           # colocalized (coloc convention)
PP4_CALL_LOOSE = 0.50     # sensitivity arm
P12_PRIORS = (1e-5, 5e-6, 1e-6)   # 1e-5 = coloc default = PRIMARY
P12_PRIMARY = 1e-5
_SEED = 13


def _bare(s: pd.Series) -> pd.Series:
    return s.astype(str).str.split(".", n=1).str[0]


def out_dir() -> Path:
    return ensure_dir(stage_out("anchoring", "coloc_modality_contrast"))


def coloc_dir(analysis: str) -> Path:
    return stage_out("anchoring", "coloc", analysis)


# --------------------------------------------------------------------------- #
# Stage: bridge
# --------------------------------------------------------------------------- #
def build_variant_bridge(lookup: Path = V8_LOOKUP, dest: Path | None = None,
                         force: bool = False) -> Path:
    """GTEx b38 `variant_id` -> dbSNP rsID, cached as one parquet per chromosome.

    GTEx v11 ships no lookup table of its own, so the v8 WGS freeze lookup is used;
    positions and alleles are identical across v8/v11 (both GRCh38). Measured coverage
    of v11 variants on chr22 is 98.6%. The dropped 1.4% are v11-only variants, and they
    are dropped from the sQTL and eQTL arms of a gene identically -- see the module
    docstring on why that costs power but cannot bias the contrast.
    """
    dest = ensure_dir(dest or BRIDGE_DIR)
    done = dest / "_COMPLETE"
    if done.exists() and not force:
        print(f"  bridge cache present ({dest}); use --force to rebuild")
        return dest
    if not lookup.exists():
        raise SystemExit(f"GTEx v8 lookup not found: {lookup}")

    # Stream to per-chromosome temp files rather than accumulating in memory: the
    # lookup is 46.6M rows and a Python-side dict of them would need >8 GB. awk writes
    # each chromosome as it goes, then each (2-4M row) chromosome is converted alone.
    print(f"  streaming {lookup.name} -> {dest}")
    tmp = ensure_dir(dest / "_tmp")
    for stale in tmp.iterdir():
        stale.unlink()
    prog = (
        'FNR==1{next} '
        '($7!="." && $7!="") {print $1"\\t"$7 > (OUT "/" $2 ".tsv")}'
    )
    zcat = subprocess.Popen(["zcat", str(lookup)], stdout=subprocess.PIPE)
    awk = subprocess.Popen(["awk", "-F", "\t", "-v", f"OUT={tmp}", prog],
                           stdin=zcat.stdout)
    zcat.stdout.close()
    awk.wait()
    zcat.wait()
    if awk.returncode or zcat.returncode:
        raise RuntimeError("variant-bridge stream failed")

    n = 0
    chroms = []
    for f in sorted(tmp.iterdir()):
        if f.suffix != ".tsv":
            continue
        df = pd.read_csv(f, sep="\t", header=None, names=["variant_id", "rsid"],
                         dtype=str)
        df.to_parquet(dest / f"{f.stem}.parquet", index=False)
        n += len(df)
        chroms.append(f.stem)
        f.unlink()
    tmp.rmdir()
    done.write_text(f"rows={n}\nchroms={len(chroms)}\nlookup={lookup}\n")
    print(f"  wrote {len(chroms)} chromosome parquets, {n:,} variant->rsID rows")
    return dest


# --------------------------------------------------------------------------- #
# Stage: prep
# --------------------------------------------------------------------------- #
def load_targets(analysis: str) -> pd.DataFrame:
    """(LOCUS_ID, gene) pairs to colocalize, from the analysis' testable loci.

    Uses the same locus->gene membership as 10.coloc_clpp.R, so the per-gene contrast
    is computed on exactly the gene/locus pairs the CLPP layer already reports on.
    """
    d = coloc_dir(analysis)
    loci = pd.read_csv(d / "susie" / "loci_testable.tsv", sep="\t")
    rows = []
    for r in loci.itertuples(index=False):
        for g in str(r.genes).split(","):
            g = g.strip()
            if g:
                rows.append((analysis, ANALYSES[analysis], r.LOCUS_ID, int(r.chr), g))
    t = pd.DataFrame(rows, columns=["analysis", "trait", "LOCUS_ID", "chr", "gene"])
    t = t.drop_duplicates()

    # Carry the module tag so a hit can be attributed to a GO-invisible module, matching
    # the CLPP table's `go_invisible` column.
    cand = pd.read_csv(d / "candidate_genes.tsv", sep="\t",
                       usecols=["gene", "module_id", "go_invisible", "symbol"])
    cand = cand.drop_duplicates("gene")
    return t.merge(cand, on="gene", how="left")


def sqtl_representative(tissues=GTEX_BRAIN, genes: set[str] | None = None) -> pd.DataFrame:
    """GTEx's own grouped-permutation representative intron per (tissue, gene).

    `*.sGenes.txt.gz` carries ONE row per gene (mean group size ~7.4 introns), with the
    representative `phenotype_id` chosen by the group-corrected permutation p-value.
    That choice is made from genotype and splicing data alone -- the GWAS plays no part
    -- which is what makes using it as the single sQTL phenotype per gene non-circular.
    """
    frames = []
    for tissue in tissues:
        f = GTEX_V11 / "GTEx_Analysis_v11_sQTL" / f"{tissue}.v11.sGenes.txt.gz"
        if not f.exists():
            raise SystemExit(f"missing sGenes file: {f}")
        d = pd.read_csv(f, sep="\t", usecols=["phenotype_id", "gene_id", "pval_beta",
                                              "num_var", "group_size", "qval"])
        d["gene"] = _bare(d["gene_id"])
        if genes is not None:
            d = d[d["gene"].isin(genes)]
        d.insert(0, "tissue", tissue)
        frames.append(d.drop(columns=["gene_id"]))
    return pd.concat(frames, ignore_index=True)


# Case counts for the two traits whose sumstats carry no per-SNP N. Both are the
# published study totals and both reconcile exactly with `fixed_n` in the trait
# registry, which is the check that they are the right numbers:
#   LBD  Chia et al. 2021   2,591 cases +   4,027 controls =   6,618 = fixed_n
#   ALS  van Rheenen 2021  27,205 cases + 110,881 controls = 138,086 = fixed_n
_FIXED_CASE_N = {"lbd": 2591, "als": 27205}


def gwas_meta(traits=("scz", "ad", "pd", "lbd", "als"), scan_rows: int = 200_000
              ) -> pd.DataFrame:
    """Per-trait GWAS N and case fraction `s`, derived from the sumstats themselves.

    Recorded for provenance only. `coloc.abf` computes its Bayes factors from beta and
    varbeta, and for a case-control dataset supplied that way the posteriors do not
    depend on N or s at all -- verified in tests/test_coloc_modality_contrast.py, which
    fits the same data at s=0.05/N=1e3 and s=0.30/N=5e4 and asserts identical summaries.
    Derived rather than hardcoded because a hand-entered case fraction is exactly the
    kind of number that goes stale silently when a trait's sumstats are swapped.
    """
    from isograph_benchmark.real_data import gwas_traits as gt
    rows = []
    for key in traits:
        spec = gt.TRAITS[key]
        with gzip.open(spec.path, "rt") as fh:
            for line in fh:
                if spec.comment and line.startswith(spec.comment):
                    continue
                hdr = line.rstrip("\n").split("\t")
                break
            if spec.n_cas_col and spec.n_con_col:
                i, j = hdr.index(spec.n_cas_col), hdr.index(spec.n_con_col)
                cas = con = 0.0
                for k, line in enumerate(fh):
                    if k >= scan_rows:
                        break
                    f = line.split("\t")
                    try:
                        cas = max(cas, float(f[i]))
                        con = max(con, float(f[j]))
                    except (ValueError, IndexError):
                        continue
                n_cas, n_con, source = cas, con, "per-SNP max"
            else:
                n_cas = float(_FIXED_CASE_N[key])
                n_con = float(spec.fixed_n) - n_cas
                source = "published totals (registry fixed_n)"
        rows.append({"trait": key, "n_cas": int(n_cas), "n_con": int(n_con),
                     "n_gwas": int(n_cas + n_con), "s": n_cas / (n_cas + n_con),
                     "source": source})
    return pd.DataFrame(rows)


def run_prep(analyses=tuple(ANALYSES), dest: Path | None = None) -> Path:
    dest = ensure_dir(dest or out_dir())
    tabs = [load_targets(a) for a in analyses]
    targets = pd.concat(tabs, ignore_index=True)
    targets.to_parquet(dest / "targets.parquet", index=False)

    genes = set(targets["gene"])
    rep = sqtl_representative(genes=genes)
    rep.to_parquet(dest / "sqtl_representative.parquet", index=False)

    gm = gwas_meta()
    gm.to_csv(dest / "gwas_meta.tsv", sep="\t", index=False)

    print(f"  {len(targets):,} (analysis, locus, gene) targets over "
          f"{targets['analysis'].nunique()} analyses, {targets['gene'].nunique():,} genes")
    print(f"  sQTL representative phenotypes: {len(rep):,} rows over "
          f"{rep['tissue'].nunique()} tissues, {rep['gene'].nunique():,} genes")
    by = targets.groupby("analysis").agg(loci=("LOCUS_ID", "nunique"),
                                         genes=("gene", "nunique"),
                                         pairs=("gene", "size"))
    print(by.to_string())
    print("\n  GWAS metadata (provenance only; coloc.abf posteriors do not use N/s):")
    print(gm.to_string(index=False))
    return dest


# --------------------------------------------------------------------------- #
# Stage: meta -- the paired statistics
# --------------------------------------------------------------------------- #
def _load_abf(src: Path) -> pd.DataFrame:
    abf = src / "abf"
    if not abf.exists():
        raise SystemExit(f"no abf/ results in {src}; run 20.coloc_modality_abf.sh first")
    parts = [pd.read_parquet(p) for p in sorted(abf.iterdir()) if p.suffix == ".parquet"]
    if not parts:
        raise SystemExit(f"no parquet files under {abf}")
    return pd.concat(parts, ignore_index=True)


def pair_cells(abf: pd.DataFrame, p12: float = P12_PRIMARY,
               min_shared: int = MIN_SHARED_SNPS) -> pd.DataFrame:
    """Wide (analysis, locus, gene, tissue) cells carrying BOTH modalities.

    The `both testable` restriction is what makes this a paired design: a cell is kept
    only when the sQTL and eQTL arms each cleared the shared-SNP floor, so neither arm
    is credited with a comparison the other could not have made.
    """
    d = abf[(abf["p12"] == p12) & (abf["nsnps"] >= min_shared)].copy()
    keys = ["analysis", "trait", "LOCUS_ID", "gene", "symbol", "module_id",
            "go_invisible", "tissue"]
    wide = d.pivot_table(index=keys, columns="modality",
                         values=["PP3", "PP4", "nsnps"], aggfunc="first")
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    need = ["PP4_sQTL", "PP4_eQTL", "PP3_sQTL", "PP3_eQTL"]
    for c in need:
        if c not in wide:
            wide[c] = np.nan
    wide = wide.dropna(subset=need)
    # Conditional posterior: P(shared causal variant | each trait has one in the
    # window). This is the primary continuous statistic -- it divides out the QTL
    # discovery-power difference between splicing and expression.
    for k in ("sQTL", "eQTL"):
        denom = wide[f"PP3_{k}"] + wide[f"PP4_{k}"]
        wide[f"cond_{k}"] = np.where(denom > 0, wide[f"PP4_{k}"] / denom, np.nan)
    return wide


def collapse_genes(wide: pd.DataFrame) -> pd.DataFrame:
    """Best tissue per (analysis, locus, gene), symmetrically for both modalities."""
    keys = ["analysis", "trait", "LOCUS_ID", "gene", "symbol", "module_id", "go_invisible"]
    agg = wide.groupby(keys, dropna=False).agg(
        n_tissue=("tissue", "nunique"),
        PP4_sQTL=("PP4_sQTL", "max"), PP4_eQTL=("PP4_eQTL", "max"),
        cond_sQTL=("cond_sQTL", "max"), cond_eQTL=("cond_eQTL", "max"),
        nsnps=("nsnps_sQTL", "median"),
    ).reset_index()
    return agg


def paired_tests(genes: pd.DataFrame, call: float = PP4_CALL) -> dict:
    """Primary binary (exact McNemar) + primary continuous (paired Wilcoxon)."""
    s_hit = genes["PP4_sQTL"] >= call
    e_hit = genes["PP4_eQTL"] >= call
    b = int((s_hit & ~e_hit).sum())     # splicing-only
    c = int((~s_hit & e_hit).sum())     # expression-only
    both = int((s_hit & e_hit).sum())
    neither = int((~s_hit & ~e_hit).sum())
    # Exact McNemar: binomial on the discordant pairs against 0.5.
    mc_p = stats.binomtest(b, b + c, 0.5).pvalue if (b + c) else np.nan

    d = genes.dropna(subset=["cond_sQTL", "cond_eQTL"])
    diff = d["cond_sQTL"] - d["cond_eQTL"]
    nz = diff[diff != 0]
    if len(nz) >= 6:
        _, w_p = stats.wilcoxon(d["cond_sQTL"], d["cond_eQTL"],
                                zero_method="wilcox", alternative="two-sided")
    else:
        w_p = np.nan
    # Unconditional PP4, reported alongside so the power asymmetry stays visible.
    du = genes.dropna(subset=["PP4_sQTL", "PP4_eQTL"])
    diffu = du["PP4_sQTL"] - du["PP4_eQTL"]
    if (diffu != 0).sum() >= 6:
        _, wu_p = stats.wilcoxon(du["PP4_sQTL"], du["PP4_eQTL"],
                                 zero_method="wilcox", alternative="two-sided")
    else:
        wu_p = np.nan
    return {
        "n_genes": int(len(genes)),
        "n_sqtl_coloc": int(s_hit.sum()), "n_eqtl_coloc": int(e_hit.sum()),
        "both": both, "splicing_only": b, "expression_only": c, "neither": neither,
        "mcnemar_p": mc_p,
        "splicing_share_discordant": (b / (b + c)) if (b + c) else np.nan,
        "median_cond_sQTL": float(d["cond_sQTL"].median()) if len(d) else np.nan,
        "median_cond_eQTL": float(d["cond_eQTL"].median()) if len(d) else np.nan,
        "median_cond_diff": float(nz.median()) if len(nz) else np.nan,
        "wilcoxon_cond_p": w_p,
        "wilcoxon_pp4_p": wu_p,
        "median_pp4_sQTL": float(du["PP4_sQTL"].median()) if len(du) else np.nan,
        "median_pp4_eQTL": float(du["PP4_eQTL"].median()) if len(du) else np.nan,
    }


def _bh(p: pd.Series) -> pd.Series:
    ok = p.notna()
    q = pd.Series(np.nan, index=p.index)
    if ok.sum():
        from statsmodels.stats.multitest import multipletests
        q.loc[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


# Sensitivity arms. `primary` is fixed in advance; every other row exists so a reviewer
# can see the claim move (or not) under the choices that could plausibly drive it.
ARMS: tuple[tuple[str, dict], ...] = (
    ("primary",           {"p12": 1e-5, "call": PP4_CALL}),
    ("p12_5e-6",          {"p12": 5e-6, "call": PP4_CALL}),
    ("p12_1e-6",          {"p12": 1e-6, "call": PP4_CALL}),
    ("pp4_call_0.5",      {"p12": 1e-5, "call": PP4_CALL_LOOSE}),
    ("min_shared_500",    {"p12": 1e-5, "call": PP4_CALL, "min_shared": 500}),
)


def run_meta(src: Path | None = None) -> Path:
    src = src or out_dir()
    abf = _load_abf(src)
    print(f"  loaded {len(abf):,} coloc.abf rows "
          f"({abf['tissue'].nunique()} tissues, {abf['analysis'].nunique()} analyses)")

    # ---- per-cell and per-gene tables for the primary arm ------------------
    wide = pair_cells(abf)
    genes = collapse_genes(wide)
    wide.to_parquet(src / "cells.parquet", index=False)
    genes.to_parquet(src / "genes.parquet", index=False)

    # ---- the contrast, per analysis x arm ----------------------------------
    rows = []
    for arm, kw in ARMS:
        w = pair_cells(abf, p12=kw["p12"], min_shared=kw.get("min_shared", MIN_SHARED_SNPS))
        g = collapse_genes(w)
        for (analysis, trait), sub in g.groupby(["analysis", "trait"]):
            rows.append({"arm": arm, "analysis": analysis, "trait": trait,
                         "p12": kw["p12"], "pp4_call": kw["call"],
                         "min_shared": kw.get("min_shared", MIN_SHARED_SNPS),
                         **paired_tests(sub, call=kw["call"])})
        # pooled across analyses, treating each (analysis, gene, locus) as one unit
        rows.append({"arm": arm, "analysis": "POOLED", "trait": "ALL",
                     "p12": kw["p12"], "pp4_call": kw["call"],
                     "min_shared": kw.get("min_shared", MIN_SHARED_SNPS),
                     **paired_tests(g, call=kw["call"])})
    contrast = pd.DataFrame(rows)
    # BH across the per-analysis cells of the primary arm only; the pooled row and the
    # sensitivity arms are not independent tests and are excluded from the correction.
    prim = (contrast["arm"] == "primary") & (contrast["analysis"] != "POOLED")
    contrast["mcnemar_q"] = np.nan
    contrast.loc[prim, "mcnemar_q"] = _bh(contrast.loc[prim, "mcnemar_p"])
    contrast.loc[prim, "wilcoxon_cond_q"] = _bh(contrast.loc[prim, "wilcoxon_cond_p"])
    contrast.to_parquet(src / "contrast.parquet", index=False)

    # ---- GO-invisible split: does the per-gene contrast localise? ----------
    go_rows = []
    for (analysis, trait, go), sub in genes.groupby(
            ["analysis", "trait", "go_invisible"], dropna=False):
        go_rows.append({"analysis": analysis, "trait": trait,
                        "go_invisible": go, **paired_tests(sub)})
    go = pd.DataFrame(go_rows)
    go.to_parquet(src / "go_split.parquet", index=False)

    _write_report(src, abf, wide, genes, contrast, go)
    _print_summary(contrast)
    return src


def _print_summary(contrast: pd.DataFrame) -> None:
    p = contrast[contrast["arm"] == "primary"]
    cols = ["analysis", "trait", "n_genes", "n_sqtl_coloc", "n_eqtl_coloc",
            "splicing_only", "expression_only", "mcnemar_p", "mcnemar_q",
            "median_cond_sQTL", "median_cond_eQTL", "wilcoxon_cond_p"]
    print("\n=== PRIMARY ARM (p12=1e-5, PP4>=0.8, best tissue, GTEx representative intron) ===")
    print(p[cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))


def _fmt(v, nd=3):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "n/a"
    return f"{v:.{nd}g}"


def _write_report(src: Path, abf: pd.DataFrame, wide: pd.DataFrame,
                  genes: pd.DataFrame, contrast: pd.DataFrame, go: pd.DataFrame) -> None:
    prim = contrast[(contrast["arm"] == "primary") & (contrast["analysis"] != "POOLED")]
    pool = contrast[(contrast["arm"] == "primary") & (contrast["analysis"] == "POOLED")]
    L: list[str] = []
    A = L.append
    A("# Per-gene sQTL-vs-eQTL colocalization contrast")
    A("")
    A("Does a disease association colocalize with a switch gene's **splicing** QTL more")
    A("than with the **same gene's expression** QTL, at the same locus, in the same GTEx")
    A("brain tissue? This puts the splicing-specificity claim on per-gene footing; the")
    A("`qtl_anchoring_meta` contrast is a set-level ratio of two odds ratios.")
    A("")
    A("Estimator: `coloc::coloc.abf` on GTEx v11 cis **all-pairs** nominal statistics")
    A("(not credible sets), against the per-locus GWAS summary statistics already built")
    A("for the CLPP layer. Because the comparison is made within a gene, the IsoGraph")
    A("module membership that selected the gene cancels between the two arms.")
    A("")
    A("## Primary result")
    A("")
    A("Primary arm: `p12 = 1e-5` (coloc default), colocalized at `PP4 >= 0.80`, best of")
    A("the 13 GTEx brain tissues, GTEx's own grouped-permutation representative intron,")
    A(f"cells with `>= {MIN_SHARED_SNPS}` shared SNPs in **both** modalities.")
    A("")
    A("| analysis | trait | genes | sQTL coloc | eQTL coloc | splicing-only | expression-only | McNemar P | q | median cond PP4 sQTL | eQTL | Wilcoxon P |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in prim.itertuples(index=False):
        A(f"| {r.analysis} | {r.trait} | {r.n_genes} | {r.n_sqtl_coloc} | {r.n_eqtl_coloc} "
          f"| {r.splicing_only} | {r.expression_only} | {_fmt(r.mcnemar_p)} | {_fmt(r.mcnemar_q)} "
          f"| {_fmt(r.median_cond_sQTL)} | {_fmt(r.median_cond_eQTL)} | {_fmt(r.wilcoxon_cond_p)} |")
    if len(pool):
        r = pool.iloc[0]
        A(f"| **POOLED** | ALL | {r.n_genes} | {r.n_sqtl_coloc} | {r.n_eqtl_coloc} "
          f"| {r.splicing_only} | {r.expression_only} | {_fmt(r.mcnemar_p)} | - "
          f"| {_fmt(r.median_cond_sQTL)} | {_fmt(r.median_cond_eQTL)} | {_fmt(r.wilcoxon_cond_p)} |")
    A("")
    A("`splicing-only` / `expression-only` are the **discordant** genes -- the ones the")
    A("McNemar test is computed on. Concordant genes (both or neither) carry no")
    A("information about which modality colocalizes and are excluded by construction.")
    A("")
    A("## Two statistics, and why both are reported")
    A("")
    A("* **Binary (McNemar).** Counts genes where exactly one modality colocalizes. It is")
    A("  the statistic a reader can check against the gene list.")
    A("* **Continuous (paired Wilcoxon on `PP4/(PP3+PP4)`).** The conditional posterior")
    A("  asks whether the two traits share a causal variant *given that each has one in")
    A("  the window*. eQTLs are better powered than sQTLs in GTEx, so unconditional PP4")
    A("  would favour expression for reasons unrelated to disease; conditioning divides")
    A("  that difference out. Unconditional PP4 is in `contrast.parquet` as")
    A("  `wilcoxon_pp4_p` / `median_pp4_*` so the size of the asymmetry stays visible.")
    A("")
    A("## Sensitivity arms")
    A("")
    A("| arm | p12 | PP4 call | min shared SNPs | pooled genes | splicing-only | expression-only | McNemar P |")
    A("|---|---|---|---|---|---|---|---|")
    for arm, _ in ARMS:
        r = contrast[(contrast["arm"] == arm) & (contrast["analysis"] == "POOLED")]
        if not len(r):
            continue
        r = r.iloc[0]
        A(f"| {arm} | {r.p12:g} | {r.pp4_call:g} | {int(r.min_shared)} | {r.n_genes} "
          f"| {r.splicing_only} | {r.expression_only} | {_fmt(r.mcnemar_p)} |")
    A("")
    A("## GO-invisible split")
    A("")
    A("Whether the per-gene contrast localises to the GO-invisible modules. The")
    A("set-level `qtl_anchoring` result does **not** localise there after the")
    A("2026-08-29 refresh, so this is a check, not a confirmation.")
    A("")
    A("| analysis | GO-invisible | genes | splicing-only | expression-only | McNemar P |")
    A("|---|---|---|---|---|---|")
    for r in go.itertuples(index=False):
        A(f"| {r.analysis} | {r.go_invisible} | {r.n_genes} | {r.splicing_only} "
          f"| {r.expression_only} | {_fmt(r.mcnemar_p)} |")
    A("")
    A("## Scope and limits")
    A("")
    A(f"* {len(abf):,} `coloc.abf` fits; {len(wide):,} paired (gene, locus, tissue) cells;")
    A(f"  {len(genes):,} (gene, locus) pairs after collapsing tissues by maximum.")
    A("* `coloc.abf` assumes a **single causal variant** per trait per window. Where two")
    A("  independent causal variants sit in one window it under-calls sharing, for both")
    A("  arms alike. The CLPP layer (`10.coloc_clpp.R`) does not make this assumption and")
    A("  remains the estimator of record for *whether* a gene colocalizes at all.")
    A("* Only loci passing the GWAS window threshold enter, so this says nothing about")
    A("  sub-threshold signal.")
    A("* The rsID bridge (GTEx v8 WGS lookup) covers 98.6% of v11 variants. Dropped")
    A("  variants are dropped from both arms of a gene identically.")
    A("")
    A("Regenerate: `05_genetic_anchoring/_h/19.coloc_modality_prep.sh`,")
    A("`20.coloc_modality_abf.sh` (array over 13 tissues), `21.coloc_modality_meta.sh`.")
    (src / "COLOC_MODALITY_CONTRAST.md").write_text("\n".join(L) + "\n")
    print(f"  report -> {src / 'COLOC_MODALITY_CONTRAST.md'}")


# --------------------------------------------------------------------------- #
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("bridge", "prep", "meta", "all"), default="all")
    ap.add_argument("--force", action="store_true",
                    help="rebuild the variant bridge cache even if present")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)

    dest = a.out or out_dir()
    if a.stage in ("bridge", "all"):
        print("[bridge] GTEx variant_id -> rsID cache")
        build_variant_bridge(force=a.force)
    if a.stage in ("prep", "all"):
        print("[prep] coloc targets + GTEx representative sQTL phenotypes")
        run_prep(dest=dest)
    if a.stage == "meta":
        print("[meta] paired sQTL-vs-eQTL contrast")
        run_meta(src=dest)


if __name__ == "__main__":
    main()
