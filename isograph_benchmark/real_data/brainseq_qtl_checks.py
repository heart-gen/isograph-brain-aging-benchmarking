"""Validation checks on the BrainSEQ switch-QTL / abundance-QTL mapping.

Two checks the validation plan requires before any BrainSEQ QTL result is read. Both run on
the outputs of `brainseq_switch_qtl` and never modify them.

  --stage signpin   `S_g` is PC1 of a within-gene CLR composition, and PC1 has no intrinsic
                    sign: `isograph.utils.stable_sign` orients it by the largest loading, and
                    a different sample set can pick a different pivot. The QTL arm
                    recomputes `S_g` on the expanded cohort, so a swQTL slope can point the
                    opposite way to the discovery fit that defined the switch. This correlates
                    the recomputed `S_g` with the discovery `feature_scores` on the libraries
                    both contain, records `n_genes_sign_flipped`, and writes pinned slopes
                    BESIDE the permutation results.

  --stage control   Positive control for the `A_g` eQTL arm: if the genotype and phenotype
                    joins are right, BrainSEQ abundance-QTL must recover the tissue-matched
                    GTEx v11 brain eGenes, with concordant direction at GTEx's lead variant. A
                    broken donor join destroys replication; an allele-coding error sends
                    direction concordance towards 0 or 0.5. The pass rule is fixed below,
                    before the check was run.

Why pinning after mapping is exact rather than approximate: tensorQTL p-values do not
depend on the sign of a phenotype; the rank inverse-normal transform is odd
(average ranks give rank(-x) = n + 1 - rank(x), and the normal quantile is symmetric); and
the hidden factors are left singular vectors of the residual matrix, which flipping the sign
of individual phenotype columns leaves unchanged. A flipped gene therefore differs only in
the sign of its slopes, and that is all the pinned table corrects.

Outputs under 05_genetic_anchoring/_m/brainseq_switch_qtl/<arm>/checks/.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import cohort_dir, ensure_dir, stage_out
from isograph_benchmark.real_data.brainseq_switch_qtl import (
    _GENO_ROOT,
    ARMS,
    REGIONS,
    _storey_pi1,
    out_dir,
)
from isograph_benchmark.real_data.locus_event_audit import GTEX_V11

# Tissue-matched GTEx v11 brain tissue per BrainSEQ region (DLPFC ~ frontal cortex BA9).
GTEX_TISSUE = {
    "caudate": "Brain_Caudate_basal_ganglia",
    "dlpfc": "Brain_Frontal_Cortex_BA9",
    "hippocampus": "Brain_Hippocampus",
}

# --- Pre-specified 2026-09-11, before either check was run --------------------------------
FDR = 0.05
MIN_SHARED_LIBRARIES = 20      # below this a per-gene correlation is not estimable
AXIS_UNSTABLE_ABS_R = 0.5      # |r| below this: the recomputed PC1 is a different axis, not a flip
# The recomputed S_g must BE the discovery switch coordinate. Recomputed from identical counts
# with the discovery transcript filter it reproduces at median |r| 1.000 on the same libraries;
# the expanded QTL sample set shifts PC1 only slightly. Below this, the arm mapped a different
# phenotype and no sign pin can repair it (the 2026-09-11 unfiltered arm sat at 0.31-0.35).
SIGNPIN_REPRODUCE_MEDIAN_ABS_R = 0.99
CONTROL_PI1_MIN = 0.5          # pi1 of BrainSEQ A_g permutation p among GTEx tissue-matched eGenes
CONTROL_CONCORDANCE_MIN = 0.90  # direction concordance at GTEx lead variants, BrainSEQ p < 1e-5
CONTROL_NOMINAL_P = 1e-5
CONTROL_MIN_STRONG_PAIRS = 20
CONTROL_TOP_N = 1000           # the strongest GTEx eGenes, reported as a recovery rate

_BED_KEYS = ("#chr", "start", "end", "phenotype_id")
_FEATURE_KEYS = ("feature_id", "gene_id", "feature_type", "n_transcripts")
_AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def _bare(s: pd.Series) -> pd.Series:
    return s.astype(str).str.split(".", n=1).str[0]


def checks_dir(arm: str) -> Path:
    return ensure_dir(stage_out("anchoring.brainseq_qtl", arm) / "checks")


def arm_check_failures(arm: str, regions: tuple[str, ...] | list[str] | None = None,
                       dest: Path | None = None) -> list[str]:
    """Why an arm's QTL results may not be read yet; empty once both checks passed everywhere.

    Downstream analyses (colocalization, SMR) call this rather than assume the checks ran. A
    recomputed S_g that does not reproduce the discovery axis, or an A_g arm that fails the
    positive control, means the mapping itself is suspect, and nothing should be built on it.
    """
    regions = list(regions or REGIONS)
    dest = dest or (stage_out("anchoring.brainseq_qtl", arm) / "checks")
    out: list[str] = []

    sp = dest / "sign_pin_summary.parquet"
    if not sp.exists():
        out.append(f"no sign-pin summary at {sp}")
    else:
        s = pd.read_parquet(sp)
        s = s[s["arm"] == arm] if "arm" in s.columns else s
        for r in regions:
            row = s[s["region"] == r]
            if row.empty:
                out.append(f"{r}: sign pin not run")
            elif not bool(row["reproduces_discovery"].astype(bool).all()):
                out.append(f"{r}: recomputed S_g does not reproduce the discovery axis "
                           f"(median |r| {float(row['median_abs_r'].iloc[0]):.3f})")

    pc = dest / "positive_control_summary.parquet"
    if not pc.exists():
        out.append(f"no positive-control summary at {pc}")
    else:
        c = pd.read_parquet(pc)
        c = c[c["arm"] == arm] if "arm" in c.columns else c
        for r in regions:
            row = c[c["region"] == r]
            if row.empty:
                out.append(f"{r}: positive control not run")
            elif (row["verdict"] != "PASS").any():
                out.append(f"{r}: positive control {row['verdict'].iloc[0]} "
                           f"({row['fail_reasons'].iloc[0]})")
    return out


# --------------------------------------------------------------------------- #
# Sign pin
# --------------------------------------------------------------------------- #
def load_discovery_switch(region: str) -> pd.DataFrame:
    """Discovery-fit S_g: index = bare gene id, columns = library (RNum)."""
    f = cohort_dir("brainseq") / region / "_m" / "isograph_vae" / "feature_scores.parquet"
    if not f.exists():
        raise SystemExit(f"missing discovery feature scores: {f}")
    d = pd.read_parquet(f)
    d = d[d["feature_type"] == "switch"]
    cols = [c for c in d.columns if c not in _FEATURE_KEYS]
    out = d[cols].set_axis(_bare(d["gene_id"]).to_numpy(), axis=0)
    return out[~out.index.duplicated()]


def load_recomputed_switch(region: str, arm: str) -> pd.DataFrame:
    """QTL-arm S_g: index = bare gene id, columns = the library (RNum) the arm used per donor."""
    src = out_dir(arm, region)
    bed = pd.read_csv(src / "phenotypes_switch.bed.gz", sep="\t")
    samples = pd.read_parquet(src / "samples.parquet")[["RNum", "BrNum"]].astype(str)
    donors = [c for c in bed.columns if c not in _BED_KEYS]
    x = bed[donors].set_axis(_bare(bed["phenotype_id"]).to_numpy(), axis=0)
    x = x.rename(columns=dict(zip(samples["BrNum"], samples["RNum"])))
    return x[~x.index.duplicated()]


def sign_pin_table(discovery: pd.DataFrame, recomputed: pd.DataFrame,
                   min_shared: int = MIN_SHARED_LIBRARIES,
                   unstable_abs_r: float = AXIS_UNSTABLE_ABS_R) -> pd.DataFrame:
    """Per gene, Pearson r between discovery and recomputed S_g over the shared libraries.

    `sign` is the factor that orients the recomputed axis like the discovery one. A gene whose
    |r| is small is flagged `axis_unstable`: its PC1 is a different axis on the expanded cohort,
    so a sign cannot rescue a directional comparison and it should not carry one.
    """
    genes = discovery.index.intersection(recomputed.index)
    cols = discovery.columns.intersection(recomputed.columns)
    D = discovery.loc[genes, cols].to_numpy(dtype=float)
    R = recomputed.loc[genes, cols].to_numpy(dtype=float)
    rows = []
    for i, g in enumerate(genes):
        ok = np.isfinite(D[i]) & np.isfinite(R[i])
        n = int(ok.sum())
        r = np.nan
        if n >= min_shared and D[i, ok].std() > 0 and R[i, ok].std() > 0:
            r = float(np.corrcoef(D[i, ok], R[i, ok])[0, 1])
        rows.append((g, n, r))
    t = pd.DataFrame(rows, columns=["gene", "n_shared", "r"])
    t["sign"] = np.where(t["r"] < 0, -1.0, np.where(t["r"] >= 0, 1.0, np.nan))
    t["flipped"] = t["r"] < 0
    t["axis_unstable"] = t["r"].abs() < unstable_abs_r
    t["pinnable"] = t["r"].notna()
    return t


def pin_slopes(perm: pd.DataFrame, pin: pd.DataFrame) -> pd.DataFrame:
    """Permutation results with `slope_pinned` on the discovery orientation. A gene with no
    estimable pin keeps its raw slope and a missing `slope_pinned`, never a guessed sign."""
    p = perm.assign(gene=_bare(perm["phenotype_id"]))
    p = p.merge(pin[["gene", "n_shared", "r", "sign", "flipped", "axis_unstable"]],
                on="gene", how="left")
    p["slope_pinned"] = p["slope"] * p["sign"]
    return p.drop(columns="gene")


def run_signpin(arm: str = "all_samples", regions: list[str] | None = None,
                fdr: float = FDR) -> pd.DataFrame:
    dest = checks_dir(arm)
    rows = []
    for region in regions or list(REGIONS):
        src = out_dir(arm, region)
        if not (src / "phenotypes_switch.bed.gz").exists():
            print(f"  {region}: no recomputed S_g for {arm}; skipped")
            continue
        disc = load_discovery_switch(region)
        rec = load_recomputed_switch(region, arm)
        t = sign_pin_table(disc, rec)
        t.insert(0, "region", region)
        t.to_parquet(dest / f"sign_pin_{region}.parquet", index=False)
        ok = t[t["pinnable"]]
        row = {
            "region": region, "arm": arm,
            "n_discovery_libraries": int(disc.shape[1]),
            "n_shared_libraries": int(len(disc.columns.intersection(rec.columns))),
            "n_genes_discovery": int(disc.shape[0]),
            "n_genes_recomputed": int(rec.shape[0]),
            "n_genes_compared": int(len(ok)),
            "n_genes_not_comparable": int((~t["pinnable"]).sum()),
            "n_genes_sign_flipped": int(ok["flipped"].sum()),
            "frac_sign_flipped": float(ok["flipped"].mean()) if len(ok) else np.nan,
            "median_abs_r": float(ok["r"].abs().median()) if len(ok) else np.nan,
            "n_axis_unstable": int(ok["axis_unstable"].sum()),
            "reproduces_discovery": bool(
                len(ok) and ok["r"].abs().median() >= SIGNPIN_REPRODUCE_MEDIAN_ABS_R),
        }
        perm_f = src / "qtl" / "cis_qtl_switch.parquet"
        if perm_f.exists():
            pinned = pin_slopes(pd.read_parquet(perm_f), t)
            pinned.to_parquet(src / "qtl" / "cis_qtl_switch.sign_pinned.parquet", index=False)
            sig = pinned[pinned["qval"] < fdr]
            row.update({
                "n_swQTL": int(len(sig)),
                "n_swQTL_pinned": int(sig["sign"].notna().sum()),
                "n_swQTL_sign_flipped": int((sig["flipped"] == True).sum()),  # noqa: E712
                "n_swQTL_axis_unstable": int((sig["axis_unstable"] == True).sum()),  # noqa: E712
            })
        rows.append(row)
        print(f"  {region}: {row['n_genes_compared']:,} genes over "
              f"{row['n_shared_libraries']} shared libraries; flipped "
              f"{row['n_genes_sign_flipped']:,}; median |r| {row['median_abs_r']:.3f}; "
              f"axis-unstable {row['n_axis_unstable']:,}")
    s = pd.DataFrame(rows)
    if len(s):
        s.to_parquet(dest / "sign_pin_summary.parquet", index=False)
    return s


# --------------------------------------------------------------------------- #
# Positive control
# --------------------------------------------------------------------------- #
def load_gtex_egenes(region: str) -> pd.DataFrame:
    f = GTEX_V11 / "GTEx_Analysis_v11_eQTL" / f"{GTEX_TISSUE[region]}.v11.eGenes.txt.gz"
    cols = ["gene_id", "gene_name", "chr", "ref", "alt", "rs_id_dbSNP157_GRCh38p14", "af",
            "pval_nominal", "slope", "pval_beta", "qval"]
    g = pd.read_csv(f, sep="\t", usecols=cols).rename(
        columns={"rs_id_dbSNP157_GRCh38p14": "variant_id"})
    g["gene"] = _bare(g["gene_id"])
    return g


def gene_replication(gtex: pd.DataFrame, bs_perm: pd.DataFrame, fdr: float = FDR,
                     top_n: int = CONTROL_TOP_N) -> dict:
    """Do GTEx tissue-matched eGenes carry an eQTL signal in BrainSEQ A_g?"""
    b = bs_perm.assign(gene=_bare(bs_perm["phenotype_id"]))[["gene", "pval_beta", "qval"]]
    e = gtex[gtex["qval"] < fdr].merge(b, on="gene", how="inner", suffixes=("_gtex", "_bs"))
    top = e.nsmallest(top_n, "pval_nominal")
    return {
        "n_gtex_egenes": int((gtex["qval"] < fdr).sum()),
        "n_egenes_tested_in_brainseq": int(len(e)),
        "pi1_brainseq_given_gtex_egene": _storey_pi1(e["pval_beta_bs"].to_numpy()),
        "pi1_brainseq_all_genes": _storey_pi1(b["pval_beta"].to_numpy()),
        "frac_brainseq_q05_given_gtex_egene": (float((e["qval_bs"] < fdr).mean())
                                               if len(e) else np.nan),
        "n_top_gtex_egenes": int(len(top)),
        "frac_brainseq_q05_top_gtex_egenes": (float((top["qval_bs"] < fdr).mean())
                                              if len(top) else np.nan),
    }


def align_effects(pairs: pd.DataFrame) -> pd.DataFrame:
    """Put the BrainSEQ slope on GTEx's ALT allele and score direction concordance.

    `pairs` carries GTEx `ref`/`alt`/`slope_gtex` and BrainSEQ `bs_ref`/`bs_alt`/`slope_bs`
    (already on the BrainSEQ ALT allele). Same orientation is kept, including palindromic
    SNPs, because both panels code REF on the GRCh38 forward strand. A swapped orientation is
    flipped unless the SNP is palindromic, where a swap cannot be told from a strand error;
    any other allele mismatch (multi-allelic rsID, indel representation) is dropped.
    """
    same = (pairs["ref"] == pairs["bs_ref"]) & (pairs["alt"] == pairs["bs_alt"])
    swap = (pairs["ref"] == pairs["bs_alt"]) & (pairs["alt"] == pairs["bs_ref"])
    pal = pd.Series([(r, a) in _AMBIGUOUS for r, a in zip(pairs["ref"], pairs["alt"])],
                    index=pairs.index)
    keep = same | (swap & ~pal)
    out = pairs[keep].copy()
    out["slope_bs_aligned"] = np.where(same[keep], out["slope_bs"], -out["slope_bs"])
    out["concordant"] = np.sign(out["slope_bs_aligned"]) == np.sign(out["slope_gtex"])
    return out


def control_verdict(pi1: float, concordance_strong: float, n_strong: int) -> tuple[str, str]:
    reasons = []
    if not (np.isfinite(pi1) and pi1 >= CONTROL_PI1_MIN):
        reasons.append(f"pi1 {pi1:.3f} < {CONTROL_PI1_MIN}")
    if n_strong < CONTROL_MIN_STRONG_PAIRS:
        reasons.append(f"only {n_strong} lead-variant pairs at BrainSEQ p < {CONTROL_NOMINAL_P:g}")
    elif not (np.isfinite(concordance_strong) and concordance_strong >= CONTROL_CONCORDANCE_MIN):
        reasons.append(f"direction concordance {concordance_strong:.3f} < "
                       f"{CONTROL_CONCORDANCE_MIN}")
    return ("PASS" if not reasons else "FAIL"), "; ".join(reasons)


def read_pvar_alleles(arm: str, chrom: int, rsids: set[str]) -> pd.DataFrame:
    """REF/ALT and the panel ALT frequency for the requested rsIDs on one chromosome."""
    f = _GENO_ROOT / arm / f"chr{chrom}.pvar"
    rows = []
    with open(f) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5 or parts[2] not in rsids:
                continue
            af = np.nan
            if len(parts) > 5:
                for kv in parts[5].split(";"):
                    if kv.startswith("AF="):
                        af = float(kv[3:])
                        break
            rows.append((parts[2], parts[3], parts[4], af))
    return pd.DataFrame(rows, columns=["variant_id", "bs_ref", "bs_alt", "panel_alt_af"])


def nominal_at(arm: str, region: str, chrom: int, keys: pd.DataFrame) -> pd.DataFrame:
    """BrainSEQ A_g nominal statistics for the requested (gene, rsID) pairs on one chromosome."""
    import pyarrow.dataset as ds
    files = sorted((out_dir(arm, region) / "qtl").glob(
        f"abundance.chr{chrom}.cis_qtl_pairs.*.parquet"))
    cols = ["phenotype_id", "variant_id", "af", "pval_nominal", "slope", "slope_se"]
    if not files or keys.empty:
        return pd.DataFrame(columns=[*cols, "gene"])
    tbl = ds.dataset([str(f) for f in files], format="parquet").to_table(
        columns=cols, filter=ds.field("variant_id").isin(sorted(set(keys["variant_id"]))))
    d = tbl.to_pandas()
    d["gene"] = _bare(d["phenotype_id"])
    return d.merge(keys[["gene", "variant_id"]].drop_duplicates(), on=["gene", "variant_id"],
                   how="inner")


def run_control(arm: str = "all_samples", regions: list[str] | None = None,
                chroms: list[int] | None = None, fdr: float = FDR) -> pd.DataFrame:
    dest = checks_dir(arm)
    regions = [r for r in (regions or list(REGIONS)) if r in GTEX_TISSUE]
    chroms = chroms or list(range(1, 23))

    gtex = {}
    for region in regions:
        g = load_gtex_egenes(region)
        g = g[(g["qval"] < fdr) & g["variant_id"].astype(str).str.startswith("rs")]
        g = g.assign(chrom=pd.to_numeric(g["chr"].str.replace("chr", "", regex=False),
                                         errors="coerce"))
        gtex[region] = g[g["chrom"].isin(chroms)]
    pvar = {}
    for ch in chroms:
        rs = set().union(*(set(g.loc[g["chrom"] == ch, "variant_id"]) for g in gtex.values()))
        if rs:
            pvar[ch] = read_pvar_alleles(arm, ch, rs)

    rows = []
    for region in regions:
        perm_f = out_dir(arm, region) / "qtl" / "cis_qtl_abundance.parquet"
        if not perm_f.exists():
            print(f"  {region}: no A_g permutation results for {arm}; skipped")
            continue
        perm = pd.read_parquet(perm_f)
        if chroms != list(range(1, 23)):
            perm = perm[perm["chrom"].isin(chroms)]
        full = load_gtex_egenes(region)
        full = full[pd.to_numeric(full["chr"].str.replace("chr", "", regex=False),
                                  errors="coerce").isin(chroms)]
        rep = gene_replication(full, perm, fdr)

        g = gtex[region]
        parts = []
        for ch in chroms:
            gk = g[g["chrom"] == ch]
            if gk.empty or ch not in pvar:
                continue
            nom = nominal_at(arm, region, ch, gk)
            if nom.empty:
                continue
            parts.append(nom.merge(gk[["gene", "variant_id", "gene_name", "ref", "alt", "af",
                                       "slope", "pval_nominal"]].rename(
                columns={"af": "af_gtex", "slope": "slope_gtex",
                         "pval_nominal": "pval_nominal_gtex"}),
                on=["gene", "variant_id"]).merge(pvar[ch], on="variant_id", how="inner"))
        pairs = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

        # Which allele does tensorQTL count? Read it from the data rather than assume it: the
        # BrainSEQ allele frequency tracks the panel's ALT frequency if the ALT is counted.
        corr_af = (float(np.corrcoef(pairs["af"], pairs["panel_alt_af"])[0, 1])
                   if len(pairs) > 2 and pairs["panel_alt_af"].notna().sum() > 2 else np.nan)
        counted_alt = not (np.isfinite(corr_af) and corr_af < 0)
        if len(pairs):
            pairs["slope_bs"] = pairs["slope"] * (1.0 if counted_alt else -1.0)
            aligned = align_effects(pairs)
        else:
            aligned = pd.DataFrame(columns=["concordant", "pval_nominal"])
        strong = aligned[aligned["pval_nominal"] < CONTROL_NOMINAL_P]
        conc_all = float(aligned["concordant"].mean()) if len(aligned) else np.nan
        conc_strong = float(strong["concordant"].mean()) if len(strong) else np.nan
        verdict, why = control_verdict(rep["pi1_brainseq_given_gtex_egene"], conc_strong,
                                       int(len(strong)))
        aligned.assign(region=region).to_parquet(dest / f"positive_control_pairs_{region}.parquet",
                                                 index=False)
        row = {"region": region, "arm": arm, "gtex_tissue": GTEX_TISSUE[region],
               "chroms": "all" if chroms == list(range(1, 23)) else ",".join(map(str, chroms)),
               **rep,
               "n_lead_pairs_found": int(len(pairs)),
               "n_lead_pairs_aligned": int(len(aligned)),
               "corr_af_brainseq_vs_panel_alt": corr_af,
               "counted_allele": "ALT" if counted_alt else "REF",
               "concordance_all_aligned": conc_all,
               "n_strong_pairs": int(len(strong)),
               "concordance_strong": conc_strong,
               "verdict": verdict, "fail_reasons": why}
        rows.append(row)
        print(f"  {region}: pi1 {row['pi1_brainseq_given_gtex_egene']:.3f} over "
              f"{row['n_egenes_tested_in_brainseq']:,} GTEx eGenes; top-{CONTROL_TOP_N} "
              f"recovery {row['frac_brainseq_q05_top_gtex_egenes']:.3f}; concordance "
              f"{conc_strong:.3f} over {len(strong):,} strong pairs "
              f"(counted allele {row['counted_allele']}, af corr {corr_af:.2f}) -> {verdict}")
    s = pd.DataFrame(rows)
    if len(s):
        s.to_parquet(dest / "positive_control_summary.parquet", index=False)
    return s


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def _f(v, nd=3):
    return "—" if v is None or (isinstance(v, float) and not np.isfinite(v)) else (
        f"{v:.{nd}f}" if isinstance(v, float) else f"{v:,}" if isinstance(v, (int, np.integer))
        else str(v))


def write_report(arm: str) -> Path:
    dest = checks_dir(arm)
    L: list[str] = []
    A = L.append
    A(f"# BrainSEQ QTL checks — `{arm}`")
    A("")
    A("Two checks required before any BrainSEQ QTL result is read. Neither modifies the mapping "
      "outputs. Thresholds were fixed in `brainseq_qtl_checks.py` before the checks were run.")
    A("")
    sp = dest / "sign_pin_summary.parquet"
    if sp.exists():
        s = pd.read_parquet(sp)
        A("## Sign pin of the recomputed switch coordinate")
        A("")
        A("`S_g` is PC1 of a within-gene CLR composition, oriented by `stable_sign` on the "
          "largest loading; recomputed on the expanded QTL cohort it can land on the opposite "
          "sign from the discovery fit. Per gene, the recomputed and discovery `S_g` are "
          f"correlated over the libraries both contain (at least {MIN_SHARED_LIBRARIES}). A "
          f"negative r is a flip; |r| < {AXIS_UNSTABLE_ABS_R} marks a gene whose PC1 is a "
          "different axis on the expanded cohort, for which no sign makes a directional "
          "comparison meaningful. Pinned slopes are in `qtl/cis_qtl_switch.sign_pinned.parquet`; "
          "p-values, q-values and hidden factors are unaffected by construction.")
        A("")
        A(f"**Pre-specified gate:** the recomputed `S_g` must reproduce the discovery coordinate "
          f"at median |r| >= {SIGNPIN_REPRODUCE_MEDIAN_ABS_R}. Below that the arm mapped a "
          "different phenotype, and no sign pin can repair it: re-map before reading any swQTL "
          "result. (The first all_samples arm omitted the discovery transcript filter and sat at "
          "0.31-0.35; recomputed with it, the same libraries reproduce at 1.000.)")
        A("")
        A("| region | shared libraries | genes compared | median abs r | reproduces discovery | "
          "sign flipped | axis-unstable | swQTL genes | swQTL flipped | swQTL axis-unstable |")
        A("|---|---|---|---|---|---|---|---|---|---|")
        for r in s.itertuples(index=False):
            rep = getattr(r, "reproduces_discovery", None)
            A(f"| {r.region} | {r.n_shared_libraries} | {_f(r.n_genes_compared)} | "
              f"{_f(r.median_abs_r)} | **{'yes' if rep else 'NO'}** | "
              f"{_f(r.n_genes_sign_flipped)} ({_f(r.frac_sign_flipped)}) | "
              f"{_f(r.n_axis_unstable)} | "
              f"{_f(getattr(r, 'n_swQTL', None))} | {_f(getattr(r, 'n_swQTL_sign_flipped', None))} | "
              f"{_f(getattr(r, 'n_swQTL_axis_unstable', None))} |")
        A("")
    pc = dest / "positive_control_summary.parquet"
    if pc.exists():
        s = pd.read_parquet(pc)
        A("## Positive control: the `A_g` eQTL arm against GTEx v11 brain eGenes")
        A("")
        A("Tissue-matched GTEx v11 eGenes (BH q < 0.05): caudate ↔ Caudate basal ganglia, DLPFC ↔ "
          "Frontal Cortex BA9, hippocampus ↔ Hippocampus. **Pre-specified pass rule**, per region: "
          f"(i) Storey pi1 of the BrainSEQ `A_g` permutation p among GTEx eGenes >= "
          f"{CONTROL_PI1_MIN}, and (ii) direction concordance at GTEx's lead variant >= "
          f"{CONTROL_CONCORDANCE_MIN} among aligned pairs with BrainSEQ nominal p < "
          f"{CONTROL_NOMINAL_P:g} (at least {CONTROL_MIN_STRONG_PAIRS} pairs). A broken donor "
          "join destroys (i); an allele-coding error drives (ii) towards 0 or 0.5.")
        A("")
        A("The counted allele is read from the data (BrainSEQ allele frequency against the "
          "panel's ALT frequency), not assumed. Palindromic SNPs are kept only when both panels "
          "agree on REF, since both code REF on the GRCh38 forward strand.")
        A("")
        A("| region | GTEx eGenes tested | pi1 | q<0.05 | top-1000 recovery | pi1 all genes | "
          "lead pairs aligned | counted allele (af r) | concordance all | strong pairs | "
          "concordance strong | verdict |")
        A("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in s.itertuples(index=False):
            A(f"| {r.region} | {_f(r.n_egenes_tested_in_brainseq)} | "
              f"{_f(r.pi1_brainseq_given_gtex_egene)} | {_f(r.frac_brainseq_q05_given_gtex_egene)} | "
              f"{_f(r.frac_brainseq_q05_top_gtex_egenes)} | {_f(r.pi1_brainseq_all_genes)} | "
              f"{_f(r.n_lead_pairs_aligned)} | {r.counted_allele} "
              f"({_f(r.corr_af_brainseq_vs_panel_alt, 2)}) | {_f(r.concordance_all_aligned)} | "
              f"{_f(r.n_strong_pairs)} | {_f(r.concordance_strong)} | **{r.verdict}** |")
        fails = s[s["verdict"] != "PASS"]
        if len(fails):
            A("")
            A("**Failures:** " + "; ".join(f"{r.region}: {r.fail_reasons}"
                                          for r in fails.itertuples(index=False)))
        if (s["chroms"] != "all").any():
            A("")
            A("*Smoke run on a chromosome subset; not the pre-specified check.*")
        A("")
    f = dest / "QTL_CHECKS.md"
    f.write_text("\n".join(L) + "\n")
    return f


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("signpin", "control", "all"), required=True)
    ap.add_argument("--arm", choices=ARMS, default="all_samples")
    ap.add_argument("--region", choices=REGIONS, action="append", default=None)
    ap.add_argument("--chrom", type=int, action="append", default=None,
                    help="control: restrict to these chromosomes (smoke tests only)")
    args = ap.parse_args(argv)
    if args.stage in ("signpin", "all"):
        print(f"== sign pin / {args.arm} ==")
        run_signpin(args.arm, args.region)
    if args.stage in ("control", "all"):
        print(f"== positive control / {args.arm} ==")
        run_control(args.arm, args.region, args.chrom)
    print(f"  report -> {write_report(args.arm)}")


if __name__ == "__main__":
    main()
