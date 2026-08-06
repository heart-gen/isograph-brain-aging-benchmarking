"""Orthogonal switch-detection concordance: IsoGraph co-switch modules vs satuRn DTU.

IsoGraph calls a gene "switching" from the correlation structure of its CLR switch
coordinate (co-switching modules). This CLI asks whether an INDEPENDENT, standard
differential-transcript-usage (DTU) caller flags the same genes as switching on the
same quantifications -- i.e. whether IsoGraph's module switches are real isoform
switches a per-gene tool would also detect.

The independent caller is **satuRn** (Gilis et al. 2022), the DTU test that
**IsoformSwitchAnalyzeR v2** runs internally (``isoformSwitchTestSatuRn``). satuRn is
already installed in the rnaseq R env; the full IsoformSwitchAnalyzeR Bioconductor
stack is not, and its extra machinery (ORF/Pfam/NMD consequence annotation) is not
needed here -- IsoGraph computes switch consequences separately
(``switch_consequence.py``). So this is the ISA switch *test*, minus the consequence
wrapper, and is honestly reported as such.

Design (mirrors ``validate_switch_splicing.gene_corroboration`` Analysis A and the
lab's age-by-ancestry satuRn pipeline, ancestry-aging-adrd-risk/
differential_expression):
  * satuRn's empirical-null FDR (its defining feature) only applies to a single-df
    contrast, so the exposure enters as one design coefficient. ``--trait age``
    models age *continuously* (linear z-scored ``age_z``, every sample -- no tertile
    dichotomy) and tests that coefficient; ``--trait dx`` compares SCZD vs Control in
    the caudate SCZD cohort (``group`` releveled to Control). NB a multi-df spline
    omnibus would bypass the empirical-null recalibration and, on satuRn's UNSCALED
    ``vcovUnsc`` (true cov = dispersion*vcovUnsc, dispersion ~ O(10)), is badly
    anti-conservative -- so it is deliberately avoided here.
  * The design is covariate-adjusted with a lean subset of the paper's covariates
    (Sex/RIN/mapping/mito [+Age for dx]) matched to IsoGraph's discovery design, so
    the two methods share an exposure model; SNP PCs / MoD / SVA are omitted (satuRn
    is an INDEPENDENT detector and design-matching matters more here than maximal
    deconfounding).
  * Per-gene satuRn evidence = ``-log10(min isoform p)`` -- a *continuous* score.
    ``switch_pos`` is membership in a trait-associated co-switch module (the
    manuscript unit; reused from ``validate_switch_splicing._module_switch_genes``).
  * PRIMARY test is threshold-free: does switch-module membership predict stronger
    independent satuRn DTU evidence? Mann-Whitney (one-sided) + rank-biserial effect
    size on that evidence, then a transcript-count-adjusted logistic OR (per unit
    -log10 p). We deliberately do NOT binarize on satuRn's empirical-null FDR: that
    null follows Efron's locfdr, which assumes a *sparse* alternative, whereas brain
    age-DTU is pervasive -- so at full genome scale the empirical null over-widens and
    collapses discoveries to 0-2 genes, leaving a binary Fisher test unpowered despite
    OR=7-12 in the right direction. The empirical-FDR positive count is still reported
    as a conservative secondary, but is never the headline.

Outputs land in ``real_data/_m/isa_concordance/<cohort>_<region>_<trait>/``:
  * ``gene_concordance.parquet`` -- per-gene switch_pos / evidence / isa_emp_q / n_tx.
  * ``saturn_transcript_results.tsv.gz`` -- raw satuRn transcript-level output.
  * ``summary.json`` -- MWU p + rank-biserial, adjusted OR/p, median evidence,
    conservative empirical-FDR count, n's, provenance.
"""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from isograph.io.artifacts import load_dataset_bundle

from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.run_models import _filter_expressed_transcripts
from isograph_benchmark.real_data.validate_switch_splicing import (
    BRAINSEQ_AGING_REGIONS,
    GTEX_REGIONS,
    _artifact_dir,
    _bundle_path,
    _load_switch_scores,
    _module_switch_genes,
    _strip_ver,
)

# satuRn lives in the rnaseq R env (the manuscript figure/rnaseq interpreter).
_RSCRIPT = "/ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript"
_R_RUNNER = Path(__file__).with_name("isa_concordance_saturn.R")

# Lean covariate subset for the satuRn design, matched to IsoGraph's discovery design
# (no SNP PCs / MoD / SVA): satuRn is an INDEPENDENT detector, so sharing an exposure
# model with IsoGraph matters more than maximal deconfounding.
_SATURN_COVARIATES = {
    ("brainseq", "age"): ["Sex", "RIN", "mapping_rate", "mito_rate"],
    ("brainseq", "dx"): ["Age", "Sex", "RIN", "mapping_rate", "mito_rate"],
    ("gtex", "age"): ["SEX", "SMRIN", "SMTSISCH", "SMMAPRT"],
}
_AGE_COL = {"brainseq": "Age", "gtex": "AGE"}
_DX_REF = "Control"  # reference level; coefficient groupSCZD = SCZD - Control


def _exposure_spec(trait: str) -> dict:
    """How the exposure enters satuRn's single-coefficient testDTU contrast.

    Both traits use one design coefficient so satuRn's empirical-null FDR applies
    (a multi-df spline omnibus would bypass that calibration). age -> continuous
    z-scored age (linear ``age_z`` coefficient, all samples); dx -> factor group
    releveled to Control (``groupSCZD`` coefficient).
    """
    if trait == "age":
        return {"col": "age_z", "kind": "numeric", "ref": None, "test_coef": "age_z"}
    hi, _ = _contrast_levels(trait)  # ("SCZD", "Control")
    return {"col": "group", "kind": "factor", "ref": _DX_REF, "test_coef": f"group{hi}"}


def _out_dir(cohort: str, region: str, trait: str):
    return ensure_dir(rel("real_data", "_m", "isa_concordance", f"{cohort}_{region}_{trait}"))


# --------------------------------------------------------------------------- #
# Exposure construction
# --------------------------------------------------------------------------- #
def _build_exposure(sample_table: pd.DataFrame, cohort: str, trait: str) -> pd.Series:
    """Return a Series (index = sample_id) of the modelled exposure, NaN where unusable.

    age -> continuous z-scored age (``age_z``), spanning all samples.
    dx  -> 'SCZD' vs 'Control' (NaN for any other Dx).
    """
    st = sample_table.copy()
    st["sample_id"] = st["sample_id"].astype(str)
    st = st.set_index("sample_id")
    if trait == "dx":
        dx = st["Dx"].astype(str)
        return pd.Series(np.where(dx.isin(["SCZD", "Control"]), dx, np.nan), index=st.index)
    age = pd.to_numeric(st[_AGE_COL[cohort]], errors="coerce")
    age_z = (age - age.mean()) / age.std(ddof=0)
    return age_z.astype(float)


def _contrast_levels(trait: str) -> tuple[str, str]:
    # (case/high level, reference/low level) -> satuRn contrast = case - reference
    return ("SCZD", "Control")


# --------------------------------------------------------------------------- #
# satuRn input staging + invocation
# --------------------------------------------------------------------------- #
def _stage_and_run_saturn(
    cohort: str, region: str, trait: str, out, sample_table: pd.DataFrame,
    tc: np.ndarray, tt: pd.DataFrame, exposure: pd.Series, cores: int, smoke: int,
) -> pd.DataFrame:
    """Write counts/colData/txInfo, invoke satuRn, return transcript-level results."""
    spec = _exposure_spec(trait)
    covs = _SATURN_COVARIATES[(cohort, trait)]
    st = sample_table.copy()
    st["sample_id"] = st["sample_id"].astype(str)
    st = st.set_index("sample_id")

    # Samples: those with a usable exposure AND complete design covariates.
    keep = exposure.dropna().index
    keep = st.loc[keep, covs].dropna().index
    keep_pos = [i for i, s in enumerate(st.index) if s in set(keep)]
    sample_ids = [st.index[i] for i in keep_pos]
    counts = tc[:, keep_pos]

    # Restrict to multi-isoform genes with >=1 isoform still expressed in this subset.
    expressed = (counts > 0).sum(axis=1) >= max(3, int(0.1 * counts.shape[1]))
    txi = tt.loc[expressed].reset_index(drop=True)
    counts = np.rint(counts[expressed]).astype(int)
    gcount = txi["gene_id"].value_counts()
    multi = set(gcount[gcount >= 2].index)
    keep_tx = txi["gene_id"].isin(multi).to_numpy()
    txi = txi.loc[keep_tx].reset_index(drop=True)
    counts = counts[keep_tx]

    if smoke:
        genes = list(dict.fromkeys(txi["gene_id"]))[:smoke]
        m = txi["gene_id"].isin(set(genes)).to_numpy()
        txi, counts = txi.loc[m].reset_index(drop=True), counts[m]

    if spec["kind"] == "numeric":
        span = f"{len(sample_ids)} samples (continuous {spec['col']}, linear)"
    else:
        span = (f"{len(sample_ids)} samples ("
                + ", ".join(f"{g}={int((exposure.loc[sample_ids] == g).sum())}"
                            for g in _contrast_levels(trait)) + ")")
    print(f"[{cohort}/{region}/{trait}] satuRn input: {counts.shape[0]} isoforms / "
          f"{txi['gene_id'].nunique()} multi-isoform genes x {span}", flush=True)

    stage = ensure_dir(out / "_saturn_inputs")
    counts_df = pd.DataFrame(counts, columns=sample_ids)
    counts_df.insert(0, "isoform_id", txi["transcript_id"].to_numpy())
    counts_path = stage / "counts.tsv.gz"
    counts_df.to_csv(counts_path, sep="\t", index=False)

    coldata = st.loc[sample_ids, covs].copy()
    coldata.insert(0, spec["col"], exposure.loc[sample_ids].to_numpy())
    coldata.insert(0, "sample_id", sample_ids)
    coldata_path = stage / "coldata.tsv"
    coldata.to_csv(coldata_path, sep="\t", index=False)

    txinfo = pd.DataFrame({"isoform_id": txi["transcript_id"].to_numpy(),
                           "gene_id": txi["gene_id"].to_numpy()})
    txinfo_path = stage / "txinfo.tsv"
    txinfo.to_csv(txinfo_path, sep="\t", index=False)

    res_path = out / "saturn_transcript_results.tsv.gz"
    cmd = [
        _RSCRIPT, str(_R_RUNNER),
        "--counts", str(counts_path), "--coldata", str(coldata_path),
        "--txinfo", str(txinfo_path), "--covariates", ",".join(covs),
        "--exposure-col", spec["col"], "--exposure-kind", spec["kind"],
        "--test-coef", spec["test_coef"], "--out", str(res_path),
        "--seed", "13", "--cores", str(cores),
    ]
    if spec["kind"] == "factor":
        cmd += ["--ref-level", spec["ref"]]
    print("[isa] running satuRn:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

    with gzip.open(res_path, "rt") as fh:
        return pd.read_csv(fh, sep="\t")


# --------------------------------------------------------------------------- #
# Concordance (mirrors validate_switch_splicing Analysis A)
# --------------------------------------------------------------------------- #
_PVAL_FLOOR = 1e-300  # cap for -log10 so a machine-zero p does not become +inf


def _gene_isa_evidence(res: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Per-gene satuRn DTU evidence for a threshold-free concordance.

    Primary readout is a *continuous* gene-level evidence score, evidence =
    -log10(min raw satuRn p over the gene's isoforms). We deliberately do NOT
    binarize on satuRn's empirical-null FDR: that null follows Efron's locfdr,
    which assumes a *sparse* alternative, whereas brain age-DTU is pervasive
    (~27% of transcripts raw-p<0.05), so the empirical null is estimated too
    wide and collapses genome-wide discoveries to 0-2 genes -- leaving the
    enrichment test unpowered despite OR=7-12 in the right direction. The
    continuous score uses the full ranked signal instead of a fragile cut.

    ``isa_pos`` (min empirical FDR <= alpha) is still computed, but only as a
    conservative secondary count reported in the summary, never as the primary
    test.
    """
    r = res.copy()
    r["gene"] = _strip_ver(r["gene_id"])
    r["pval"] = pd.to_numeric(r["pval"], errors="coerce")
    emp = pd.to_numeric(
        r["empirical_FDR"].where(r["empirical_FDR"].notna(), r["regular_FDR"]),
        errors="coerce",
    )
    r["emp_q"] = emp
    g = r.dropna(subset=["pval"]).groupby("gene")
    per_gene = g["pval"].min().rename("isa_min_p").reset_index()
    per_gene["evidence"] = -np.log10(per_gene["isa_min_p"].clip(lower=_PVAL_FLOOR))
    emp_min = r.dropna(subset=["emp_q"]).groupby("gene")["emp_q"].min()
    per_gene["isa_emp_q"] = per_gene["gene"].map(emp_min)
    per_gene["isa_pos"] = per_gene["isa_emp_q"] <= alpha  # conservative secondary only
    return per_gene


def concordance(
    cohort: str, region: str, variant: str, trait: str, alpha: float,
    sample_table: pd.DataFrame, tc: np.ndarray, tt: pd.DataFrame, cores: int, smoke: int,
    out,
) -> tuple[pd.DataFrame, dict]:
    exposure = _build_exposure(sample_table, cohort, trait)
    res = _stage_and_run_saturn(cohort, region, trait, out, sample_table, tc, tt, exposure, cores, smoke)
    isa = _gene_isa_evidence(res, alpha)

    # IsoGraph switch-scored gene universe (+ n_transcripts for the adjusted logistic).
    sample_ids = list(sample_table["sample_id"].astype(str))
    sw_meta, _, _ = _load_switch_scores(cohort, region, variant, trait, sample_ids)
    sw = sw_meta[["gene", "n_transcripts"]].drop_duplicates("gene").copy()
    mod_genes = _module_switch_genes(cohort, region, variant, trait, alpha)
    sw["switch_pos"] = sw["gene"].isin(mod_genes)

    genes = sw.merge(isa, on="gene", how="inner")  # genes IsoGraph-scored AND satuRn-tested

    # ---- PRIMARY: threshold-free continuous concordance --------------------- #
    # Does co-switch-module membership predict stronger *independent* satuRn DTU
    # evidence? Mann-Whitney (one-sided) + rank-biserial effect size on the
    # gene-level evidence score; then a transcript-count-adjusted logistic OR.
    sw_ev = genes.loc[genes["switch_pos"], "evidence"].to_numpy()
    bg_ev = genes.loc[~genes["switch_pos"], "evidence"].to_numpy()
    n_sw, n_bg = int(sw_ev.size), int(bg_ev.size)
    mwu_u = mwu_p = rank_biserial = float("nan")
    if n_sw > 0 and n_bg > 0:
        mwu_u, mwu_p = stats.mannwhitneyu(sw_ev, bg_ev, alternative="greater")
        rank_biserial = float(2.0 * mwu_u / (n_sw * n_bg) - 1.0)  # +1 => switch >> bg
    median_ev_switch = float(np.median(sw_ev)) if n_sw else float("nan")
    median_ev_background = float(np.median(bg_ev)) if n_bg else float("nan")

    adj_or = adj_p = float("nan")
    try:
        if n_sw == 0 or n_bg == 0:
            raise ValueError("no switch (or no background) genes; concordance is vacuous")
        import statsmodels.api as sm
        m = genes.dropna(subset=["n_transcripts", "evidence"]).copy()
        X = pd.DataFrame({
            "const": 1.0,
            "evidence": m["evidence"].astype(float),
            "log_n_tx": np.log1p(m["n_transcripts"].astype(float)),
        })
        fit = sm.Logit(m["switch_pos"].astype(float), X).fit(disp=0)
        adj_or = float(np.exp(fit.params["evidence"]))  # OR of switch membership / unit -log10(p)
        adj_p = float(fit.pvalues["evidence"])
    except Exception as exc:  # pragma: no cover - defensive
        print(f"[{region}] adjusted logistic skipped: {exc}", flush=True)

    # ---- SECONDARY (conservative): empirical-FDR binary count, reported only - #
    n_isa_emp_pos = int(genes["isa_pos"].sum())
    n_switch_isa_emp = int((genes["switch_pos"] & genes["isa_pos"]).sum())

    summary = {
        "cohort": cohort, "region": region, "variant": variant, "trait": trait,
        "caller": "satuRn (DTU test used by IsoformSwitchAnalyzeR v2)",
        "test": (f"testDTU coefficient '{_exposure_spec(trait)['test_coef']}'"
                 + (" (continuous age, linear)" if trait == "age"
                    else f" ({'-'.join(_contrast_levels(trait))})")),
        "primary_test": ("Mann-Whitney (one-sided) + rank-biserial on gene-level "
                         "satuRn evidence = -log10(min isoform p); "
                         "transcript-count-adjusted logistic OR"),
        "n_genes": int(len(genes)),
        "n_switch_genes": n_sw,
        "n_background_genes": n_bg,
        "median_evidence_switch": median_ev_switch,
        "median_evidence_background": median_ev_background,
        "mwu_u": float(mwu_u), "mwu_p": float(mwu_p),
        "rank_biserial": rank_biserial,
        "adjusted_or": adj_or, "adjusted_p": adj_p,
        # conservative secondary (empirical-null FDR is over-conservative for
        # pervasive age-DTU; reported for transparency, not the headline):
        "empirical_fdr_note": ("empirical-null FDR assumes a sparse alternative "
                               "(Efron locfdr); pervasive brain DTU violates it, "
                               "so this binary count is a conservative lower bound"),
        "n_isa_empirical_pos": n_isa_emp_pos,
        "n_switch_and_isa_empirical_pos": n_switch_isa_emp,
        "alpha": alpha, "smoke": int(smoke),
    }
    return genes, summary


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def run_region(cohort: str, region: str, variant: str, trait: str, alpha: float,
               cores: int, smoke: int) -> dict:
    bundle = load_dataset_bundle(_bundle_path(cohort, region, trait))
    sample_table = bundle.sample_table
    tc, tt = _filter_expressed_transcripts(
        bundle.matrices["transcript_counts"], bundle.feature_tables["transcript"])
    del bundle
    out = _out_dir(cohort, region, trait)
    genes, summary = concordance(cohort, region, variant, trait, alpha,
                                 sample_table, tc, tt, cores, smoke, out)
    genes.to_parquet(out / "gene_concordance.parquet", index=False, compression="zstd")
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    s = summary
    print(f"[{cohort}/{region}/{trait}] ISA(satuRn) evidence: switch median "
          f"{s['median_evidence_switch']:.2f} vs background "
          f"{s['median_evidence_background']:.2f} (-log10 p) | rank-biserial "
          f"{s['rank_biserial']:.3f} MWU p={s['mwu_p']:.2e} | adj OR/unit="
          f"{s['adjusted_or']:.2f} p={s['adjusted_p']:.2e} | empFDR-pos "
          f"{s['n_isa_empirical_pos']} -> {out}", flush=True)
    return summary


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cohort", choices=["brainseq", "gtex"], default="brainseq")
    p.add_argument("--region", action="append", dest="regions",
                   help="Region(s); repeatable. Default: BrainSEQ aging regions / caudate (dx) / all GTEx.")
    p.add_argument("--trait", choices=["age", "dx"], default="age",
                   help="age -> continuous linear age_z testDTU coef; dx -> SCZD vs Control (caudate, brainseq only).")
    p.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--cores", type=int, default=1, help="satuRn BiocParallel workers.")
    p.add_argument("--smoke", type=int, default=0,
                   help="If >0, subset to this many genes for a fast login-node validation.")
    args = p.parse_args()

    if args.cohort == "gtex":
        if args.trait == "dx":
            raise SystemExit("--trait dx is BrainSEQ-only (GTEx has no case/control cohort).")
        regions = args.regions or list(GTEX_REGIONS)
    elif args.trait == "dx":
        regions = args.regions or ["caudate"]
        if regions != ["caudate"]:
            raise SystemExit("--trait dx is caudate-only (the SCZD cohort is caudate).")
    else:
        regions = args.regions or list(BRAINSEQ_AGING_REGIONS)

    for region in regions:
        run_region(args.cohort, region, args.variant, args.trait, args.alpha,
                   args.cores, args.smoke)


if __name__ == "__main__":
    main()
