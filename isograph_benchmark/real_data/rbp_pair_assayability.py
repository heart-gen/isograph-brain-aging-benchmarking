"""Score the prespecified IsoGraph switch pairs of an RBP panel for wet-lab assayability.

`rbp_target_panel.py` freezes *which* transcript pairs a perturbation experiment should read
out (`rbp_target_switch_pairs.parquet`: the motif-differential pairs of the panel genes). It
does not say whether those pairs can actually be measured at the bench, and a pair that cannot
move is not a test of anything. This CLI answers the two bench-facing questions separately,
because they fail for different reasons and have different fixes:

  MEASURABLE  -- can an isoform-ratio assay resolve the two transcripts, and are both of them
                 actually present in brain at a ratio that is not pinned to 0 or 1?

                 * structural: each member needs a feature the other lacks. A splice junction
                   unique to one transcript is the cheap readout (a primer sitting across it
                   amplifies that isoform only), so `assay_class == "junction"` is the design
                   we want. Failing that, a long enough unique exonic stretch supports a
                   probe/primer inside it (`unique_segment`). Pairs whose shorter member is
                   nested inside the longer one at a terminus (`nested_terminal`) have nothing
                   unique to amplify and need a 5'/3'-specific design, so they do not count.
                 * expression: both transcripts must clear `--min-tpm` on the median GTEx
                   brain sample, and the carrier fraction must sit inside
                   [`--min-frac`, 1 - `--min-frac`]. A pair whose minor isoform is 2% of the
                   gene has no headroom for a knockdown to shift.

  RESPONSIVE  -- does that ratio demonstrably move in human brain? A pair can be perfectly
                 measurable and still be a constitutive ratio that nothing perturbs. We use
                 the aging axis IsoGraph itself was built on: logit(carrier fraction) is
                 regressed on age with the standard GTEx technical covariates, and the pair is
                 responsive if the age term survives BH within the RBP *and* the fraction's
                 interquartile range clears `--min-iqr`. The IQR gate is the bench-relevant
                 half -- a significant but 1-point shift is not a plate-readable effect.

The two gates are reported independently rather than collapsed into one score, so a pair that
fails can be triaged: `structurally_measurable == False` is an assay-design problem (switch to
a 3'-end assay, or drop the pair), while `responsive == False` on a measurable pair means the
ratio is real but static and the pair is a poor primary endpoint.

Age is the only perturbation we can observe in human tissue; it is a proxy for "this ratio is
regulatable", not evidence that the RBP is what regulates it. That is the experiment's job.

Inputs
  08_integration/_m/rbp_target_panel/<RBP>/rbp_target_switch_pairs.parquet   (the frozen pairs)
  inputs/processed/gtex_v11/<region>/transcript_tpm.parquet             (RSEM TPM)
  inputs/bundles/gtex_v11_brain/<region>/samples.parquet                (AGE + covariates)
  GENCODE GTF (cached parquet)                                          (exon structure)

Outputs (08_integration/_m/rbp_target_panel/<RBP>/)
  rbp_pair_assayability.parquet/tsv         one row per pair, summarised across regions
  rbp_pair_assayability_by_region.parquet/tsv   one row per pair x region (the fitted detail)

Only GTEx contributes the expression arm: it is one uniform RSEM TPM quantification across 13
brain regions, whereas the BrainSEQ tables are un-normalised transcript counts whose within-
pair ratio is confounded by transcript length. Pairs called only in BrainSEQ regions therefore
carry structural verdicts but no expression or age verdict, and are reported as such.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from isograph.explain.structure import parse_gtf
from isograph_benchmark.paths import ensure_dir, rel, stage_out
from isograph_benchmark.real_data.interpret_modules import (DEFAULT_GTF_CACHE,
                                                            DEFAULT_GTF_PATH)

# GTEx regions carrying an RSEM TPM table. The BrainSEQ entries of the panel's region list
# have no length-normalised transcript quantification, so they cannot enter the ratio arm.
_GTEX_REGIONS = ("amygdala", "anterior_cingulate_cortex_ba24", "caudate_basal_ganglia",
                 "cerebellar_hemisphere", "cerebellum", "cortex", "frontal_cortex_ba9",
                 "hippocampus", "hypothalamus", "nucleus_accumbens_basal_ganglia",
                 "putamen_basal_ganglia", "spinal_cord_cervical_c_1", "substantia_nigra")

# GTEx technical covariates carried alongside age, matching configs/real_data.yaml.
_COVARIATES = ("SEX", "SMRIN", "SMTSISCH")

_PANEL_DIR = stage_out("integration", "rbp_target_panel")


# --------------------------------------------------------------------------- structure
def _introns(exons: list[tuple[int, int]]) -> set[tuple[int, int]]:
    """Genomic intron intervals. Strand-agnostic: we only ever compare these as sets."""
    e = sorted(exons)
    return {(e[i][1], e[i + 1][0]) for i in range(len(e) - 1)}


def _subtract(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    """Exonic bp present in `a` and absent from `b`."""
    covered = np.zeros(0, dtype=bool)
    if not a:
        return 0
    lo = min(s for s, _ in a)
    hi = max(e for _, e in a)
    covered = np.zeros(hi - lo, dtype=bool)
    for s, e in a:
        covered[s - lo:e - lo] = True
    for s, e in b:
        s, e = max(s, lo), min(e, hi)
        if e > s:
            covered[s - lo:e - lo] = False
    return int(covered.sum())


def _classify(rec1, rec2, min_unique_bp: int) -> dict:
    """Structural verdict for one pair.

    An isoform-ratio assay needs each member to carry a feature the other lacks, and there are
    exactly two kinds of feature a primer can sit on: a splice junction only that transcript
    makes, or a stretch of sequence only that transcript retains. Either one identifies its
    isoform, so the pair is discriminable when *both* members have at least one of them --
    a retained-intron pair qualifies (spliced side owns the junction, retained side owns the
    intron body) even though only one side owns a junction.

    The class then records which design the pair actually supports, because they cost
    different things at the bench. `junction` is the target: both isoforms are read with a
    junction-spanning primer. `junction_and_segment` needs one junction primer and one
    internal amplicon. `unique_segment` needs two internal amplicons and no junction primer at
    all, which is the least specific option.

    `nested_terminal` is the failure worth naming: the transcripts share every junction and one
    is a terminal extension of the other, so nothing is unique to the shorter isoform. Its
    ratio is only reachable with a proximal-versus-distal amplicon design, which measures a
    different quantity (the extension's usage), so it is not counted as measurable here.
    """
    i1, i2 = _introns(rec1.exons), _introns(rec2.exons)
    u1, u2 = i1 - i2, i2 - i1
    bp1 = _subtract(rec1.exons, rec2.exons)
    bp2 = _subtract(rec2.exons, rec1.exons)

    e1, e2 = sorted(rec1.exons), sorted(rec2.exons)
    # Can each isoform be identified on its own? A unique junction or a long enough unique
    # stretch will each do it; which one decides the assay design, not whether it is possible.
    id1 = bool(u1) or bp1 >= min_unique_bp
    id2 = bool(u2) or bp2 >= min_unique_bp

    if id1 and id2:
        if u1 and u2:
            cls = "junction"
        elif u1 or u2:
            cls = "junction_and_segment"
        else:
            cls = "unique_segment"
    elif not u1 and not u2 and (e1[0] != e2[0] or e1[-1] != e2[-1]):
        cls = "nested_terminal"
    else:
        cls = "not_discriminable"

    return {"n_exons_1": len(e1), "n_exons_2": len(e2),
            "n_unique_junctions_1": len(u1), "n_unique_junctions_2": len(u2),
            "unique_exonic_bp_1": bp1, "unique_exonic_bp_2": bp2,
            "assay_class": cls,
            "structurally_measurable": cls in ("junction", "junction_and_segment",
                                               "unique_segment")}


def structural_table(pairs: pd.DataFrame, gtf_path: Path, gtf_cache: Path,
                     min_unique_bp: int) -> pd.DataFrame:
    """Structural verdicts for every distinct transcript pair in `pairs`."""
    recs = parse_gtf(gtf_path, cache=gtf_cache)
    rows = []
    for r in pairs.itertuples():
        a, b = recs.get(r.transcript_id_1), recs.get(r.transcript_id_2)
        if a is None or b is None:
            rows.append({"transcript_id_1": r.transcript_id_1,
                         "transcript_id_2": r.transcript_id_2,
                         "assay_class": "not_in_annotation",
                         "structurally_measurable": False})
            continue
        rows.append({"transcript_id_1": r.transcript_id_1,
                     "transcript_id_2": r.transcript_id_2,
                     **_classify(a, b, min_unique_bp)})
    return pd.DataFrame(rows).drop_duplicates(["transcript_id_1", "transcript_id_2"])


# --------------------------------------------------------------------------- expression
def _load_region(region: str, wanted: set[str]) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """TPM rows for `wanted` transcripts plus the aligned sample table for one GTEx region."""
    tpm_p = rel("inputs", "processed", "gtex_v11", region, "transcript_tpm.parquet")
    samp_p = rel("inputs", "bundles", "gtex_v11_brain", region, "samples.parquet")
    if not tpm_p.exists() or not samp_p.exists():
        return None
    tpm = pd.read_parquet(tpm_p)
    tpm = tpm[tpm["transcript_id"].isin(wanted)]
    if tpm.empty:
        return None
    samples = pd.read_parquet(samp_p)
    keep = [c for c in tpm.columns if c in set(samples["sample_id"])]
    samples = samples.set_index("sample_id").loc[keep]
    return tpm.set_index("transcript_id")[keep], samples


def _fit_age(frac: np.ndarray, samples: pd.DataFrame) -> dict:
    """logit(carrier fraction) ~ age + technical covariates, on the samples where the pair is on.

    The logit keeps a bounded proportion from being fitted as if it were unbounded; the shift
    away from the open interval is the standard empirical-logit guard and is applied to every
    pair identically so it cannot favour one.
    """
    eps = 1e-3
    y = np.log(np.clip(frac, eps, 1 - eps) / (1 - np.clip(frac, eps, 1 - eps)))
    X = pd.DataFrame({"age": samples["AGE"].to_numpy(float)}, index=samples.index)
    for c in _COVARIATES:
        if c in samples.columns:
            X[c] = pd.to_numeric(samples[c], errors="coerce").to_numpy(float)
    X = X.dropna(axis=1, how="all")
    ok = np.isfinite(y) & X.notna().all(axis=1).to_numpy()
    if ok.sum() < 20 or np.nanstd(y[ok]) == 0:
        return {"age_beta": np.nan, "age_p": np.nan, "n_samples_fit": int(ok.sum())}
    fit = sm.OLS(y[ok], sm.add_constant(X[ok], has_constant="add")).fit()
    return {"age_beta": float(fit.params["age"]), "age_p": float(fit.pvalues["age"]),
            "n_samples_fit": int(ok.sum())}


def region_table(pairs: pd.DataFrame, min_tpm: float, min_frac: float) -> pd.DataFrame:
    """Per (pair, region) expression profile and age association of the carrier fraction."""
    wanted = set(pairs["transcript_id_1"]) | set(pairs["transcript_id_2"])
    rows = []
    for region in _GTEX_REGIONS:
        loaded = _load_region(region, wanted)
        if loaded is None:
            continue
        tpm, samples = loaded
        # Only score a pair in a region where IsoGraph actually called that pair.
        sub = pairs[[region in set(r.split(",")) for r in pairs["regions"]]]
        for r in sub.itertuples():
            if r.transcript_id_1 not in tpm.index or r.transcript_id_2 not in tpm.index:
                continue
            t1 = tpm.loc[r.transcript_id_1].to_numpy(float)
            t2 = tpm.loc[r.transcript_id_2].to_numpy(float)
            carrier = t1 if r.motif_carrier == r.transcript_id_1 else t2
            total = t1 + t2
            on = total >= min_tpm
            frac = np.divide(carrier, total, out=np.full_like(total, np.nan), where=total > 0)
            row = {"gene": r.gene, "transcript_id_1": r.transcript_id_1,
                   "transcript_id_2": r.transcript_id_2, "region": region,
                   "median_tpm_1": float(np.median(t1)), "median_tpm_2": float(np.median(t2)),
                   "frac_samples_on": float(on.mean()),
                   "median_carrier_frac": float(np.nanmedian(frac[on])) if on.any() else np.nan,
                   "iqr_carrier_frac": float(np.nanpercentile(frac[on], 75)
                                             - np.nanpercentile(frac[on], 25))
                   if on.sum() > 4 else np.nan}
            row.update(_fit_age(frac[on], samples[on]) if on.sum() >= 20
                       else {"age_beta": np.nan, "age_p": np.nan, "n_samples_fit": int(on.sum())})
            mid = min_frac <= (row["median_carrier_frac"] or -1) <= 1 - min_frac
            row["expression_measurable"] = bool(
                row["median_tpm_1"] >= min_tpm and row["median_tpm_2"] >= min_tpm and mid)
            rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- assembly
def summarise(struct: pd.DataFrame, byreg: pd.DataFrame, pairs: pd.DataFrame,
              min_iqr: float, q: float) -> pd.DataFrame:
    """Collapse the per-region fits to one verdict per pair.

    A pair counts as responsive if it clears both gates in *any* region it was called in: the
    perturbation is run in one cell model, and one region with a moving, measurable ratio is
    enough to justify assaying the pair there.
    """
    out = pairs.merge(struct, on=["transcript_id_1", "transcript_id_2"], how="left")
    key = ["gene", "transcript_id_1", "transcript_id_2"]
    if byreg.empty:
        out["n_regions_scored"] = 0
        for c in ("expression_measurable", "responsive"):
            out[c] = False
        for c in ("best_region", "median_carrier_frac", "iqr_carrier_frac", "age_beta",
                  "age_q", "median_tpm_1", "median_tpm_2"):
            out[c] = np.nan
        out["measurable"] = out["structurally_measurable"].fillna(False)
        return out

    b = byreg.copy()
    ok = b["age_p"].notna()
    b["age_q"] = np.nan
    if ok.any():
        b.loc[ok, "age_q"] = multipletests(b.loc[ok, "age_p"], method="fdr_bh")[1]
    b["responsive_here"] = (b["expression_measurable"] & (b["age_q"] <= q)
                            & (b["iqr_carrier_frac"] >= min_iqr)).fillna(False)

    # Report the region that best supports the assay: a responsive one first, then the widest
    # dynamic range. This is the region the pair should actually be assayed against.
    b = b.sort_values(["responsive_here", "iqr_carrier_frac"], ascending=False)
    best = b.drop_duplicates(key)
    agg = (b.groupby(key, as_index=False)
           .agg(n_regions_scored=("region", "nunique"),
                expression_measurable=("expression_measurable", "max"),
                responsive=("responsive_here", "max")))
    cols = ["region", "median_tpm_1", "median_tpm_2", "median_carrier_frac",
            "iqr_carrier_frac", "age_beta", "age_q", "n_samples_fit"]
    out = (out.merge(agg, on=key, how="left")
           .merge(best[key + cols].rename(columns={"region": "best_region"}), on=key, how="left"))
    for c in ("expression_measurable", "responsive", "structurally_measurable"):
        out[c] = out[c].fillna(False).astype(bool)
    out["n_regions_scored"] = out["n_regions_scored"].fillna(0).astype(int)
    out["measurable"] = out["structurally_measurable"] & out["expression_measurable"]
    out["verdict"] = np.where(out["measurable"] & out["responsive"], "measurable+responsive",
                     np.where(out["measurable"], "measurable_static",
                     np.where(out["structurally_measurable"], "structure_only",
                              "not_assayable")))
    return out.sort_values(["verdict", "gene", "iqr_carrier_frac"],
                           ascending=[True, True, False]).reset_index(drop=True)


def run(rbp: str, min_tpm: float, min_frac: float, min_iqr: float, q: float,
        min_unique_bp: int, gtf_path: Path, gtf_cache: Path) -> pd.DataFrame:
    src = _PANEL_DIR / rbp / "rbp_target_switch_pairs.parquet"
    if not src.exists():
        raise SystemExit(f"no frozen switch pairs for {rbp}: {src} -- run rbp_target_panel first")
    pairs = pd.read_parquet(src)
    struct = structural_table(pairs, gtf_path, gtf_cache, min_unique_bp)
    byreg = region_table(pairs, min_tpm, min_frac)
    out = summarise(struct, byreg, pairs, min_iqr, q)

    d = ensure_dir(_PANEL_DIR / rbp)
    out.to_parquet(d / "rbp_pair_assayability.parquet", index=False, compression="zstd")
    out.to_csv(d / "rbp_pair_assayability.tsv", sep="\t", index=False)
    if not byreg.empty:
        byreg.to_parquet(d / "rbp_pair_assayability_by_region.parquet", index=False,
                         compression="zstd")
        byreg.to_csv(d / "rbp_pair_assayability_by_region.tsv", sep="\t", index=False)
    n = out["verdict"].value_counts().to_dict()
    print(f"[{rbp}] {len(out)} pairs over {out['gene'].nunique()} genes: "
          + ", ".join(f"{k}={v}" for k, v in sorted(n.items())))
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--rbp", action="append", required=True,
                   help="panel RBP to score; repeat for several")
    p.add_argument("--min-tpm", type=float, default=1.0,
                   help="median TPM each isoform must reach to be callable (default 1)")
    p.add_argument("--min-frac", type=float, default=0.10,
                   help="carrier fraction must sit in [f, 1-f] (default 0.10)")
    p.add_argument("--min-iqr", type=float, default=0.05,
                   help="interquartile range of the carrier fraction required (default 0.05)")
    p.add_argument("--q", type=float, default=0.05, help="BH threshold on the age term")
    p.add_argument("--min-unique-bp", type=int, default=80,
                   help="unique exonic bp needed to place a primer/probe (default 80)")
    p.add_argument("--gtf", default=str(DEFAULT_GTF_PATH))
    p.add_argument("--gtf-cache", default=str(DEFAULT_GTF_CACHE))
    a = p.parse_args()
    for rbp in a.rbp:
        run(rbp, a.min_tpm, a.min_frac, a.min_iqr, a.q, a.min_unique_bp,
            Path(a.gtf), Path(a.gtf_cache))


if __name__ == "__main__":
    main()
