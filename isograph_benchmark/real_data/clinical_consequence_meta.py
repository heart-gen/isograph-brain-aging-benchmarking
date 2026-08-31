"""Cross-region rollup of the switched-exon clinical-consequence test.

Aggregates each region's `clinical_consequence/` outputs into a per-stratum summary: how many
regions show the switched/background ClinVar P/LP density ratio > 1 at empirical p < 0.05, the
median ratio, a Fisher-combined permutation p, and the pooled gnomAD LOEUF contrast. Answers
whether switched exons carry more disease signal than constitutive exons across the switch layer.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues

from isograph_benchmark.paths import cohort_dir, stage_out

_ROOTS = [cohort_dir("brainseq"), cohort_dir("gtex")]


def _collect(name: str) -> pd.DataFrame:
    frames = []
    for root in _ROOTS:
        for f in root.glob(f"*/_m/isograph_vae/clinical_consequence/{name}.parquet"):
            frames.append(pd.read_parquet(f))
    if not frames:
        raise SystemExit(f"no per-region {name}.parquet found; run clinical_consequence first.")
    return pd.concat(frames, ignore_index=True)


def run() -> pd.DataFrame:
    cc = _collect("clinical_consequence")
    con = _collect("constraint_summary")
    scope_col = "scope" if "scope" in cc.columns else None
    keys = ["stratum", "scope"] if scope_col else ["stratum"]
    rows = []
    for key, sub in cc.groupby(keys):
        stratum = key[0] if isinstance(key, tuple) else key
        scope = key[1] if isinstance(key, tuple) else "all_exons"
        sub = sub.dropna(subset=["ratio", "p_emp"])
        if sub.empty:
            continue
        p = np.clip(sub["p_emp"].to_numpy(), 1e-6, 1.0)
        _, p_comb = combine_pvalues(p, method="fisher")
        cs = con[con["stratum"] == stratum].dropna(subset=["mwu_p_more_constrained"])
        _, loeuf_p = (combine_pvalues(np.clip(cs["mwu_p_more_constrained"], 1e-12, 1.0),
                                      method="fisher") if not cs.empty else (np.nan, np.nan))
        rows.append({
            "stratum": stratum,
            "scope": scope,
            "n_regions": int(sub["region"].nunique()),
            "n_ratio_gt1_p05": int(((sub["ratio"] > 1) & (sub["p_emp"] < 0.05)).sum()),
            "n_ratio_lt1_p05": int(((sub["ratio"] < 1) & (sub["p_emp"] < 0.05)).sum()),
            "median_ratio": float(sub["ratio"].median()),
            "median_switched_per_kb": float(sub["switched_plp_per_kb"].median()),
            "median_bg_per_kb": float(sub["bg_plp_per_kb"].median()),
            "fisher_p": float(p_comb),
            "median_loeuf_switch": float(cs["median_loeuf_switch"].median()) if not cs.empty else np.nan,
            "loeuf_fisher_p": float(loeuf_p) if cs is not None and not cs.empty else np.nan,
        })
    meta = pd.DataFrame(rows).sort_values(["scope", "stratum"] if "scope" in
                                          pd.DataFrame(rows).columns else ["stratum"])
    out = stage_out("mechanism", "clinical_consequence_meta.parquet")
    out.parent.mkdir(parents=True, exist_ok=True)
    meta.to_parquet(out, index=False)
    _write_report(meta)
    print(f"clinical-consequence meta over {cc['region'].nunique()} regions -> {out}")
    return meta


def _write_report(meta: pd.DataFrame) -> None:
    lines = [
        "# Switched-exon clinical consequence — cross-region rollup", "",
        "**Primary anchor = gnomAD LOEUF**: are switch genes more loss-of-function constrained "
        "than genome-wide (lower median LOEUF; Fisher-combined MWU p). The exon-level ClinVar "
        "columns are a direction-neutral secondary readout: how many regions have switched exons "
        "with higher (ratio>1) vs lower (ratio<1) P/LP density than constitutive exons at "
        "two-sided within-gene permutation p < 0.05, and the median ratio. `scope` = all "
        "switch-pair exons vs coding (CDS-overlapping) exons only. Alt-spliced exons are usually "
        "less constrained, so ratio < 1 is the expected baseline; the CDS scope is the fairer "
        "coding-vs-coding contrast.", "",
        "| scope | stratum | regions | median LOEUF (switch) | LOEUF p | ratio>1 (p<.05) | ratio<1 (p<.05) | median ratio | Fisher p |",
        "|-------|---------|---------|-----------------------|---------|-----------------|-----------------|--------------|----------|",
    ]
    for r in meta.itertuples():
        lines.append(
            f"| {getattr(r, 'scope', 'all_exons')} | {r.stratum} | {r.n_regions} | "
            f"{r.median_loeuf_switch:.3f} | {r.loeuf_fisher_p:.2e} | {r.n_ratio_gt1_p05} | "
            f"{r.n_ratio_lt1_p05} | {r.median_ratio:.2f} | {r.fisher_p:.2e} |")
    (stage_out("mechanism", "CLINICAL_CONSEQUENCE_META.md")).write_text("\n".join(lines) + "\n")


def main() -> None:
    argparse.ArgumentParser(description="Cross-region clinical-consequence rollup.").parse_args()
    run()


if __name__ == "__main__":
    main()
