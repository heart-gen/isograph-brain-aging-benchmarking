"""S-LDSC summary: parse partitioned-heritability results into a tidy table + writeup.

Reads every `<trait>_<model>.results` under 05_genetic_anchoring/_m/ldsc/<annotation>/results/ and
assembles the genetic-anchoring-beyond-MAGMA statistics for both cases:

  * disease  — SCZ heritability partitioned on the SCZD switch layer's QTL annotations.
  * aging    — AD / PD / LBD / ALS[/ FTD] heritability partitioned on the pooled aging
               switch layer's QTL annotations.

Four models per trait, each on baselineLD v2.2:
  sqtl_only / eqtl_only / cis_only — single-annotation Finucane enrichment + coefficient
    (tau) conditional on baselineLD. This is the size-robust "is the switch-QTL layer
    genetically anchored" test (the step beyond MAGMA's size-confounded gene-set test).
  joint — baseline + sqtl_switch + eqtl_switch, for the head-to-head splicing-vs-
    expression conditional contrast.

The coefficient z-score gives a one-sided p (tau > 0). Writes, under 05_genetic_anchoring/_m/ldsc/:
  ldsc_partitioned.parquet — tidy (trait, case, annotation, model, annot, prop_snps,
                             prop_h2, enrichment, enrichment_p, coef, coef_se, coef_z,
                             coef_p)
  LDSC_SUMMARY.md          — Manubot writeup
"""
from __future__ import annotations

import argparse

import pandas as pd
from scipy.stats import norm

from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data.gwas_traits import TRAITS

# custom annotations appended (in order) per model
_MODEL_ANNOTS = {
    "sqtl_only": ["sqtl_switch"],
    "eqtl_only": ["eqtl_switch"],
    "cis_only": ["cis_switch"],
    "joint": ["sqtl_switch", "eqtl_switch"],
}
# annotation dir -> case label
_CASE = {"brainseq-sczd": "disease", "aging": "aging"}


def _parse(path, trait, case, annotation, model, annots) -> list[dict]:
    d = pd.read_csv(path, sep="\t")
    tail = d.tail(len(annots)).reset_index(drop=True)  # appended custom annotations
    rows = []
    for i, annot in enumerate(annots):
        r = tail.iloc[i]
        z = float(r["Coefficient_z-score"])
        rows.append({
            "trait": trait, "case": case, "annotation": annotation, "model": model,
            "annot": annot,
            "prop_snps": float(r["Prop._SNPs"]),
            "prop_h2": float(r["Prop._h2"]),
            "enrichment": float(r["Enrichment"]),
            "enrichment_p": float(r["Enrichment_p"]),
            "coef": float(r["Coefficient"]),
            "coef_se": float(r["Coefficient_std_error"]),
            "coef_z": z,
            "coef_p": float(norm.sf(z)),
        })
    return rows


def collect() -> pd.DataFrame:
    base = stage_out("anchoring.ldsc")
    rows = []
    for annotation, case in _CASE.items():
        rdir = base / annotation / "results"
        if not rdir.exists():
            continue
        for trait in TRAITS:
            for model, annots in _MODEL_ANNOTS.items():
                f = rdir / f"{trait}_{model}.results"
                if f.exists():
                    rows.append(pd.DataFrame(_parse(f, trait, case, annotation, model, annots)))
    if not rows:
        raise SystemExit("No S-LDSC .results files found; run 05_genetic_anchoring/_h/03d.ldsc_munge_h2.sh first.")
    return pd.concat(rows, ignore_index=True)


def _fmt(df, model, annot, col):
    r = df[(df.model == model) & (df.annot == annot)]
    return r[col].iloc[0] if len(r) else float("nan")


def _write_report(tidy: pd.DataFrame, out_dir) -> None:
    lines = [
        "# Stratified LD-score regression: genetic anchoring of the switch layer",
        "",
        "Partitioned SNP-heritability of the switch layer's GTEx brain sQTL / eQTL "
        "annotations on top of baselineLD v2.2 — the size-robust step beyond MAGMA. "
        "`coef_p` is the one-sided p that the annotation's per-SNP heritability "
        "coefficient (tau) exceeds 0, conditional on baselineLD (single-annotation "
        "models) or on baselineLD + the other QTL layer (joint model).",
        "",
    ]
    for annotation, case in _CASE.items():
        sub = tidy[tidy.annotation == annotation]
        if sub.empty:
            continue
        lines += [f"## {case.capitalize()} case (annotation: {annotation})", ""]
        cols = ["trait", "model", "annot", "prop_snps", "enrichment", "enrichment_p",
                "coef_z", "coef_p"]
        show = sub[cols].copy()
        for c in ["prop_snps", "enrichment", "enrichment_p", "coef_z", "coef_p"]:
            show[c] = show[c].map(lambda x: f"{x:.3g}")
        hdr = "| " + " | ".join(cols) + " |"
        sep = "| " + " | ".join("---" for _ in cols) + " |"
        body = ["| " + " | ".join(str(getattr(r, c)) for c in cols) + " |"
                for r in show.itertuples(index=False)]
        lines += [hdr, sep, *body, ""]

    # headline reading (single-annotation coefficient tests)
    lines += ["## Reading", ""]
    for annotation, case in _CASE.items():
        sub = tidy[tidy.annotation == annotation]
        for trait in sub.trait.unique():
            t = sub[sub.trait == trait]
            label = TRAITS[trait].label
            sq_e = _fmt(t, "sqtl_only", "sqtl_switch", "enrichment")
            sq_p = _fmt(t, "sqtl_only", "sqtl_switch", "coef_p")
            eq_e = _fmt(t, "eqtl_only", "eqtl_switch", "enrichment")
            eq_p = _fmt(t, "eqtl_only", "eqtl_switch", "coef_p")
            ci_e = _fmt(t, "cis_only", "cis_switch", "enrichment")
            ci_p = _fmt(t, "cis_only", "cis_switch", "coef_p")
            js_p = _fmt(t, "joint", "sqtl_switch", "coef_p")
            je_p = _fmt(t, "joint", "eqtl_switch", "coef_p")
            lines.append(
                f"- **{label}** ({case}): switch-layer sQTL enrichment {sq_e:.2f}× "
                f"(tau p={sq_p:.3g}); eQTL {eq_e:.2f}× (p={eq_p:.3g}); total cis "
                f"{ci_e:.2f}× (p={ci_p:.3g}). Joint sQTL-vs-eQTL conditional tau p: "
                f"sQTL={js_p:.3g}, eQTL={je_p:.3g}.")
    lines += [
        "",
        "S-LDSC is robust to the gene-size confound that inflates MAGMA on giant "
        "modules, so a positive coefficient is heritability genuinely concentrated in "
        "the switch layer's cis-regulatory variants, not an artifact of annotation "
        "size. Because a gene's sQTL and eQTL SNPs overlap, the joint model splits "
        "signal between them and understates each; the single-annotation coefficients "
        "are the primary enrichment test and the joint model is only the head-to-head "
        "contrast. GTEx brain QTLs are bulk-tissue, so cell-type-specific splicing "
        "(e.g. microglial for AD, dopaminergic for PD) is under-sampled — enrichment "
        "here is a floor, not a ceiling.",
    ]
    (out_dir / "LDSC_SUMMARY.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    argparse.ArgumentParser(description="Summarize S-LDSC partitioned heritability.").parse_args()
    out_dir = ensure_dir(stage_out("anchoring.ldsc"))
    tidy = collect()
    tidy.to_parquet(out_dir / "ldsc_partitioned.parquet", index=False, compression="zstd")
    _write_report(tidy, out_dir)
    print(tidy.to_string(index=False))
    print(f"\nWrote {out_dir/'ldsc_partitioned.parquet'} and LDSC_SUMMARY.md")


if __name__ == "__main__":
    main()
