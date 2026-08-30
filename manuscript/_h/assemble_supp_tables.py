"""Assemble real-data supplementary tables from the analysis parquet ledgers.

Reads the committed result parquets under 04_module_characterization/_m and 03_module_trust/_m/stability,
emits one clean CSV per supplementary table under manuscript/_m/supp_tables/, and a
machine-checkable manifest. Presentation only; every number is copied verbatim from
the source ledgers -- regenerate, do not hand-edit.

Run: python manuscript/_h/assemble_supp_tables.py   (login node, no SLURM)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from isograph_benchmark.paths import ensure_dir, region_store, stage_out  # noqa: E402

OUT = ensure_dir(stage_out("manuscript", "supp_tables"))
QTL = stage_out("anchoring", "qtl_anchoring_meta")
BASE = stage_out("characterize", "baseline_comparison")
TRUST = stage_out("trust.stability", "module_trust")
GATE = region_store("brainseq", "caudate_sczd")


def write(df: pd.DataFrame, name: str) -> None:
    df.to_csv(OUT / name, index=False)
    print(f"  wrote {name:32s} {df.shape[0]:>4d} x {df.shape[1]}")


def round_num(df: pd.DataFrame, n: int = 4) -> pd.DataFrame:
    """Round numeric columns to n decimals, but p-values to 3 significant figures
    so small p's (e.g. 8e-6) survive instead of collapsing to 0.0."""
    df = df.copy()
    p_like = ("p_value", "p_fe", "p_re", "pheno_fdr")
    for col in df.select_dtypes("number").columns:
        if any(tok in col for tok in p_like):
            df[col] = df[col].apply(
                lambda x: float(f"{x:.3g}") if pd.notna(x) else x)
        else:
            df[col] = df[col].round(n)
    return df


# --- S1 / S2  three-baseline comparison ------------------------------------
def baseline_tables() -> None:
    pooled = pd.read_parquet(BASE / "baseline_comparison_pooled.parquet")
    pooled = pooled[
        ["method", "features", "n_regions", "median_n_modules", "median_module_size",
         "mean_frac_pheno_sig", "mean_frac_both", "mean_frac_go_enriched",
         "total_pheno_sig", "total_both"]
    ]
    write(round_num(pooled), "tableS1_baseline_pooled.csv")

    per = pd.read_parquet(BASE / "baseline_comparison.parquet")
    per = per.sort_values(["cohort", "region", "method"]).reset_index(drop=True)
    write(round_num(per), "tableS2_baseline_per_region.csv")


# --- S3 / S4  QTL splicing-specificity contrast ----------------------------
def _tidy_contrast(df: pd.DataFrame) -> pd.DataFrame:
    df = df[["graph_method", "module_set", "k", "ratio_fe", "ratio_fe_low",
             "ratio_fe_high", "p_fe", "I2"]].copy()
    df.columns = ["method", "module_set", "n_analyses", "sqtl_eqtl_ratio",
                  "ratio_ci_low", "ratio_ci_high", "p_value", "I2"]
    order = {"all_modules": 0, "pheno_sig_modules": 1, "go_invisible_modules": 2,
             "go_visible_modules": 3}
    mord = {"isograph": 0, "wgcna_switch_only": 1, "wgcna_multiplex": 2}
    df["_o"] = df.module_set.map(order)
    df["_m"] = df.method.map(mord)
    df = df.sort_values(["_o", "_m"]).drop(columns=["_o", "_m"]).reset_index(drop=True)
    return df


def qtl_tables() -> None:
    contrast = _tidy_contrast(pd.read_parquet(QTL / "qtl_anchoring_meta_contrast.parquet"))
    write(round_num(contrast), "tableS3_qtl_specificity_contrast.csv")

    common = _tidy_contrast(pd.read_parquet(QTL / "qtl_anchoring_meta_contrast_common.parquet"))
    write(round_num(common), "tableS4_qtl_specificity_matched_baseline.csv")

    raw = pd.read_parquet(QTL / "qtl_anchoring_meta.parquet")
    raw = raw[raw.k > 0][["graph_method", "xqtl_kind", "module_set", "k",
                          "n_fg_total", "or_fe", "or_fe_low", "or_fe_high",
                          "p_fe", "I2"]].copy()
    raw.columns = ["method", "xqtl_kind", "module_set", "n_analyses", "n_fg_total",
                   "or", "or_ci_low", "or_ci_high", "p_value", "I2"]
    write(round_num(raw), "tableS5_qtl_raw_enrichment_or.csv")


# --- S6  GO-invisible disease switch modules -------------------------------
def gate_table() -> None:
    import json
    g = pd.read_parquet(GATE / "go_invisible_gate.parquet").copy()
    bg = json.load(open(GATE / "go_invisible_gate_background.json"))
    cols = ["module", "n_genes", "pheno_fdr", "go_invisible", "n_go_terms",
            "genes_with_real_switch", "max_switch_strength", "n_sig_switch_tx",
            "frac_cds_changed", "frac_coding_status_change", "frac_biotype_switch",
            "frac_utr_changed", "top_switch_genes"]
    g = g[cols]
    bg_row = {c: bg.get(c, "") for c in cols}
    bg_row["module"] = "_background (all switch tx)"
    g = pd.concat([g, pd.DataFrame([bg_row])], ignore_index=True)
    write(round_num(g), "tableS6_go_invisible_gate.csv")


# --- S7  per-region module trust funnel ------------------------------------
def trust_table() -> None:
    rows = []
    for stab in sorted(TRUST.iterdir()):
        if not stab.name.startswith("module_stability__") or "isograph" not in stab.name:
            continue
        # module_stability__<cohort>__<region>__isograph.parquet
        parts = stab.stem.split("__")
        cohort, region = parts[1], parts[2]
        s = pd.read_parquet(stab)
        n_mod = len(s)
        n_trust = int(s.trusted.sum())

        wfile = TRUST / f"within_cohort__{cohort}__{region}__isograph.parquet"
        rho_med = pos_frac = None
        if wfile.exists():
            w = pd.read_parquet(wfile)
            rho = w.driver_load_rho.dropna()
            if len(rho):
                rho_med = float(rho.median())
                pos_frac = float((rho > 0).mean())

        cfile = TRUST / f"module_complementarity__{cohort}__{region}__isograph.parquet"
        dtu_med = wgage_med = None
        if cfile.exists():
            c = pd.read_parquet(cfile)
            dtu_med = float(c.frac_dtu_without_dge.median())
            wgage_med = float(c.frac_in_wgcna_age_modules.median())

        n_pairs = n_rep = None
        if cohort == "brainseq":
            # replication keys BrainSEQ DLPFC as dlpfc_ba9 (paired GTEx region name)
            rep_region = "dlpfc_ba9" if region == "dlpfc" else region
            rfile = TRUST / f"module_aging_replication__{rep_region}__isograph.parquet"
            if rfile.exists():
                r = pd.read_parquet(rfile)
                n_pairs = len(r)
                n_rep = int(r.replicates.sum())

        rows.append(dict(
            cohort=cohort, region=region,
            n_modules=n_mod, n_trusted=n_trust,
            frac_trusted=round(n_trust / n_mod, 3) if n_mod else None,
            median_driver_rho=round(rho_med, 3) if rho_med is not None else None,
            frac_positive_rho=round(pos_frac, 3) if pos_frac is not None else None,
            n_replication_pairs=n_pairs, n_concordant=n_rep,
            median_frac_dtu_without_dge=round(dtu_med, 4) if dtu_med is not None else None,
            median_frac_in_wgcna_age=round(wgage_med, 3) if wgage_med is not None else None,
        ))
    df = pd.DataFrame(rows).sort_values(["cohort", "region"]).reset_index(drop=True)
    write(df, "tableS7_module_trust_funnel.csv")


def main() -> None:
    print(f"Writing supplementary tables to {OUT}")
    baseline_tables()
    qtl_tables()
    gate_table()
    trust_table()
    print("done.")


if __name__ == "__main__":
    main()
