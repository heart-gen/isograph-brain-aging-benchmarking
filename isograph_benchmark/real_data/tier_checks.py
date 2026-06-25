"""Six validation checks comparing the IsoGraph multiplex tiers (pilot + fan-out).

Consumes the per-tier artifacts written by ``project_tiers`` (modules / edges /
module_gene_roles) plus the shared fit's ``feature_scores`` and the WGCNA baseline,
and emits one table per check so the full-multiplex promotion can be defended as a
characterization rather than asserted. feature_scores is identical across the three
IsoGraph tiers (one VAE fit, three edge projections), so it is loaded once from the
source fit dir while ``modules``/``roles`` are read per tier.

Checks (per the promotion brief):
  1. Channel composition  -- switch vs abundance role mix per module (are modules
                             genuinely multiplex or just one channel re-badged?).
  2. WGCNA overlap        -- best-match Jaccard of each tier module against the
                             gene-abundance baseline (novel structure vs re-derivation).
  3. Driver transcripts   -- top switch drivers per module (interpretability).
  4. QC / degradation     -- module eigengene vs RIN / mito / rRNA / PMI / TIN; shows
                             modules are not simply degradation/3'-bias proxies.
  5. Trait decomposition  -- module-trait association attributable to switching,
                             abundance, or both (reuses incremental_association).
  6. Ablation roll-up     -- tier_projection_summary + per-check tallies, i.e. WHY the
                             full multiplex is selected over the switch-restricted tiers.

Usage (project root, isograph env):
    python -m isograph_benchmark.real_data.tier_checks brainseq-aging --region caudate
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from isograph.io.artifacts import load_dataset_bundle
from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.incremental_association import (
    _channel_matrix,
    _eigengene_table,
    _sample_cols,
    module_level_incremental,
)
from isograph_benchmark.real_data.project_tiers import (
    TIER_DIRS,
    PRIMARY_TIER,
    _bundle_path,
    _region_dir,
)

# Degradation / library-quality proxies. A module that merely tracks 3' degradation
# would correlate strongly with these; we want trusted modules to NOT.
QC_PROXIES = ["RIN", "mito_rate", "r_rna_rate", "mapping_rate", "PMI"]
DRIVER_TOPK = 15
QC_ALIGN_THRESHOLD = 0.5  # |corr| above which a module looks degradation-driven


def _tin_median(analysis: str, region: str | None) -> pd.Series | None:
    """Per-sample median TIN (a direct 3'-degradation proxy) if extracted."""
    if analysis != "brainseq-aging" or region is None:
        return None
    p = rel("inputs", "tin", f"brainseq__{region}__sample_median_tin.csv")
    if not p.exists():
        return None
    df = pd.read_csv(p)
    scol = "sample_id" if "sample_id" in df.columns else df.columns[0]
    vcol = [c for c in df.columns if c != scol][0]
    return df.set_index(df[scol].astype(str))[vcol].astype(float)


# --------------------------------------------------------------------------- #
# 1. Channel composition
# --------------------------------------------------------------------------- #
def check_channel_composition(roles: pd.DataFrame, tier: str) -> pd.DataFrame:
    if roles.empty:
        return pd.DataFrame()
    rows = []
    for mid, g in roles.groupby("module_id"):
        sw = int(g["switch_active"].sum()) if "switch_active" in g else 0
        ab = int(g["abundance_active"].sum()) if "abundance_active" in g else 0
        dual = int((g.get("switch_active", False) & g.get("abundance_active", False)).sum())
        n = len(g)
        rows.append({
            "tier": tier, "module_id": mid, "n_genes": n,
            "n_switch_active": sw, "n_abundance_active": ab, "n_dual_active": dual,
            "switch_frac": round(sw / n, 4) if n else 0.0,
            "abundance_frac": round(ab / n, 4) if n else 0.0,
            "dual_frac": round(dual / n, 4) if n else 0.0,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# 2. WGCNA overlap
# --------------------------------------------------------------------------- #
def check_wgcna_overlap(modules: pd.DataFrame, wgcna: pd.DataFrame, tier: str) -> pd.DataFrame:
    if modules.empty or wgcna is None or wgcna.empty:
        return pd.DataFrame()
    wg = {mid: set(g["gene_id"]) for mid, g in wgcna.groupby("module_id")}
    rows = []
    for mid, g in modules.groupby("module_id"):
        genes = set(g["gene_id"])
        best_j, best_w, best_ov = 0.0, None, 0
        for wmid, wgenes in wg.items():
            inter = len(genes & wgenes)
            if inter == 0:
                continue
            j = inter / len(genes | wgenes)
            if j > best_j:
                best_j, best_w, best_ov = j, wmid, inter
        rows.append({
            "tier": tier, "module_id": mid, "n_genes": len(genes),
            "best_wgcna_module": best_w, "best_jaccard": round(best_j, 4),
            "overlap_genes": best_ov,
            "frac_in_best_wgcna": round(best_ov / len(genes), 4) if genes else 0.0,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# 3. Driver transcripts
# --------------------------------------------------------------------------- #
def check_drivers(roles: pd.DataFrame, tier: str, topk: int = DRIVER_TOPK) -> pd.DataFrame:
    if roles.empty or "switch_r" not in roles:
        return pd.DataFrame()
    sw = roles[roles.get("switch_active", False)].copy()
    if sw.empty:
        return pd.DataFrame()
    sw["abs_switch_r"] = sw["switch_r"].abs()
    out = []
    for mid, g in sw.groupby("module_id"):
        top = g.sort_values("abs_switch_r", ascending=False).head(topk)
        for rank, (_, r) in enumerate(top.iterrows(), 1):
            out.append({
                "tier": tier, "module_id": mid, "rank": rank,
                "gene_id": r["gene_id"], "switch_r": round(float(r["switch_r"]), 4),
                "module_role": r.get("module_role", None),
            })
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- #
# 4. QC / degradation alignment
# --------------------------------------------------------------------------- #
def check_qc_alignment(fs: pd.DataFrame, modules: pd.DataFrame, bundle,
                       tin: pd.Series | None, tier: str) -> pd.DataFrame:
    if modules.empty:
        return pd.DataFrame()
    st = bundle.sample_table
    sample_ids = set(st["sample_id"].astype(str))
    samp = _sample_cols(fs, sample_ids)
    sti = st.set_index(st["sample_id"].astype(str)).loc[samp]

    proxies = {}
    for c in QC_PROXIES:
        if c in sti.columns:
            proxies[c] = pd.to_numeric(sti[c], errors="coerce").to_numpy(float)
    if tin is not None:
        proxies["median_TIN"] = tin.reindex(samp).to_numpy(float)
    if not proxies:
        return pd.DataFrame()

    # _eigengene_table -> module_id-indexed rows with one column per sample.
    eg_sw = _eigengene_table(_channel_matrix(fs, "switch", samp), modules, samp).set_index("module_id")[samp]
    eg_ab = _eigengene_table(_channel_matrix(fs, "abundance", samp), modules, samp).set_index("module_id")[samp]

    def corr(y, x):
        m = np.isfinite(y) & np.isfinite(x)
        if m.sum() < 10 or np.std(y[m]) == 0 or np.std(x[m]) == 0:
            return np.nan
        return float(np.corrcoef(y[m], x[m])[0, 1])

    rows = []
    for mid in modules["module_id"].unique():
        ys = eg_sw.loc[mid].to_numpy(float) if mid in eg_sw.index else None
        ya = eg_ab.loc[mid].to_numpy(float) if mid in eg_ab.index else None
        rec = {"tier": tier, "module_id": mid}
        max_abs = 0.0
        worst = None
        for name, x in proxies.items():
            rs = corr(ys, x) if ys is not None else np.nan
            ra = corr(ya, x) if ya is not None else np.nan
            best = max([abs(v) for v in (rs, ra) if np.isfinite(v)], default=np.nan)
            rec[f"r_switch_{name}"] = round(rs, 4) if np.isfinite(rs) else np.nan
            rec[f"r_abund_{name}"] = round(ra, 4) if np.isfinite(ra) else np.nan
            if np.isfinite(best) and best > max_abs:
                max_abs, worst = best, name
        rec["max_abs_qc_corr"] = round(max_abs, 4)
        rec["worst_proxy"] = worst
        rec["degradation_flag"] = bool(max_abs >= QC_ALIGN_THRESHOLD)
        rows.append(rec)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def _load_tier(rdir: Path, tier: str):
    d = rdir / TIER_DIRS[tier]
    mods = pd.read_parquet(d / "modules.parquet") if (d / "modules.parquet").exists() else pd.DataFrame()
    rp = d / "module_gene_roles.parquet"
    roles = pd.read_parquet(rp) if rp.exists() else pd.DataFrame()
    return mods, roles


def run(analysis: str, region: str | None, source_subdir: str = "isograph_vae") -> dict:
    label = f"{analysis}/{region}" if region else analysis
    rdir = _region_dir(analysis, region)
    src = rdir / source_subdir

    fs = pd.read_parquet(src / "feature_scores.parquet")
    bundle = load_dataset_bundle(_bundle_path(analysis, region))
    tin = _tin_median(analysis, region)
    wgcna_path = rdir / "wgcna_gene" / "modules.parquet"
    wgcna = pd.read_parquet(wgcna_path) if wgcna_path.exists() else None

    out = ensure_dir(rdir / "tier_checks")
    c1, c2, c3, c4, c5 = [], [], [], [], []
    headline = {}

    for tier in TIER_DIRS:
        mods, roles = _load_tier(rdir, tier)
        if mods.empty:
            print(f"[{label}] WARN: {tier} has no modules.parquet -- skipping", flush=True)
            continue
        comp = check_channel_composition(roles, tier)
        ov = check_wgcna_overlap(mods, wgcna, tier)
        dr = check_drivers(roles, tier)
        qc = check_qc_alignment(fs, mods, bundle, tin, tier)
        dec = module_level_incremental(analysis, fs, mods, bundle)
        dec.insert(0, "tier", tier)

        c1.append(comp); c2.append(ov); c3.append(dr); c4.append(qc); c5.append(dec)
        cat = dec["category"].value_counts().to_dict() if not dec.empty else {}
        headline[tier] = {
            "n_modules": int(mods["module_id"].nunique()),
            "median_switch_frac": round(float(comp["switch_frac"].median()), 3) if not comp.empty else None,
            "median_best_jaccard_vs_wgcna": round(float(ov["best_jaccard"].median()), 3) if not ov.empty else None,
            "n_degradation_flagged": int(qc["degradation_flag"].sum()) if not qc.empty else 0,
            "n_qc_tested": int(len(qc)),
            "trait_composition_unique": int(cat.get("composition_unique", 0)),
            "trait_abundance_unique": int(cat.get("abundance_unique", 0)),
            "trait_both": int(cat.get("both", 0)),
            "trait_neither": int(cat.get("neither", 0)),
        }
        h = headline[tier]
        print(f"[{label}] {tier}: modules={h['n_modules']} sw_frac={h['median_switch_frac']} "
              f"medJ_wgcna={h['median_best_jaccard_vs_wgcna']} deg_flag={h['n_degradation_flagged']}/{h['n_qc_tested']} "
              f"trait[sw={h['trait_composition_unique']} ab={h['trait_abundance_unique']} both={h['trait_both']}]",
              flush=True)

    def _concat_write(parts, name):
        df = pd.concat([p for p in parts if not p.empty], ignore_index=True) if any(not p.empty for p in parts) else pd.DataFrame()
        if not df.empty:
            df.to_parquet(out / name, index=False, compression="zstd")
        return df

    _concat_write(c1, "check1_channel_composition.parquet")
    _concat_write(c2, "check2_wgcna_overlap.parquet")
    _concat_write(c3, "check3_drivers.parquet")
    _concat_write(c4, "check4_qc_alignment.parquet")
    _concat_write(c5, "check5_trait_decomposition.parquet")

    # 6. Ablation roll-up: tier projection summary + per-check tallies.
    tps_path = rdir / "tier_projection_summary.parquet"
    ablation = pd.read_parquet(tps_path) if tps_path.exists() else pd.DataFrame()
    summary = {"analysis": analysis, "region": region or "caudate_sczd",
               "primary_tier": PRIMARY_TIER, "tiers": headline}
    (out / "tier_checks_summary.json").write_text(json.dumps(summary, indent=2))
    if not ablation.empty:
        ablation.to_parquet(out / "check6_ablation_summary.parquet", index=False, compression="zstd")

    print(f"\n[{label}] tier checks written to {out}")
    return summary


def main() -> None:
    p = argparse.ArgumentParser(description="Six validation checks across IsoGraph tiers.")
    p.add_argument("analysis", choices=["brainseq-sczd", "brainseq-aging", "gtex-aging"])
    p.add_argument("--region", default=None)
    p.add_argument("--source-subdir", default="isograph_vae")
    args = p.parse_args()
    run(args.analysis, args.region, source_subdir=args.source_subdir)


if __name__ == "__main__":
    main()
