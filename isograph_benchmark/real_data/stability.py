"""Within-cohort split-half stability test (method reproducibility, no cross-quantifier confound).

Separates the two explanations for IsoGraph's weak BrainSEQ<->GTEx replication:
  (a) isoform-switch co-expression is genuinely low-reproducibility biology, vs.
  (b) the cross-cohort test is confounded by the Salmon (BrainSEQ) vs RSEM (GTEx)
      transcript-quantifier mismatch, which scrambles within-gene isoform ratios
      (the switch signal) while leaving gene abundance — and thus WGCNA — intact.

Within ONE cohort+region the quantifier and preprocessing are fixed, so a split-half
refit isolates pure estimation/clustering stability. We randomly partition the
samples 50/50 (several seeds), refit the SAME method on each half with the SAME
feature set, and measure partition agreement (ARI, NMI) over the genes assigned in
both halves — the same metrics used for the cross-cohort comparison.

Interpretation:
  IsoGraph within-cohort ARI >> its cross-cohort ARI  -> the cross-cohort failure is
      the quantifier/preprocessing confound (the method is stable given consistent
      input); the weak BrainSEQ<->GTEx switch replication is a data-comparability
      limitation, not a method flaw.
  IsoGraph within-cohort ARI ~ cross-cohort ARI (~0)  -> method instability at this
      sample size (the switch network is hard to estimate), independent of quantifier.
  WGCNA within-cohort ARI gives the abundance-network reference ceiling.

Subcommands:
  fit-isograph  fit IsoGraph on both halves for each seed; write per-half partitions
  aggregate     read all per-half partitions (both methods) -> ARI/NMI -> summary

WGCNA partitions are produced by ``real_data/stability/_h/stability_wgcna.R`` (R/WGCNA)
and written to the same ``partitions/`` dir, so ``aggregate`` treats both methods
uniformly. Splits are drawn independently per method with seed = k (k = 0..seeds-1);
we compare the mean +/- SD over seeds, not a paired per-split test, so identical
cross-method splits are not required.

Usage::

    python -m isograph_benchmark.real_data.stability fit-isograph --cohort brainseq --region caudate --seeds 5
    python -m isograph_benchmark.real_data.stability aggregate
"""
from __future__ import annotations

import argparse
import gc
import json
import math
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from isograph.io.artifacts import load_dataset_bundle
from isograph.models.vae import VaeNetworkModel
from isograph.workflow.config import VaeModelConfig
from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.run_models import _filter_expressed_transcripts
from isograph_benchmark.real_data.replication import REGION_PAIRS

SEED_BASE = 1000  # split seeds are SEED_BASE + k; VAE init seed is fixed (below)
VAE_SEED = 13     # fixed across halves: variation comes from the sample split, not init

# Per-cohort fitting spec, mirroring run_models.run_brainseq_region / run_gtex_region.
# Keep these in sync with run_models if the production configs change. "covariates" is
# the discovery-only residualization set (mirrors run_models.*_DISCOVERY_COVARIATES):
# technical/topology confounds only, so the split-half modules match the shipped fit.
COHORTS = {
    "brainseq": {
        "bundle_root": ("inputs", "bundles", "brainseq_v1"),
        "regions": ["caudate", "hippocampus", "dlpfc"],
        "covariates": ["RIN", "mapping_rate", "mito_rate",
                       "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5"],
        "age_col": "Age",
        "filter_transcripts": True,   # run_brainseq_region filters; GTEx bundles are pre-filtered
        "lr": None,                   # default (1e-3)
    },
    "gtex": {
        "bundle_root": ("inputs", "bundles", "gtex_v11_brain"),
        "regions": ["caudate_basal_ganglia", "hippocampus", "frontal_cortex_ba9"],
        "covariates": ["SMRIN", "SMTSISCH", "SMMAPRT"],
        "age_col": "AGE",
        "filter_transcripts": False,
        "lr": None,                   # promoted single LR: default 1e-3 + grad_clip_norm=1.0
                                      # (gate B.2) replaces the old hand-tuned GTEx lr=3e-4.
    },
}


def _partitions_dir():
    # Optional sandbox override so an A/B can be run in isolation (fresh baseline +
    # candidate arms) without overwriting the committed production partitions. The
    # aggregate output dir is derived from this dir's parent, so a sandbox stays
    # self-contained. Unset -> the default production location (behavior unchanged).
    override = os.environ.get("STABILITY_PARTITIONS_DIR")
    if override:
        return ensure_dir(Path(override))
    return ensure_dir(rel("real_data", "stability", "_m", "partitions"))


def _split_indices(n: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Random 50/50 partition of sample positions (half A gets floor(n/2))."""
    perm = np.random.default_rng(seed).permutation(n)
    h = n // 2
    return perm[:h], perm[h:]


def _vae_config(spec: dict, consensus_runs: int = 1, reliability: bool = False,
                min_minor_usage: float = 0.1,
                reliability_floor: float = 0.0,
                extra_covariates: list[str] | None = None,
                max_module_frac: float | None = None,
                leiden_giant_frac: float | None = None,
                lr_override: float | None = None,
                grad_clip_norm: float | None = None) -> VaeModelConfig:
    covariates = list(spec["covariates"]) + list(extra_covariates or [])
    kw = dict(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        residualize_covariates=covariates,
        min_module_size=20, trait_columns=[spec["age_col"]], random_state=VAE_SEED,
        allow_abundance_abundance=False, alpha_switch=0.5, leiden_resolution=5.0,
        # Promoted production defaults (mirror run_models._PROMOTED_VAE, 2026-06-24):
        # the split-half baseline must equal the shipped config so the trust funnel
        # validates the ACTUAL production modules. grad_clip_norm=1.0 (single-LR gate
        # B.2) + estimability switch-reliability (the only positive stability lever).
        grad_clip_norm=1.0,
        switch_reliability_weighting=True,
        switch_reliability_source="estimability",
        switch_estimability_min_minor_usage=min_minor_usage,
    )
    if consensus_runs and consensus_runs >= 2:
        kw["consensus_runs"] = consensus_runs
    if reliability:
        # NOTE (2026-06-24): estimability is now baked into the baseline kw above
        # (it was promoted to production), so this flag is redundant for the source
        # and only adds the '_reliability' method tag. Retained for back-reference;
        # the production trust-funnel re-run uses the plain 'isograph' baseline.
        kw["switch_reliability_weighting"] = True
        kw["switch_reliability_source"] = "estimability"
        kw["switch_estimability_min_minor_usage"] = min_minor_usage
    if reliability_floor > 0:
        # Cap the per-gene downweight at this floor (reliability in [floor, 1]) so a
        # noisy/degraded gene keeps a minimum switch contribution instead of being
        # fully pruned -- guards n_common against over-aggressive edge removal.
        kw["switch_reliability_floor"] = reliability_floor
    if max_module_frac is not None:
        # Giant-module cap (post-hoc split, the *negative* lever): any community over
        # this fraction of assigned genes is recursively re-clustered at escalating
        # resolution. Kept for back-reference; superseded by leiden_giant_frac.
        kw["max_module_frac"] = max_module_frac
    if leiden_giant_frac is not None:
        # Collapse fix C (resolution-sweep cap): select the smallest Leiden resolution
        # whose largest community is <= this fraction of genes, instead of splitting a
        # giant post-hoc. Distinct mechanism from max_module_frac; the sweep grid is
        # seeded at leiden_resolution (2.0 -> 2,4,8,16,32,64). Tags partitions _gcapNN.
        kw["leiden_max_giant_frac"] = leiden_giant_frac
    if grad_clip_norm is not None:
        # B.2 single-LR validation: global grad-norm clip (the S1 lever) is what lets one
        # fixed LR train every region without the per-cohort lr=3e-4 babysitting.
        kw["grad_clip_norm"] = grad_clip_norm
    if lr_override is not None:
        # B.2: force ONE learning rate across all regions, ignoring the per-cohort spec.
        kw["lr"] = lr_override
    elif spec["lr"] is not None:
        kw["lr"] = spec["lr"]
    return VaeModelConfig(**kw)


def _write_partition(modules: pd.DataFrame, cohort, region, method, seed, half) -> None:
    df = modules[["gene_id", "module_id"]].copy()
    df["gene_id"] = df["gene_id"].astype(str)
    df["method"], df["cohort"], df["region"] = method, cohort, region
    df["seed"], df["half"] = int(seed), half
    fname = f"{method}__{cohort}__{region}__seed{seed}__{half}.parquet"
    df.to_parquet(_partitions_dir() / fname, index=False, compression="zstd")


def fit_isograph(cohort: str, region: str, seeds: int, only_seed: int | None = None,
                 only_half: str | None = None, consensus_runs: int = 1,
                 reliability: bool = False, min_minor_usage: float = 0.1,
                 reliability_floor: float = 0.0,
                 max_module_frac: float | None = None,
                 leiden_giant_frac: float | None = None) -> None:
    """Fit IsoGraph on both split-halves for every seed in ``range(seeds)``, or for a
    single ``only_seed`` / ``only_half`` when given.

    With ``consensus_runs >= 2`` the VAE's Leiden step is stabilised by consensus over
    that many seeded, edge-weighted runs (IsoGraph ``consensus_runs`` config); partitions
    are tagged ``isograph_consensus`` so ``aggregate`` reports them as a separate method
    for a clean A/B against the baseline ``isograph`` partitions.

    Each VAE fit on the full transcript matrix has a large memory footprint (~18k genes
    -> a gene-gene similarity matrix held in several transient copies), and the
    allocations are not fully reclaimed within one Python process, so the running total
    grows fit-by-fit and OOM-kills a multi-fit loop. The SLURM driver therefore invokes
    one process per fit (``--seed k --half A``); peak is then a single fit, and process
    exit reclaims everything. Use ``--seed`` alone (both halves) for the smaller regions.
    """
    if cohort not in COHORTS:
        raise SystemExit(f"unknown cohort {cohort!r} (expected one of {list(COHORTS)})")
    spec = COHORTS[cohort]
    method = "isograph"
    if consensus_runs and consensus_runs >= 2:
        method += "_consensus"
    if reliability:
        method += "_reliability"
    if reliability_floor > 0:
        method += f"_f{int(round(reliability_floor * 100)):02d}"
    if max_module_frac is not None:
        method += f"_cap{int(round(max_module_frac * 100)):02d}"
    if leiden_giant_frac is not None:
        method += f"_gcap{int(round(leiden_giant_frac * 100)):02d}"
    bundle = load_dataset_bundle(rel(*spec["bundle_root"], region))
    sample_table = bundle.sample_table.reset_index(drop=True)
    tc = bundle.matrices["transcript_counts"]
    tt = bundle.feature_tables["transcript"]
    if spec["filter_transcripts"]:
        tc, tt = _filter_expressed_transcripts(tc, tt)
    else:
        # ensure a standalone (non-view) array we can column-slice cheaply
        tc = np.asarray(tc)
    del bundle
    n = sample_table.shape[0]

    cfg = _vae_config(spec, consensus_runs=consensus_runs,
                      reliability=reliability, min_minor_usage=min_minor_usage,
                      reliability_floor=reliability_floor,
                      max_module_frac=max_module_frac,
                      leiden_giant_frac=leiden_giant_frac)
    ks = [only_seed] if only_seed is not None else list(range(seeds))
    print(f"[{cohort}/{region}] {n} samples, {tc.shape[0]} transcripts | "
          f"method={method} | seeds {ks}", flush=True)

    for k in ks:
        a_idx, b_idx = _split_indices(n, SEED_BASE + k)
        halves = (("A", a_idx), ("B", b_idx))
        if only_half is not None:
            halves = tuple(h for h in halves if h[0] == only_half)
        for half, idx in halves:
            out = _partitions_dir() / f"{method}__{cohort}__{region}__seed{k}__{half}.parquet"
            if out.exists():
                print(f"  seed{k} {half}: exists, skipping", flush=True)
                continue
            t0 = time.time()
            try:
                st = sample_table.iloc[idx].reset_index(drop=True)
                art = VaeNetworkModel(cfg).fit(
                    transcript_counts=tc[:, idx], transcript_table=tt, sample_table=st,
                )
                mods = art.module_table[["gene_id", "module_id"]]
                _write_partition(mods, cohort, region, method, k, half)
                print(f"  seed{k} {half}: {mods['module_id'].nunique()} modules, "
                      f"{len(idx)} samples, {time.time() - t0:.0f}s", flush=True)
                del art, mods
                gc.collect()
            except Exception as e:  # one bad split should not lose the whole region
                print(f"  seed{k} {half}: FAILED ({type(e).__name__}: {e})", flush=True)


# ---------------------------------------------------------------------------
# Aggregation: ARI/NMI per (method, cohort, region, seed) + within-vs-cross summary
# ---------------------------------------------------------------------------
def _agreement(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    ca = a.set_index("gene_id")["module_id"]
    cb = b.set_index("gene_id")["module_id"]
    common = sorted(set(ca.index) & set(cb.index))
    if len(common) < 2:
        return {"n_common": len(common), "ari": np.nan, "nmi": np.nan}
    la, lb = ca.loc[common].values, cb.loc[common].values
    return {
        "n_common": len(common),
        "ari": float(adjusted_rand_score(la, lb)),
        "nmi": float(normalized_mutual_info_score(la, lb)),
        "n_mod_a": int(a["module_id"].nunique()),
        "n_mod_b": int(b["module_id"].nunique()),
    }


def _load_partitions() -> pd.DataFrame:
    pdir = _partitions_dir()
    files = [p for p in pdir.iterdir() if p.suffix == ".parquet"]
    if not files:
        raise SystemExit(f"no partition files in {pdir} — run fit-isograph / stability_wgcna.R first")
    return files


def _cross_cohort_rows() -> list[dict]:
    """Cross-cohort ARI/NMI from the saved full-data fits, for side-by-side contrast."""
    name = {"isograph": "isograph_vae", "wgcna": "wgcna_gene"}
    rows = []
    for method, mdir in name.items():
        for bs, gt, label in REGION_PAIRS:
            bp = rel("real_data", "brainseq", bs, "_m", mdir, "modules.parquet")
            gp = rel("real_data", "gtex", gt, "_m", mdir, "modules.parquet")
            if not (bp.exists() and gp.exists()):
                continue
            a = pd.read_parquet(bp)[["gene_id", "module_id"]]; a["gene_id"] = a["gene_id"].astype(str)
            b = pd.read_parquet(gp)[["gene_id", "module_id"]]; b["gene_id"] = b["gene_id"].astype(str)
            ag = _agreement(a, b)
            rows.append({"method": method, "region": label, "comparison": "cross_cohort",
                         "ari": ag["ari"], "nmi": ag["nmi"], "n_common": ag["n_common"]})
    return rows


def aggregate() -> None:
    out_dir = ensure_dir(_partitions_dir().parent)
    files = _load_partitions()
    parts = {}  # (method, cohort, region, seed) -> {half: df}
    for f in files:
        df = pd.read_parquet(f)
        key = (df["method"].iloc[0], df["cohort"].iloc[0], df["region"].iloc[0], int(df["seed"].iloc[0]))
        parts.setdefault(key, {})[df["half"].iloc[0]] = df

    pair_rows = []
    for (method, cohort, region, seed), halves in sorted(parts.items()):
        if "A" not in halves or "B" not in halves:
            print(f"  incomplete pair {method}/{cohort}/{region} seed{seed}: halves {list(halves)} — skipping")
            continue
        ag = _agreement(halves["A"], halves["B"])
        pair_rows.append({"method": method, "cohort": cohort, "region": region,
                          "seed": seed, **ag})
    pairs = pd.DataFrame(pair_rows)
    pairs.to_parquet(out_dir / "stability_pairs.parquet", index=False, compression="zstd")

    # within-cohort summary: mean +/- SD over seeds
    summ = (pairs.groupby(["method", "cohort", "region"])
                  .agg(n_seeds=("seed", "nunique"),
                       mean_ari=("ari", "mean"), sd_ari=("ari", "std"),
                       mean_nmi=("nmi", "mean"), sd_nmi=("nmi", "std"),
                       mean_n_common=("n_common", "mean"))
                  .reset_index())
    summ["comparison"] = "within_cohort_split_half"

    cross = pd.DataFrame(_cross_cohort_rows())
    combined = {
        "within_cohort": summ.to_dict(orient="records"),
        "cross_cohort": cross.to_dict(orient="records") if not cross.empty else [],
    }
    summ.to_parquet(out_dir / "stability_summary.parquet", index=False, compression="zstd")
    (out_dir / "stability_summary.json").write_text(json.dumps(combined, indent=2, default=str))

    print("\n=== within-cohort split-half stability (mean +/- SD over seeds) ===")
    print(summ.to_string(index=False))
    if not cross.empty:
        print("\n=== cross-cohort (full-data fits, for contrast) ===")
        print(cross.to_string(index=False))
    print(f"\nWrote {out_dir}/stability_summary.parquet + .json and stability_pairs.parquet")


# ---------------------------------------------------------------------------
# B.2 single-LR validation: ONE full-data fit per region at a single fixed LR
# ---------------------------------------------------------------------------
# Acceptance band for reconstruction RMSE. The HARD gate is no-divergence/no-OOM (that is
# what grad_clip_norm + the divergence guard actually fix); the band is a soft "quiet
# degradation" guard so a single LR that silently mis-fits one region is still caught.
# Calibrated from the B.2 single-LR run (2026-06-22): all 7 regions land in 0.58-0.81
# (BrainSEQ 0.58-0.62, GTEx 0.77-0.81) at one fixed lr=1e-3 + grad_clip=1.0 with no
# divergence — the original pre-registration guess of (1.00, 1.12) was simply too high.
# The band brackets the observed healthy range with margin.
RMSE_BAND = (0.45, 1.00)


def _rmse_dir():
    override = os.environ.get("STABILITY_RMSE_DIR")
    if override:
        return ensure_dir(Path(override))
    return ensure_dir(rel("real_data", "stability", "_m", "lr_validation"))


def fit_rmse(cohort: str, region: str, lr: float = 1e-3,
             grad_clip_norm: float | None = 1.0) -> None:
    """Validation gate B.2: ONE full-data IsoGraph fit at a single, fixed learning rate
    (no per-region tuning) with gradient clipping, recording reconstruction RMSE and
    whether training diverged.

    The gate: a single documented LR/optimizer config must train all 6 trust-funnel
    regions AND the GTEx region that diverged at the BrainSEQ default lr=1e-3
    (``nucleus_accumbens``) without divergence and with no OOM. This re-tests whether the
    merged ``grad_clip_norm`` + divergence guard remove the need for per-dataset LR
    babysitting (the production ``COHORTS`` spec still hard-codes GTEx ``lr=3e-4``). One
    JSON row per region is written to ``_rmse_dir()`` for ``aggregate-rmse`` to tabulate.
    """
    if cohort not in COHORTS:
        raise SystemExit(f"unknown cohort {cohort!r} (expected one of {list(COHORTS)})")
    spec = COHORTS[cohort]
    bundle = load_dataset_bundle(rel(*spec["bundle_root"], region))
    sample_table = bundle.sample_table.reset_index(drop=True)
    tc = bundle.matrices["transcript_counts"]
    tt = bundle.feature_tables["transcript"]
    if spec["filter_transcripts"]:
        tc, tt = _filter_expressed_transcripts(tc, tt)
    else:
        tc = np.asarray(tc)
    del bundle
    n = sample_table.shape[0]

    cfg = _vae_config(spec, lr_override=lr, grad_clip_norm=grad_clip_norm)
    print(f"[{cohort}/{region}] {n} samples, {tc.shape[0]} transcripts | "
          f"lr={lr} grad_clip_norm={grad_clip_norm}", flush=True)
    rec = {"cohort": cohort, "region": region, "n_samples": int(n),
           "n_transcripts": int(tc.shape[0]), "lr": float(lr),
           "grad_clip_norm": (float(grad_clip_norm) if grad_clip_norm is not None else None)}
    t0 = time.time()
    try:
        art = VaeNetworkModel(cfg).fit(
            transcript_counts=tc, transcript_table=tt, sample_table=sample_table)
        cal = art.calibration or {}
        rmse = cal.get("reconstruction_rmse")
        elbo = cal.get("vae_final_elbo")
        rmse_ok = rmse is not None and math.isfinite(float(rmse))
        rec.update({
            "reconstruction_rmse": (float(rmse) if rmse_ok else None),
            "vae_final_elbo": (float(elbo) if elbo is not None and math.isfinite(float(elbo)) else None),
            "vae_best_epoch": cal.get("vae_best_epoch"),
            "vae_early_stopped": cal.get("vae_early_stopped"),
            "vae_n_epochs_trained": cal.get("vae_n_epochs_trained"),
            "n_modules": int(art.module_table["module_id"].nunique()),
            "diverged": not rmse_ok,
            "seconds": round(time.time() - t0, 1),
            "status": "ok",
        })
        print(f"  rmse={rec['reconstruction_rmse']} elbo={rec['vae_final_elbo']} "
              f"early_stopped={rec['vae_early_stopped']} {rec['seconds']:.0f}s", flush=True)
        del art
        gc.collect()
    except Exception as e:  # OOM / runtime error is itself a gate failure for this LR
        rec.update({"status": f"FAILED:{type(e).__name__}", "error": str(e),
                    "diverged": True, "reconstruction_rmse": None,
                    "seconds": round(time.time() - t0, 1)})
        print(f"  FAILED ({type(e).__name__}: {e})", flush=True)
    out = _rmse_dir() / f"lrval__{cohort}__{region}.json"
    out.write_text(json.dumps(rec, indent=2, default=str))
    print(f"  wrote {out}", flush=True)


def aggregate_rmse() -> None:
    rdir = _rmse_dir()
    files = sorted(p for p in rdir.iterdir()
                   if p.suffix == ".json" and p.name.startswith("lrval__"))
    if not files:
        raise SystemExit(f"no lr-validation rows in {rdir} — run fit-rmse first")
    rows = [json.loads(p.read_text()) for p in files]
    df = pd.DataFrame(rows)
    lo, hi = RMSE_BAND

    def _verdict(r) -> str:
        if r.get("status") != "ok" or r.get("diverged"):
            return "DIVERGED/FAILED"
        rm = r.get("reconstruction_rmse")
        if rm is None or not (lo <= float(rm) <= hi):
            return "OUT-OF-BAND"
        return "PASS"

    df["verdict"] = df.apply(_verdict, axis=1)
    gate_pass = bool((df["verdict"] == "PASS").all())
    lr_used = float(df["lr"].iloc[0]) if "lr" in df.columns and len(df) else None

    out = _rmse_dir() / "lr_validation_summary"
    df.to_parquet(out.with_suffix(".parquet"), index=False, compression="zstd")
    summary = {"single_lr_gate_pass": gate_pass, "lr": lr_used,
               "rmse_band": list(RMSE_BAND), "rows": df.to_dict(orient="records")}
    out.with_suffix(".json").write_text(json.dumps(summary, indent=2, default=str))

    show = [c for c in ["cohort", "region", "n_samples", "reconstruction_rmse",
                        "vae_final_elbo", "vae_best_epoch", "vae_early_stopped",
                        "n_modules", "diverged", "verdict", "seconds", "status"]
            if c in df.columns]
    print("\n=== single-LR full-data RMSE validation (gate B.2) ===")
    print(df[show].sort_values(["cohort", "region"]).to_string(index=False))
    print(f"\nSINGLE-LR GATE: {'PASS' if gate_pass else 'FAIL'}  "
          f"(lr={lr_used}; all regions in RMSE band {RMSE_BAND}, no divergence/OOM)")
    print(f"Wrote {out}.parquet + .json")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    fi = sub.add_parser("fit-isograph", help="fit IsoGraph on both halves per seed")
    fi.add_argument("--cohort", required=True, choices=list(COHORTS))
    fi.add_argument("--region", required=True)
    fi.add_argument("--seeds", type=int, default=5,
                    help="number of seeds (0..seeds-1) when --seed is not given")
    fi.add_argument("--seed", type=int, default=None,
                    help="fit only this single seed (one process per seed avoids OOM)")
    fi.add_argument("--half", choices=["A", "B"], default=None,
                    help="fit only this half (with --seed: one process per fit, lowest peak)")
    fi.add_argument("--consensus", type=int, default=1, metavar="N",
                    help="consensus Leiden over N seeded runs (>=2 enables; tags "
                         "partitions 'isograph_consensus' for A/B vs baseline)")
    fi.add_argument("--reliability", action="store_true",
                    help="covariate-free isoform-estimability downweighting of "
                         "switch-switch edges (tags partitions 'isograph_reliability')")
    fi.add_argument("--min-minor-usage", type=float, default=0.1, metavar="U",
                    help="minor-isoform usage floor for --reliability (default 0.1)")
    fi.add_argument("--reliability-floor", type=float, default=0.0, metavar="F",
                    help="floor for the per-gene reliability weight (reliability in "
                         "[F,1]); caps downweighting. Tags partitions '_fNN'.")
    fi.add_argument("--max-module-frac", type=float, default=None, metavar="C",
                    help="giant-module cap: recursively re-cluster any module over this "
                         "fraction of assigned genes (e.g. 0.15). Tags partitions '_capNN'.")
    fi.add_argument("--leiden-giant-frac", type=float, default=None, metavar="G",
                    help="collapse fix C: select the smallest Leiden resolution whose "
                         "largest module is <= this fraction of genes (e.g. 0.15), instead "
                         "of post-hoc splitting. Tags partitions '_gcapNN'.")
    sub.add_parser("aggregate", help="compute ARI/NMI over all partitions and summarize")

    fr = sub.add_parser("fit-rmse",
                        help="B.2: one full-data fit at a single LR; record reconstruction RMSE")
    fr.add_argument("--cohort", required=True, choices=list(COHORTS))
    fr.add_argument("--region", required=True)
    fr.add_argument("--lr", type=float, default=1e-3, metavar="LR",
                    help="single learning rate applied to every region (default 1e-3 = "
                         "the BrainSEQ default that diverges on some GTEx without clipping)")
    fr.add_argument("--grad-clip-norm", type=float, default=1.0, metavar="C",
                    help="global grad-norm clip (default 1.0; <=0 disables)")
    sub.add_parser("aggregate-rmse",
                   help="tabulate lr-validation rows -> single-LR gate (PASS/FAIL)")
    args = ap.parse_args()

    if args.cmd == "fit-isograph":
        fit_isograph(args.cohort, args.region, args.seeds,
                     only_seed=args.seed, only_half=args.half,
                     consensus_runs=args.consensus,
                     reliability=args.reliability,
                     min_minor_usage=args.min_minor_usage,
                     reliability_floor=args.reliability_floor,
                     max_module_frac=args.max_module_frac,
                     leiden_giant_frac=args.leiden_giant_frac)
    elif args.cmd == "fit-rmse":
        gc_norm = args.grad_clip_norm if args.grad_clip_norm and args.grad_clip_norm > 0 else None
        fit_rmse(args.cohort, args.region, lr=args.lr, grad_clip_norm=gc_norm)
    elif args.cmd == "aggregate-rmse":
        aggregate_rmse()
    else:
        aggregate()


if __name__ == "__main__":
    main()
