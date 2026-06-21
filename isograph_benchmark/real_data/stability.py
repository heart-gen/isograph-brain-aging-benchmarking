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
import time

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
# Keep these in sync with run_models if the production configs change.
COHORTS = {
    "brainseq": {
        "bundle_root": ("inputs", "bundles", "brainseq_v1"),
        "regions": ["caudate", "hippocampus", "dlpfc"],
        "covariates": ["Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
                       "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5"],
        "age_col": "Age",
        "filter_transcripts": True,   # run_brainseq_region filters; GTEx bundles are pre-filtered
        "lr": None,                   # default (1e-3)
    },
    "gtex": {
        "bundle_root": ("inputs", "bundles", "gtex_v11_brain"),
        "regions": ["caudate_basal_ganglia", "hippocampus", "frontal_cortex_ba9"],
        "covariates": ["SEX", "SMRIN", "SMTSISCH", "SMMAPRT"],
        "age_col": "AGE",
        "filter_transcripts": False,
        "lr": 3e-4,                   # GTEx needs lr=3e-4 (default diverges at this scale)
    },
}


def _partitions_dir():
    return ensure_dir(rel("real_data", "stability", "_m", "partitions"))


def _split_indices(n: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Random 50/50 partition of sample positions (half A gets floor(n/2))."""
    perm = np.random.default_rng(seed).permutation(n)
    h = n // 2
    return perm[:h], perm[h:]


def _vae_config(spec: dict, consensus_runs: int = 1, reliability: bool = False,
                min_minor_usage: float = 0.1, tin: bool = False,
                reliability_floor: float = 0.0,
                extra_covariates: list[str] | None = None,
                max_module_frac: float | None = None) -> VaeModelConfig:
    covariates = list(spec["covariates"]) + list(extra_covariates or [])
    kw = dict(
        hidden_dim=256, latent_dim=32, n_epochs=500,
        residualize_covariates=covariates,
        min_module_size=20, trait_columns=[spec["age_col"]], random_state=VAE_SEED,
        allow_abundance_abundance=False, alpha_switch=0.5, leiden_resolution=2.0,
    )
    if consensus_runs and consensus_runs >= 2:
        kw["consensus_runs"] = consensus_runs
    if reliability:
        # Covariate-free isoform-estimability downweighting of switch-switch edges:
        # genes whose minor isoform lacks read support carry a noise switch
        # coordinate that flips between split-halves; downweighting their edges
        # should stabilise the surviving module structure.
        kw["switch_reliability_weighting"] = True
        kw["switch_reliability_source"] = "estimability"
        kw["switch_estimability_min_minor_usage"] = min_minor_usage
    if tin:
        # Per-gene differential transcript-integrity (TIN) downweighting: genes whose
        # isoform composition tracks within-gene differential degradation are switch
        # artifacts -> downweight their switch-switch edges (needs transcript_tin).
        kw["switch_reliability_weighting"] = True
        kw["switch_reliability_source"] = "tin_differential"
    if reliability_floor > 0:
        # Cap the per-gene downweight at this floor (reliability in [floor, 1]) so a
        # noisy/degraded gene keeps a minimum switch contribution instead of being
        # fully pruned -- guards n_common against over-aggressive edge removal.
        kw["switch_reliability_floor"] = reliability_floor
    if max_module_frac is not None:
        # Giant-module cap: any community over this fraction of assigned genes is
        # recursively re-clustered at escalating resolution (IsoGraph core).
        kw["max_module_frac"] = max_module_frac
    if spec["lr"] is not None:
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
                 tin: bool = False, median_tin_covariate: bool = False,
                 reliability_floor: float = 0.0,
                 max_module_frac: float | None = None) -> None:
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
    if tin:
        method += "_tin"
    if median_tin_covariate:
        method += "_mediantin"
    if reliability_floor > 0:
        method += f"_f{int(round(reliability_floor * 100)):02d}"
    if max_module_frac is not None:
        method += f"_cap{int(round(max_module_frac * 100)):02d}"
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

    # Optional TIN inputs (brainseq caudate pilot). tin_mat is aligned row-for-row to
    # the (filtered) transcript table and column-for-column to sample_table order, so
    # it can be sliced by the same half index as tc.
    tin_mat = None
    extra_covariates: list[str] = []
    if tin or median_tin_covariate:
        from isograph_benchmark.real_data.tin import load_tin_aligned, load_sample_median_tin
        sample_ids = sample_table["sample_id"].astype(str).tolist()
        if tin:
            tin_mat = load_tin_aligned(cohort, region,
                                       tt["transcript_id"].astype(str).tolist(), sample_ids)
        if median_tin_covariate:
            sample_table = sample_table.copy()
            sample_table["median_tin"] = load_sample_median_tin(cohort, region, sample_ids)
            extra_covariates.append("median_tin")

    cfg = _vae_config(spec, consensus_runs=consensus_runs,
                      reliability=reliability, min_minor_usage=min_minor_usage,
                      tin=tin, reliability_floor=reliability_floor,
                      extra_covariates=extra_covariates,
                      max_module_frac=max_module_frac)
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
                tin_half = tin_mat[:, idx] if tin_mat is not None else None
                art = VaeNetworkModel(cfg).fit(
                    transcript_counts=tc[:, idx], transcript_table=tt, sample_table=st,
                    transcript_tin=tin_half,
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
    out_dir = ensure_dir(rel("real_data", "stability", "_m"))
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
    fi.add_argument("--tin", action="store_true",
                    help="per-gene differential-TIN switch-edge downweighting "
                         "(needs cached TIN; tags partitions 'isograph_tin')")
    fi.add_argument("--median-tin-covariate", action="store_true",
                    help="add per-sample median TIN to residualization covariates "
                         "(tags partitions 'isograph_mediantin')")
    fi.add_argument("--reliability-floor", type=float, default=0.0, metavar="F",
                    help="floor for the per-gene reliability weight (reliability in "
                         "[F,1]); caps downweighting. Tags partitions '_fNN'.")
    fi.add_argument("--max-module-frac", type=float, default=None, metavar="C",
                    help="giant-module cap: recursively re-cluster any module over this "
                         "fraction of assigned genes (e.g. 0.15). Tags partitions '_capNN'.")
    sub.add_parser("aggregate", help="compute ARI/NMI over all partitions and summarize")
    args = ap.parse_args()

    if args.cmd == "fit-isograph":
        fit_isograph(args.cohort, args.region, args.seeds,
                     only_seed=args.seed, only_half=args.half,
                     consensus_runs=args.consensus,
                     reliability=args.reliability,
                     min_minor_usage=args.min_minor_usage,
                     tin=args.tin,
                     median_tin_covariate=args.median_tin_covariate,
                     reliability_floor=args.reliability_floor,
                     max_module_frac=args.max_module_frac)
    else:
        aggregate()


if __name__ == "__main__":
    main()
