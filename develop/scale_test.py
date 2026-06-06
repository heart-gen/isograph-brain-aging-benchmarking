"""Fast IsoGraph scale-test for troubleshooting at real-data-comparable sizes.

Generates synthetic bundles at two scales (intermediate and large) and fits IsoGraph
with multiple configurations, reporting module quality and runtime. NOT part of the
publication benchmark — use this for quick parameter validation before running full
real-data or full-benchmark jobs.

Scales:
  intermediate:  5,000 genes / 250 samples
  large:        17,000 genes / 300 samples  (matches BrainSEQ)

Configs tested per scale:
  A: allow_abundance_abundance=False  (current default)
  B: allow_abundance_abundance=True + alpha_abundance_grid  (re-enable test)

Usage (from project root, with the isograph conda env active):
    python develop/scale_test.py                        # all scales, n_tx=2
    python develop/scale_test.py --scale intermediate   # only 5k-gene runs
    python develop/scale_test.py --n-tx 5               # variable isoforms
    python develop/scale_test.py --leiden-resolution 3.0
    python develop/scale_test.py --replicates 1         # single fast check

Results written to: develop/_m/scale_test_results.parquet
"""

from __future__ import annotations

import argparse
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── project root so relative imports work when called directly ───────────────
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from isograph.evaluation.metrics import module_recovery_score
from isograph.models.vae import VaeNetworkModel
from isograph.workflow.config import VaeModelConfig
from isograph_benchmark.benchmark.synthetic_data import build_synthetic_bundle
from isograph_benchmark.paths import ensure_dir

DEVELOP_DIR = Path(__file__).resolve().parent
DATA_DIR = DEVELOP_DIR / "_data"
RESULTS_DIR = DEVELOP_DIR / "_m"

SCALES = {
    "intermediate": {"n_genes": 5_000, "n_samples": 250},
    "large":        {"n_genes": 17_000, "n_samples": 300},
}

ALPHA_ABUNDANCE_GRID = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]


def _peak_rss_mb() -> float:
    """Peak resident set size in MB (Linux)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _build_config(
    allow_abundance: bool,
    leiden_resolution: float,
    n_genes: int,
) -> VaeModelConfig:
    hidden_dim = 128 if n_genes <= 1000 else 256
    base = dict(
        hidden_dim=hidden_dim,
        latent_dim=32,
        n_epochs=300,
        patience=35,
        min_module_size=20,
        random_state=13,
        alpha_switch=0.5,
        leiden_resolution=leiden_resolution,
    )
    if allow_abundance:
        return VaeModelConfig(
            **base,
            allow_abundance_abundance=True,
            alpha_abundance_grid=ALPHA_ABUNDANCE_GRID,
        )
    return VaeModelConfig(**base, allow_abundance_abundance=False)


def run_one(
    scale_name: str,
    n_tx: int,
    allow_abundance: bool,
    leiden_resolution: float,
    replicate: int,
    switching_fraction: float = 0.25,
    noise_sd: float = 0.15,
) -> dict:
    dims = SCALES[scale_name]
    n_genes = dims["n_genes"]
    n_samples = dims["n_samples"]
    config_label = "with_abundance" if allow_abundance else "no_abundance"

    row = pd.Series({
        "dataset_id": f"scale_{scale_name}_{config_label}_ntx{n_tx}_rep{replicate}",
        "scenario": "scale_test",
        "seed": replicate * 100 + 7,
        "n_genes": n_genes,
        "n_samples": n_samples,
        "n_transcripts_per_gene": n_tx,
        "switching_fraction": switching_fraction,
        "noise_sd": noise_sd,
        "abundance_fraction": 0.0,
    })

    print(
        f"\n[{scale_name}|ntx={n_tx}|{config_label}|rep={replicate}] "
        f"Generating bundle ({n_genes}g×{n_samples}s) ...",
        flush=True,
    )
    t_gen = time.time()
    bundle = build_synthetic_bundle(row)
    print(f"  generated in {time.time()-t_gen:.1f}s", flush=True)

    cfg = _build_config(allow_abundance, leiden_resolution, n_genes)

    print(f"  Fitting IsoGraph (leiden_resolution={leiden_resolution}) ...", flush=True)
    t_fit = time.time()
    artifacts = VaeNetworkModel(cfg).fit(
        transcript_counts=bundle.matrices["transcript_counts"],
        transcript_table=bundle.feature_tables["transcript"],
        sample_table=bundle.sample_table,
    )
    elapsed = time.time() - t_fit
    peak_mb = _peak_rss_mb()

    truth_modules = bundle.truth_tables.get("truth_modules.parquet", pd.DataFrame())
    recovery = float("nan")
    if not artifacts.module_table.empty and not truth_modules.empty:
        try:
            recovery = module_recovery_score(artifacts.module_table, truth_modules)
        except Exception as exc:
            print(f"  WARNING: module_recovery_score failed: {exc}")

    n_modules = 0
    giant_fraction = 0.0
    if not artifacts.module_table.empty:
        sizes = artifacts.module_table.groupby("module_id").size().sort_values(ascending=False)
        n_modules = int(len(sizes))
        n_assigned = len(artifacts.module_table)
        giant_fraction = float(sizes.iloc[0] / n_assigned) if n_assigned > 0 else 0.0

    calibration = artifacts.calibration or {}

    result = {
        "scale": scale_name,
        "n_genes": n_genes,
        "n_samples": n_samples,
        "n_transcripts_per_gene": n_tx,
        "config": config_label,
        "leiden_resolution": leiden_resolution,
        "replicate": replicate,
        "n_modules": n_modules,
        "giant_fraction": round(giant_fraction, 4),
        "module_recovery": round(recovery, 4) if not np.isnan(recovery) else float("nan"),
        "alpha_switch_selected": calibration.get("alpha_switch"),
        "alpha_abundance_selected": calibration.get("alpha_abundance"),
        "elapsed_fit_s": round(elapsed, 1),
        "peak_rss_mb": round(peak_mb, 1),
    }

    print(
        f"  DONE: n_modules={n_modules}, giant={giant_fraction:.1%}, "
        f"recovery={recovery:.3f}, elapsed={elapsed:.0f}s, peak={peak_mb:.0f}MB",
        flush=True,
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Scale test for IsoGraph parameter troubleshooting.")
    parser.add_argument(
        "--scale", choices=list(SCALES), default=None,
        help="Run only the named scale (default: both intermediate and large).",
    )
    parser.add_argument(
        "--n-tx", type=int, default=None,
        help="Fix n_transcripts_per_gene (default: test both 2 and 5).",
    )
    parser.add_argument(
        "--leiden-resolution", type=float, default=2.0,
        help="Leiden resolution for module detection (default: 2.0; update after sweep_leiden).",
    )
    parser.add_argument(
        "--replicates", type=int, default=3,
        help="Number of random replicates per config (default: 3).",
    )
    parser.add_argument(
        "--no-abundance", action="store_true",
        help="Skip the allow_abundance_abundance=True variant.",
    )
    args = parser.parse_args()

    scales = [args.scale] if args.scale else list(SCALES)
    n_tx_values = [args.n_tx] if args.n_tx else [2, 5]
    configs = [False] if args.no_abundance else [False, True]

    results = []
    for scale in scales:
        for n_tx in n_tx_values:
            for allow_abundance in configs:
                for rep in range(1, args.replicates + 1):
                    try:
                        row = run_one(
                            scale_name=scale,
                            n_tx=n_tx,
                            allow_abundance=allow_abundance,
                            leiden_resolution=args.leiden_resolution,
                            replicate=rep,
                        )
                        results.append(row)
                    except Exception as exc:
                        print(f"ERROR in {scale}/ntx={n_tx}/abundance={allow_abundance}/rep={rep}: {exc}")
                        results.append({
                            "scale": scale, "n_genes": SCALES[scale]["n_genes"],
                            "n_samples": SCALES[scale]["n_samples"],
                            "n_transcripts_per_gene": n_tx,
                            "config": "with_abundance" if allow_abundance else "no_abundance",
                            "leiden_resolution": args.leiden_resolution,
                            "replicate": rep, "error": str(exc),
                        })

    df = pd.DataFrame(results)
    print("\n=== SCALE TEST SUMMARY ===")
    cols = ["scale", "n_transcripts_per_gene", "config", "leiden_resolution",
            "replicate", "n_modules", "giant_fraction", "module_recovery", "elapsed_fit_s", "peak_rss_mb"]
    available = [c for c in cols if c in df.columns]
    print(df[available].to_string(index=False))

    ensure_dir(RESULTS_DIR)
    out_path = RESULTS_DIR / "scale_test_results.parquet"
    df.to_parquet(out_path, index=False, compression="zstd")
    print(f"\nResults written to {out_path}")


if __name__ == "__main__":
    main()
