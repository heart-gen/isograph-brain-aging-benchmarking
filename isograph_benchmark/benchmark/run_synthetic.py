from __future__ import annotations

import hashlib
from itertools import product

import pandas as pd

from isograph_benchmark.config import load_yaml
from isograph_benchmark.paths import ensure_dir, rel


RESOURCE_DEFAULTS = {
    "cpu_short": {"requested_cpus": 16, "requested_gpus": 0, "mem_per_cpu_gb": 2, "requested_mem_gb": 32, "target_minutes": 90},
    "vae": {"requested_cpus": 32, "requested_gpus": 0, "mem_per_cpu_gb": 2, "requested_mem_gb": 64, "target_minutes": 90},
    "wgcna_cpu": {"requested_cpus": 50, "requested_gpus": 0, "mem_per_cpu_gb": 2, "requested_mem_gb": 100, "target_minutes": 90},
    "scale": {"requested_cpus": 64, "requested_gpus": 0, "mem_per_cpu_gb": 2, "requested_mem_gb": 128, "target_minutes": 240},
    "gpu": {"requested_cpus": 8, "requested_gpus": 1, "mem_per_cpu_gb": 8, "requested_mem_gb": 64, "target_minutes": 90},
    "gpu_scale": {"requested_cpus": 8, "requested_gpus": 1, "mem_per_cpu_gb": 16, "requested_mem_gb": 128, "target_minutes": 240},
}

# Scenarios that use scale resources and scale_methods.
# scale_realistic is kept separate from scale so existing run hashes are
# preserved while adding the BrainSEQ-scale (16k genes / 300 samples) point.
_SCALE_SCENARIOS = frozenset({"scale", "scale_realistic"})


def resource_class(method: str, scenario: str) -> str:
    if method in ("isograph_gpu_latent", "isograph_vae_gpu"):
        return "gpu_scale" if scenario in _SCALE_SCENARIOS else "gpu"
    if scenario in _SCALE_SCENARIOS:
        return "scale"
    if method == "wgcna_gene":
        return "wgcna_cpu"
    if method in ("isograph_vae", "isograph_vae_multiplex", "isograph_cpu_latent"):
        return "vae" if method in ("isograph_vae", "isograph_vae_multiplex") else "cpu_short"
    return "cpu_short"


def stable_id(values: dict[str, object], keys: list[str]) -> str:
    payload = "|".join(f"{key}={values.get(key, '')}" for key in keys)
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:16]


_SCENARIO_DATASET_KEYS = [
    "scenario",
    "n_genes",
    "n_samples",
    "switching_fraction",
    "noise_sd",
    "abundance_imbalance",
    "count_dispersion",
    "interaction_strength",
    "interaction_fraction",
    "seed",
]

# abundance_switch_mixed datasets additionally hash on abundance_fraction so that
# datasets with the same base params but different abundance_fraction get distinct IDs.
# Existing scenario hashes are intentionally unchanged (backward compatible).
_MULTIPLEX_SCENARIO_DATASET_KEYS = _SCENARIO_DATASET_KEYS + ["abundance_fraction"]


def _scenario_methods(cfg: dict, scenario: str) -> list[str]:
    if scenario in _SCALE_SCENARIOS:
        return cfg["scale_methods"]
    if scenario == "abundance_switch_mixed":
        return cfg.get("multiplex_methods", cfg["methods"])
    return cfg["methods"]


def expand_grid() -> pd.DataFrame:
    cfg = load_yaml("configs/synthetic_grid.yaml")
    all_known_methods = list(dict.fromkeys(cfg["methods"] + cfg.get("multiplex_methods", [])))
    rows: list[dict[str, object]] = []
    for scenario, params in cfg["scenarios"].items():
        scenario_seed_count = params.get("seed_count")
        grid_params = {key: value for key, value in params.items() if key != "seed_count"}
        keys = list(grid_params)
        for values in product(*[grid_params[key] for key in keys]):
            seed_count = int(
                scenario_seed_count
                if scenario_seed_count is not None
                else (cfg["seed_count_scale"] if scenario in _SCALE_SCENARIOS else cfg["seed_count_core"])
            )
            parameter_values = dict(zip(keys, values, strict=True))
            methods = _scenario_methods(cfg, scenario)
            for seed_idx in range(seed_count):
                dataset_seed = int(cfg["base_seed"]) + seed_idx
                dataset = {
                    **parameter_values,
                    "scenario": scenario,
                    "seed": dataset_seed,
                    "replicate": seed_idx,
                }
                hash_keys = _MULTIPLEX_SCENARIO_DATASET_KEYS if scenario == "abundance_switch_mixed" else _SCENARIO_DATASET_KEYS
                dataset_id = stable_id(dataset, hash_keys)
                for method in all_known_methods:
                    if method not in methods:
                        continue
                    row = dict(dataset)
                    rc = resource_class(method, scenario)
                    defaults = RESOURCE_DEFAULTS[rc]
                    row.update(
                        {
                            "method": method,
                            "dataset_id": dataset_id,
                            "run_id": stable_id({**dataset, "method": method}, list(dataset) + ["method"]),
                            "resource_class": rc,
                            **defaults,
                        }
                    )
                    rows.append(row)
    grid = pd.DataFrame(rows)
    first_cols = [
        "run_id",
        "dataset_id",
        "scenario",
        "method",
        "resource_class",
        "seed",
        "replicate",
        "n_genes",
        "n_samples",
        "switching_fraction",
        "noise_sd",
    ]
    ordered = [col for col in first_cols if col in grid.columns]
    ordered.extend(col for col in grid.columns if col not in ordered)
    return grid.loc[:, ordered].sort_values(["scenario", "replicate", "dataset_id", "method"]).reset_index(drop=True)


def main() -> None:
    out = rel("benchmark", "00_design", "_m", "synthetic_run_grid.parquet")
    ensure_dir(out.parent)
    grid = expand_grid()
    grid.to_parquet(out, index=False, compression="zstd")
    print(f"Wrote {len(grid):,} synthetic run definitions to {out}")


if __name__ == "__main__":
    main()
