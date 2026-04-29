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
}


def resource_class(method: str, scenario: str) -> str:
    if scenario == "scale":
        return "scale"
    if method == "wgcna_gene":
        return "wgcna_cpu"
    if method == "isograph_vae":
        return "vae"
    return "cpu_short"


def stable_id(values: dict[str, object], keys: list[str]) -> str:
    payload = "|".join(f"{key}={values.get(key, '')}" for key in keys)
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:16]


def expand_grid() -> pd.DataFrame:
    cfg = load_yaml("configs/synthetic_grid.yaml")
    rows: list[dict[str, object]] = []
    for scenario, params in cfg["scenarios"].items():
        scenario_seed_count = params.get("seed_count")
        grid_params = {key: value for key, value in params.items() if key != "seed_count"}
        keys = list(grid_params)
        for values in product(*[grid_params[key] for key in keys]):
            seed_count = int(
                scenario_seed_count
                if scenario_seed_count is not None
                else (cfg["seed_count_scale"] if scenario == "scale" else cfg["seed_count_core"])
            )
            parameter_values = dict(zip(keys, values, strict=True))
            methods = cfg["scale_methods"] if scenario == "scale" else cfg["methods"]
            for seed_idx in range(seed_count):
                dataset_seed = int(cfg["base_seed"]) + seed_idx
                dataset = {
                    **parameter_values,
                    "scenario": scenario,
                    "seed": dataset_seed,
                    "replicate": seed_idx,
                }
                dataset_id = stable_id(
                    dataset,
                    [
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
                    ],
                )
                for method in cfg["methods"]:
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
