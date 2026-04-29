from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from isograph.io.artifacts import (
    DatasetBundle,
    build_feature_spec,
    build_matrix_spec,
    save_dataset_bundle,
)
from isograph.validation import DatasetManifest


def _as_float(row: pd.Series, key: str, default: float) -> float:
    value = row.get(key, default)
    if pd.isna(value):
        return default
    return float(value)


def _as_int(row: pd.Series, key: str, default: int) -> int:
    value = row.get(key, default)
    if pd.isna(value):
        return default
    return int(value)


def dataset_dir(row: pd.Series, root: Path) -> Path:
    return root / str(row["dataset_id"])


def build_synthetic_bundle(row: pd.Series) -> DatasetBundle:
    rng = np.random.default_rng(_as_int(row, "seed", 13))
    n_genes = _as_int(row, "n_genes", 400)
    n_samples = _as_int(row, "n_samples", 160)
    switching_fraction = np.clip(_as_float(row, "switching_fraction", 0.25), 0.0, 1.0)
    noise_sd = max(_as_float(row, "noise_sd", 0.20), 0.0)
    abundance_imbalance = max(_as_float(row, "abundance_imbalance", 1.0), 1.0)
    count_dispersion = max(_as_float(row, "count_dispersion", 15.0), 1.0)
    interaction_strength = max(_as_float(row, "interaction_strength", 0.0), 0.0)
    interaction_fraction = np.clip(_as_float(row, "interaction_fraction", 0.0), 0.0, 1.0)

    n_modules = max(2, min(8, n_genes // 50))
    n_switching = max(1, int(round(n_genes * switching_fraction)))
    gene_ids = [f"G{i:05d}" for i in range(n_genes)]
    sample_ids = [f"S{i:04d}" for i in range(n_samples)]
    switching_mask = np.zeros(n_genes, dtype=bool)
    switching_mask[:n_switching] = True
    rng.shuffle(switching_mask)

    module_index = np.full(n_genes, -1, dtype=int)
    switching_modules = np.arange(n_switching) % n_modules
    rng.shuffle(switching_modules)
    module_index[switching_mask] = switching_modules

    dx = np.where(np.arange(n_samples) < n_samples / 2, "Control", "SCZD")
    age = np.linspace(25, 85, n_samples) + rng.normal(0, 3, n_samples)
    sex = np.where(np.arange(n_samples) % 2 == 0, "F", "M")
    sample_table = pd.DataFrame({"sample_id": sample_ids, "Dx": dx, "Age": age, "Sex": sex})

    module_latent = rng.normal(size=(n_modules, n_samples))
    module_latent += ((age - age.mean()) / age.std())[None, :] * rng.normal(0.15, 0.05, (n_modules, 1))
    module_latent += (dx == "SCZD")[None, :] * rng.normal(0.15, 0.05, (n_modules, 1))

    signal = rng.normal(0, noise_sd, size=(n_genes, n_samples))
    for gene_idx in np.where(switching_mask)[0]:
        signal[gene_idx] += module_latent[module_index[gene_idx]]

    interaction_genes = np.where(switching_mask)[0]
    if interaction_strength > 0 and len(interaction_genes) > 0:
        n_interaction = int(round(len(interaction_genes) * interaction_fraction))
        chosen = rng.choice(interaction_genes, size=n_interaction, replace=False)
        interaction_term = module_latent[0] * module_latent[1 % n_modules]
        interaction_term = (interaction_term - interaction_term.mean()) / interaction_term.std()
        signal[chosen] += interaction_strength * interaction_term

    baseline_logit = np.log(abundance_imbalance)
    p1 = 1.0 / (1.0 + np.exp(-(signal + baseline_logit)))
    p1[~switching_mask] = 1.0 / (1.0 + abundance_imbalance)
    p1 = np.clip(p1 + rng.normal(0, noise_sd * 0.15, p1.shape), 1e-4, 1 - 1e-4)

    gene_means = rng.lognormal(mean=np.log(80), sigma=0.8, size=(n_genes, 1))
    sample_depth = rng.lognormal(mean=0.0, sigma=0.25, size=(1, n_samples))
    expected_total = gene_means * sample_depth
    gamma_shape = count_dispersion
    gamma_scale = expected_total / gamma_shape
    total_rate = rng.gamma(shape=gamma_shape, scale=gamma_scale)
    totals = rng.poisson(np.maximum(total_rate, 1.0)).astype(float)

    transcript_counts = np.zeros((n_genes * 2, n_samples), dtype=float)
    transcript_rows: list[dict[str, object]] = []
    for gene_idx, gene_id in enumerate(gene_ids):
        tx1 = rng.binomial(totals[gene_idx].astype(int), p1[gene_idx]).astype(float)
        tx2 = np.maximum(totals[gene_idx] - tx1, 0.0)
        transcript_counts[gene_idx * 2] = tx1
        transcript_counts[gene_idx * 2 + 1] = tx2
        transcript_rows.extend(
            [
                {"transcript_id": f"{gene_id}_T1", "gene_id": gene_id, "length": 1000},
                {"transcript_id": f"{gene_id}_T2", "gene_id": gene_id, "length": 900},
            ]
        )

    gene_counts = transcript_counts.reshape(n_genes, 2, n_samples).sum(axis=1)
    psi = p1
    gene_table = pd.DataFrame({"gene_id": gene_ids})
    transcript_table = pd.DataFrame(transcript_rows)
    psi_table = pd.DataFrame({"psi_uid": [f"PSI_{gene_id}" for gene_id in gene_ids], "gene_id": gene_ids})
    truth_modules = pd.DataFrame(
        {
            "gene_id": [gene_ids[i] for i in np.where(switching_mask)[0]],
            "module_id": [f"M{module_index[i]:03d}" for i in np.where(switching_mask)[0]],
        }
    )
    truth_switch = pd.DataFrame({"gene_id": gene_ids, "has_switch": switching_mask})

    manifest = DatasetManifest(
        dataset_name=str(row["dataset_id"]),
        suite_name="isograph_brain_aging_synthetic",
        description=f"Synthetic IsoGraph benchmark dataset for {row['scenario']}",
        sample_table="samples.parquet",
        feature_tables=[
            build_feature_spec("gene", "genes.parquet", gene_table),
            build_feature_spec("transcript", "transcripts.parquet", transcript_table),
            build_feature_spec("psi", "psi.parquet", psi_table),
            build_feature_spec("truth_module", "truth_modules.parquet", truth_modules),
            build_feature_spec("truth_switch", "truth_switch.parquet", truth_switch),
        ],
        matrices=[
            build_matrix_spec("gene_counts", "gene_counts.npz", gene_counts),
            build_matrix_spec("transcript_counts", "transcript_counts.npz", transcript_counts),
            build_matrix_spec("psi", "psi.npz", psi),
        ],
        provenance={
            "generator": "isograph_brain_aging_benchmark_v1",
            "scenario": str(row["scenario"]),
            "seed": str(row["seed"]),
        },
        truth_tables=["truth_modules.parquet", "truth_switch.parquet"],
    )
    return DatasetBundle(
        manifest=manifest,
        sample_table=sample_table,
        feature_tables={
            "gene": gene_table,
            "transcript": transcript_table,
            "psi": psi_table,
            "truth_module": truth_modules,
            "truth_switch": truth_switch,
        },
        matrices={"gene_counts": gene_counts, "transcript_counts": transcript_counts, "psi": psi},
        truth_tables={"truth_modules.parquet": truth_modules, "truth_switch.parquet": truth_switch},
    )


def ensure_dataset(row: pd.Series, root: Path) -> Path:
    out = dataset_dir(row, root)
    if (out / "manifest.json").exists():
        return out
    save_dataset_bundle(build_synthetic_bundle(row), out)
    return out
