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
    n_transcripts_per_gene = _as_int(row, "n_transcripts_per_gene", 2)
    switching_fraction = np.clip(_as_float(row, "switching_fraction", 0.25), 0.0, 1.0)
    # fraction of module genes driven by abundance rather than switching [0,1]
    abundance_fraction = np.clip(_as_float(row, "abundance_fraction", 0.0), 0.0, 1.0)
    # Dual-signal module genes: a fraction of switch-driven module genes that ALSO
    # carry abundance signal for the same module (realistic co-regulation, where a
    # gene both switches isoforms and shifts total expression). Under 3' degradation
    # the switch channel is corrupted but the gene's own abundance keeps it in the
    # module -- the basis of the multiplex abundance fallback. Default 0.0 draws no
    # extra RNG, leaving existing scenarios byte-identical.
    dual_signal_fraction = np.clip(_as_float(row, "dual_signal_fraction", 0.0), 0.0, 1.0)
    noise_sd = max(_as_float(row, "noise_sd", 0.20), 0.0)
    abundance_imbalance = max(_as_float(row, "abundance_imbalance", 1.0), 1.0)
    count_dispersion = max(_as_float(row, "count_dispersion", 15.0), 1.0)
    interaction_strength = max(_as_float(row, "interaction_strength", 0.0), 0.0)
    interaction_fraction = np.clip(_as_float(row, "interaction_fraction", 0.0), 0.0, 1.0)
    # Gap #6 real-data confounds. All default to 0.0, so any scenario that does not
    # set them draws no extra RNG and produces byte-identical data to before.
    degradation_3p_bias = max(_as_float(row, "degradation_3p_bias", 0.0), 0.0)
    cell_composition_cv = max(_as_float(row, "cell_composition_cv", 0.0), 0.0)
    batch_effect_sd = max(_as_float(row, "batch_effect_sd", 0.0), 0.0)
    library_depth_cv = max(_as_float(row, "library_depth_cv", 0.0), 0.0)

    n_modules = max(2, min(8, n_genes // 50))
    n_module_genes = max(1, int(round(n_genes * switching_fraction)))
    gene_ids = [f"G{i:05d}" for i in range(n_genes)]
    sample_ids = [f"S{i:04d}" for i in range(n_samples)]

    # Partition module genes into switch-driven vs abundance-driven
    module_mask = np.zeros(n_genes, dtype=bool)
    module_mask[:n_module_genes] = True
    rng.shuffle(module_mask)

    n_abundance_genes = int(round(n_module_genes * abundance_fraction))
    n_switch_genes = n_module_genes - n_abundance_genes

    module_indices = np.where(module_mask)[0]
    rng.shuffle(module_indices)
    abundance_gene_set = set(module_indices[:n_abundance_genes].tolist())

    switching_mask = module_mask.copy()
    abundance_mask = np.zeros(n_genes, dtype=bool)
    for idx in module_indices[:n_abundance_genes]:
        switching_mask[idx] = False
        abundance_mask[idx] = True

    # Dual-signal genes: give a fraction of the switch-driven module genes an
    # abundance signal for the SAME module too (they keep their switch signal).
    if dual_signal_fraction > 0:
        switch_only_idx = module_indices[n_abundance_genes:]
        n_dual = int(round(len(switch_only_idx) * dual_signal_fraction))
        if n_dual > 0:
            dual_idx = rng.choice(switch_only_idx, size=n_dual, replace=False)
            abundance_mask[dual_idx] = True

    module_index = np.full(n_genes, -1, dtype=int)
    all_module_assignments = np.arange(n_module_genes) % n_modules
    rng.shuffle(all_module_assignments)
    for i, gene_idx in enumerate(module_indices):
        module_index[gene_idx] = all_module_assignments[i]

    dx = np.where(np.arange(n_samples) < n_samples / 2, "Control", "SCZD")
    age = np.linspace(25, 85, n_samples) + rng.normal(0, 3, n_samples)
    sex = np.where(np.arange(n_samples) % 2 == 0, "F", "M")
    # sample_table is assembled later, after confound covariates and library size
    # are known (see below).

    module_latent = rng.normal(size=(n_modules, n_samples))
    module_latent += ((age - age.mean()) / age.std())[None, :] * rng.normal(0.15, 0.05, (n_modules, 1))
    module_latent += (dx == "SCZD")[None, :] * rng.normal(0.15, 0.05, (n_modules, 1))

    # PSI signal: only switch genes get latent signal in PSI proportions
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

    # ---- Gap #6: inject real-data confounds into the PSI signal ----
    # Each confound adds structured nuisance variation that is observable via a
    # recorded covariate (RIN, neuron_frac, batch, library_size). isograph_vae_residual
    # regresses these covariates out; isograph_vae does not. That contrast is the
    # WITH/WITHOUT residualization ablation that motivates this gap. Covariates default
    # to constants when their confound is inactive, so non-confound scenarios are
    # unaffected and draw no extra RNG.
    rin = np.full(n_samples, 8.0)
    neuron_frac = np.full(n_samples, 0.5)
    batch_label = np.zeros(n_samples, dtype=int)
    depth_factor: np.ndarray | None = None

    if degradation_3p_bias > 0:
        # Faithful 3' coverage bias (postmortem RNA degradation). Low-RIN samples
        # lose 5' coverage, so each gene's longer/5' isoform loses compositional
        # share. The shift is (a) gene-specific via an all-positive per-gene 3'
        # sensitivity (NOT a single shared linear axis) and (b) one-sided (only
        # degraded, low-RIN samples shift). Because it is not a single RIN axis it
        # is NOT cleanly removed by regressing RIN out post-PC1 -- it genuinely
        # corrupts the switch coordinate, motivating the degradation-aware
        # reliability weighting / multiplex abundance fallback. RIN is recorded so
        # those features (and residualization) can attempt to recover it.
        rin = rng.normal(7.5, 1.2, n_samples)
        deg = -(rin - rin.mean()) / (rin.std() + 1e-8)
        deg_pos = np.clip(deg, 0.0, None)
        sens = np.abs(rng.normal(0.0, 1.0, n_genes))
        signal -= degradation_3p_bias * sens[:, None] * deg_pos[None, :]

    if cell_composition_cv > 0:
        # Bulk brain is a neuron/glia mixture; proportions shift between samples.
        # The centered fraction (sd ~ cell_composition_cv) drives a shared axis.
        neuron_frac = np.clip(rng.normal(0.5, cell_composition_cv, n_samples), 0.05, 0.95)
        comp = neuron_frac - neuron_frac.mean()
        loadings = rng.normal(0.0, 1.0, n_genes)
        signal += 4.0 * loadings[:, None] * comp[None, :]

    if batch_effect_sd > 0:
        # Technical processing-batch structure: a per-gene, per-batch additive
        # offset. Residualizing the batch dummies removes the per-batch means.
        batch_label = rng.integers(0, 4, n_samples)
        batch_offset = rng.normal(0.0, batch_effect_sd, size=(n_genes, 4))
        signal += batch_offset[:, batch_label]

    if library_depth_cv > 0:
        # Per-sample sequencing-depth variation. depth_factor is the log library
        # size; it both biases the signal and sets the count depth below.
        depth_factor = rng.normal(0.0, library_depth_cv, n_samples)
        loadings = rng.normal(0.0, 1.0, n_genes)
        signal += loadings[:, None] * (depth_factor - depth_factor.mean())[None, :]

    baseline_logit = np.log(abundance_imbalance)
    p1 = 1.0 / (1.0 + np.exp(-(signal + baseline_logit)))
    # Abundance-driven genes get PSI ≈ 0.5 (no isoform switching)
    p1[abundance_mask] = 0.5
    p1[~module_mask & ~abundance_mask] = 1.0 / (1.0 + abundance_imbalance)
    p1 = np.clip(p1 + rng.normal(0, noise_sd * 0.15, p1.shape), 1e-4, 1 - 1e-4)

    gene_means = rng.lognormal(mean=np.log(80), sigma=0.8, size=(n_genes, 1))
    if depth_factor is not None:
        sample_depth = np.exp(depth_factor)[None, :]
    else:
        sample_depth = rng.lognormal(mean=0.0, sigma=0.25, size=(1, n_samples))
    library_size = np.log(sample_depth).ravel()
    expected_total = gene_means * sample_depth

    sample_table = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "Dx": dx,
            "Age": age,
            "Sex": sex,
            "RIN": rin,
            "neuron_frac": neuron_frac,
            "batch": batch_label.astype(str),
            "library_size": library_size,
        }
    )

    # Abundance-driven genes: total count modulated by module latent (signal strength ~0.5 SD)
    abundance_multiplier = np.ones((n_genes, n_samples))
    for gene_idx in np.where(abundance_mask)[0]:
        latent = module_latent[module_index[gene_idx]]
        latent_z = (latent - latent.mean()) / (latent.std() + 1e-8)
        abundance_multiplier[gene_idx] = np.exp(0.5 * latent_z)

    gamma_shape = count_dispersion
    gamma_scale = (expected_total * abundance_multiplier) / gamma_shape
    total_rate = rng.gamma(shape=gamma_shape, scale=gamma_scale)
    totals = rng.poisson(np.maximum(total_rate, 1.0)).astype(float)

    n_tx = max(2, n_transcripts_per_gene)
    transcript_counts = np.zeros((n_genes * n_tx, n_samples), dtype=float)
    transcript_rows: list[dict[str, object]] = []
    for gene_idx, gene_id in enumerate(gene_ids):
        # Distribute total counts across n_tx transcripts.
        # For the first two transcripts, use the PSI-driven split; remaining get
        # uniform shares of a small residual so that multi-isoform genes have
        # realistic sparsity without adding independent switching signal.
        total = totals[gene_idx].astype(int)
        tx1 = rng.binomial(total, p1[gene_idx]).astype(float)
        remainder = np.maximum(total - tx1, 0)
        if n_tx == 2:
            tx_counts = [tx1, remainder.astype(float)]
        else:
            # Distribute remainder across the n_tx-1 secondary transcripts via a
            # multinomial draw with equal probabilities.  Each transcript gets an
            # independent vector of per-sample counts that sum to remainder.
            probs = np.full(n_tx - 1, 1.0 / (n_tx - 1))
            secondary = rng.multinomial(remainder.astype(int), probs).T.astype(float)
            # secondary shape: (n_tx-1, n_samples)
            tx_counts = [tx1] + [secondary[i] for i in range(n_tx - 1)]
        for t_idx, tx_vec in enumerate(tx_counts):
            transcript_counts[gene_idx * n_tx + t_idx] = tx_vec
            transcript_rows.append(
                {"transcript_id": f"{gene_id}_T{t_idx + 1}", "gene_id": gene_id, "length": 1000 - t_idx * 50}
            )

    gene_counts = transcript_counts.reshape(n_genes, n_tx, n_samples).sum(axis=1)
    psi = p1
    gene_table = pd.DataFrame({"gene_id": gene_ids})
    transcript_table = pd.DataFrame(transcript_rows)
    psi_table = pd.DataFrame({"psi_uid": [f"PSI_{gene_id}" for gene_id in gene_ids], "gene_id": gene_ids})
    truth_modules = pd.DataFrame(
        {
            "gene_id": [gene_ids[i] for i in module_indices],
            "module_id": [f"M{module_index[i]:03d}" for i in module_indices],
        }
    )
    truth_switch = pd.DataFrame({"gene_id": gene_ids, "has_switch": switching_mask})
    truth_abundance = pd.DataFrame({"gene_id": gene_ids, "has_abundance": abundance_mask})

    extra_feature_specs = [
        build_feature_spec("truth_abundance", "truth_abundance.parquet", truth_abundance),
    ]
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
        ] + extra_feature_specs,
        matrices=[
            build_matrix_spec("gene_counts", "gene_counts.npz", gene_counts),
            build_matrix_spec("transcript_counts", "transcript_counts.npz", transcript_counts),
            build_matrix_spec("psi", "psi.npz", psi),
        ],
        provenance={
            "generator": "isograph_brain_aging_benchmark_v2",
            "scenario": str(row["scenario"]),
            "seed": str(row["seed"]),
            "abundance_fraction": str(_as_float(row, "abundance_fraction", 0.0)),
            "degradation_3p_bias": str(degradation_3p_bias),
            "cell_composition_cv": str(cell_composition_cv),
            "batch_effect_sd": str(batch_effect_sd),
            "library_depth_cv": str(library_depth_cv),
        },
        truth_tables=["truth_modules.parquet", "truth_switch.parquet", "truth_abundance.parquet"],
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
            "truth_abundance": truth_abundance,
        },
        matrices={"gene_counts": gene_counts, "transcript_counts": transcript_counts, "psi": psi},
        truth_tables={
            "truth_modules.parquet": truth_modules,
            "truth_switch.parquet": truth_switch,
            "truth_abundance.parquet": truth_abundance,
        },
    )


def ensure_dataset(row: pd.Series, root: Path) -> Path:
    out = dataset_dir(row, root)
    if (out / "manifest.json").exists():
        return out
    save_dataset_bundle(build_synthetic_bundle(row), out)
    return out
