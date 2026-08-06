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
    # A1: cis-genetic anchoring of co-switching modules. All default to 0.0/inactive,
    # so any scenario that does not set them draws no extra RNG and produces
    # byte-identical data to before. genetic_effect = per-allele shift of an anchored
    # module's switch latent; genetic_module_fraction = fraction of modules anchored;
    # genetic_maf = minor-allele frequency of the planted variants.
    genetic_effect = max(_as_float(row, "genetic_effect", 0.0), 0.0)
    genetic_module_fraction = np.clip(_as_float(row, "genetic_module_fraction", 0.0), 0.0, 1.0)
    genetic_maf = float(np.clip(_as_float(row, "genetic_maf", 0.25), 0.01, 0.5))

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

    # ---- A1: cis-genetic anchoring of a subset of co-switching modules ----
    # For each anchored module, plant one bi-allelic cis-variant whose per-sample
    # dosage (0/1/2 under Hardy-Weinberg at genetic_maf) shifts that module's SWITCH
    # latent by genetic_effect per allele. The genotype therefore enters only through
    # the module's switch genes' PSI -- the same channel IsoGraph's switch coordinate
    # denoises -- and a method recovers the planted variant->module link only if it
    # first recovers the co-switching module (an eigengene mQTL). We also draw an equal
    # number of UNLINKED "null" variants (added to no module) so a method's spurious
    # anchoring rate is measurable. All RNG is gated on genetic_effect>0, so inactive
    # scenarios draw nothing extra and stay byte-identical.
    genotypes: np.ndarray | None = None
    genetics_meta: list[dict[str, object]] = []
    if genetic_effect > 0 and genetic_module_fraction > 0:
        n_anchor = max(1, int(round(n_modules * genetic_module_fraction)))
        anchored = list(range(n_anchor))  # deterministic first-k modules
        anchor_geno = rng.binomial(2, genetic_maf, size=(n_anchor, n_samples)).astype(float)
        for a, m in enumerate(anchored):
            g = anchor_geno[a]
            gz = (g - g.mean()) / (g.std() + 1e-8)
            module_latent[m] = module_latent[m] + genetic_effect * gz
            genetics_meta.append(
                {"snp_id": f"rsSYN{m:03d}", "module_id": f"M{m:03d}", "anchored": True}
            )
        # Unlinked null variants: same MAF, added to no module (false-positive control).
        null_geno = rng.binomial(2, genetic_maf, size=(n_anchor, n_samples)).astype(float)
        for k in range(n_anchor):
            genetics_meta.append(
                {"snp_id": f"rsNULL{k:03d}", "module_id": "", "anchored": False}
            )
        genotypes = np.vstack([anchor_geno, null_geno])

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
    median_tin: np.ndarray | None = None

    if degradation_3p_bias > 0:
        # Faithful 3' coverage bias (postmortem RNA degradation). A latent per-sample
        # degradation severity drives 5' coverage loss, so each gene's longer/5'
        # isoform loses compositional share. The shift is (a) gene-specific via an
        # all-positive per-gene 3' sensitivity and (b) one-sided (only degraded
        # samples shift). It is NOT a single recorded axis, so regressing RIN out
        # post-PC1 cannot cleanly remove it -- it genuinely corrupts the switch
        # coordinate, motivating degradation-aware reliability / abundance fallback.
        deg_latent = rng.normal(0.0, 1.0, n_samples)
        deg_pos = np.clip(deg_latent - deg_latent.mean(), 0.0, None)
        sens = np.abs(rng.normal(0.0, 1.0, n_genes))
        signal -= degradation_3p_bias * sens[:, None] * deg_pos[None, :]
        # Two OBSERVED covariates of the same latent severity, differing in fidelity:
        #   RIN        - blunt bench RNA-integrity score; weak, noisy anti-correlation
        #                with the true coverage artifact (the usual postmortem case).
        #   median_tin - RSeQC median Transcript Integrity Number, a direct gene-body
        #                3' coverage metric; a SHARPER observation of the same artifact.
        # median_tin lets switch-reliability localize degradation better than RIN while
        # staying a single sample-level vector (approach #3). Both are drawn from the
        # main rng so a fixed seed reproduces them exactly.
        rin = np.clip(7.5 - 1.2 * deg_latent + rng.normal(0.0, 1.1, n_samples), 1.0, 10.0)
        median_tin = 75.0 - 12.0 * deg_pos + rng.normal(0.0, 2.0, n_samples)

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
    # median_tin is recorded only under active 3' degradation, so non-degradation
    # scenarios keep an unchanged sample_table schema (byte-identical parquet).
    if median_tin is not None:
        sample_table["median_tin"] = median_tin

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

    # Transcript-level switch-event truth (for interpretation accuracy benchmarking).
    # By construction transcript T1 carries the PSI signal p1 (driven by the gene's
    # module_latent); the remaining transcripts split a residual. So T1 is the
    # ground-truth switch-driver isoform, and the switch amplitude is the spread of
    # p1 across samples. Computed purely from already-drawn quantities (no RNG draws),
    # so the count/PSI matrices and every existing dataset stay byte-identical; this
    # table is only DISCRIMINATIVE when n_transcripts_per_gene >= 3 (at n_tx=2 the two
    # isoforms are mirror images). switch_transcript_id matches the transcript_table
    # ids ("{gene_id}_T1").
    switch_gene_idx = np.where(switching_mask)[0]
    truth_switch_event = pd.DataFrame(
        {
            "gene_id": [gene_ids[i] for i in switch_gene_idx],
            "module_id": [f"M{module_index[i]:03d}" for i in switch_gene_idx],
            "switch_transcript_id": [f"{gene_ids[i]}_T1" for i in switch_gene_idx],
            "n_transcripts": [n_tx] * len(switch_gene_idx),
            "true_delta_psi": [float(p1[i].max() - p1[i].min()) for i in switch_gene_idx],
            "true_psi_std": [float(p1[i].std()) for i in switch_gene_idx],
        }
    )

    extra_feature_specs = [
        build_feature_spec("truth_abundance", "truth_abundance.parquet", truth_abundance),
        build_feature_spec("truth_switch_event", "truth_switch_event.parquet", truth_switch_event),
    ]
    matrix_specs = [
        build_matrix_spec("gene_counts", "gene_counts.npz", gene_counts),
        build_matrix_spec("transcript_counts", "transcript_counts.npz", transcript_counts),
        build_matrix_spec("psi", "psi.npz", psi),
    ]
    provenance = {
        "generator": "isograph_brain_aging_benchmark_v2",
        "scenario": str(row["scenario"]),
        "seed": str(row["seed"]),
        "abundance_fraction": str(_as_float(row, "abundance_fraction", 0.0)),
        "degradation_3p_bias": str(degradation_3p_bias),
        "cell_composition_cv": str(cell_composition_cv),
        "batch_effect_sd": str(batch_effect_sd),
        "library_depth_cv": str(library_depth_cv),
    }
    truth_table_names = [
        "truth_modules.parquet", "truth_switch.parquet", "truth_abundance.parquet",
        "truth_switch_event.parquet",
    ]
    feature_tables = {
        "gene": gene_table,
        "transcript": transcript_table,
        "psi": psi_table,
        "truth_module": truth_modules,
        "truth_switch": truth_switch,
        "truth_abundance": truth_abundance,
        "truth_switch_event": truth_switch_event,
    }
    matrices = {"gene_counts": gene_counts, "transcript_counts": transcript_counts, "psi": psi}
    truth_tables = {
        "truth_modules.parquet": truth_modules,
        "truth_switch.parquet": truth_switch,
        "truth_abundance.parquet": truth_abundance,
        "truth_switch_event.parquet": truth_switch_event,
    }

    # A1: genetics artifacts are attached only when active, so all pre-existing
    # scenarios keep an unchanged manifest/file set (byte-identical datasets). Row i
    # of the genotypes matrix corresponds to row i of truth_genetics (anchored SNPs
    # first, then the unlinked null SNPs).
    if genotypes is not None:
        truth_genetics = pd.DataFrame(
            [
                {
                    **meta,
                    "maf": genetic_maf,
                    "per_allele_effect": genetic_effect,
                    "n_module_genes": (
                        int((module_index == int(str(meta["module_id"])[1:])).sum())
                        if meta["anchored"] else 0
                    ),
                    "n_module_switch_genes": (
                        int(((module_index == int(str(meta["module_id"])[1:])) & switching_mask).sum())
                        if meta["anchored"] else 0
                    ),
                }
                for meta in genetics_meta
            ]
        )
        # truth_genetics is persisted purely via the truth_tables mechanism
        # (save_dataset_bundle writes bundle.truth_tables by filename; load reads them
        # from manifest.truth_tables). It needs no FeatureTableSpec -- whose `kind` is a
        # closed Literal in the isograph package -- so A1 stays decoupled from the
        # package schema. The genotypes matrix uses a free-form MatrixSpec assay name.
        matrix_specs.append(build_matrix_spec("genotypes", "genotypes.npz", genotypes))
        provenance.update(
            {
                "genetic_effect": str(genetic_effect),
                "genetic_module_fraction": str(genetic_module_fraction),
                "genetic_maf": str(genetic_maf),
            }
        )
        truth_table_names.append("truth_genetics.parquet")
        matrices["genotypes"] = genotypes
        truth_tables["truth_genetics.parquet"] = truth_genetics

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
        matrices=matrix_specs,
        provenance=provenance,
        truth_tables=truth_table_names,
    )
    return DatasetBundle(
        manifest=manifest,
        sample_table=sample_table,
        feature_tables=feature_tables,
        matrices=matrices,
        truth_tables=truth_tables,
    )


def ensure_dataset(row: pd.Series, root: Path) -> Path:
    out = dataset_dir(row, root)
    if (out / "manifest.json").exists():
        return out
    save_dataset_bundle(build_synthetic_bundle(row), out)
    return out
