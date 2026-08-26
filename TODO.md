# Analysis TODO — Cell Genomics

Reconciled against the repository 2026-08-26. Every `[x]` names the committed CLI
and wrapper that closed it, so the claim can be checked rather than trusted.
Status vocabulary: `[x]` done, `[~]` partially done (the gap is stated), `[ ]` open.

Superseded framing note: the old "Additional analysis for Nature Communications"
section is gone — the target is Cell Genomics (retargeted 2026-07-16) and all three
of its items are now covered under Mechanistic validation below.

## Done — do not redo

* [x] **Validate major IsoGraph switches using junction counts or PSI.**
  `validate_switch_splicing.py` + `real_data/brainseq/_h/16.validate_switch_splicing.sh`.
  Caudate marginal OR ~105 vs module OR ~2.3.
* [x] **Adjust human-cohort analyses for estimated cell-type composition.**
  `celltype_composition.py` — joins the committed MuSiC BrainSEQ deconvolution and
  supplies a marker-depletion cut; consumed by `incremental_association`.
* [x] **Reconcile and validate the benchmark recovery-metric definition.**
  `benchmark/partition_metrics.py` + `backfill_metrics.py`, 15 unit tests,
  `benchmark/03_metrics/_m/PARTITION_METRICS.md`. Best-match Jaccard reproduces the
  published `module_recovery_score` to 2.22e-16 over all 13,060 runs, and the
  count-preserving Jaccard null shows WGCNA's negative-control recovery (0.500) is
  fully explained by its module-count structure (z=0.47, perm p=0.746).
* [x] **Direct genetic analysis of module eigengenes.**
  `module_genetic_anchoring.py` + `real_data/brainseq/_h/21.module_anchoring.sh` —
  per-module splicing-specificity contrast against a permutation null of random
  gene sets, closing the "pooling does not show any individual module is anchored"
  objection. Plus S-LDSC + coloc (`ldsc_annot_prep.py`, `coloc_*.py`).
* [x] **Ablations: VAE denoising, estimability weighting, multiplex integration,
  abundance-channel inclusion.** Carried as first-class methods in the synthetic
  grid — `isograph_baseline`, `isograph_spearman_leiden` (Leiden without the VAE),
  `isograph_vae`, `isograph_vae_reliability` (estimability), `isograph_vae_multiplex`,
  `isograph_vae_residual` — plus `project_tiers.py`, which derives the switch_only /
  switch_primary / multiplex edge tiers from ONE fit so the latent is held constant.
  Residualization cost on unconfounded data: `benchmark/03_metrics/_m/RESIDUAL_COST.md`.
* [x] **Sensitivity to counts versus TPM.** N/A — no method ever consumes TPM.
  BrainSEQ uses Salmon counts, GTEx RSEM expected counts, synthetic simulated counts;
  TPM is archived as a data product (`build_parquet.convert_gtex_transcript_tpm`) and
  never fed to a model (`build_bundles.py:230,249`). Manuscript methods corrected.

## Partially done — the gap is the work

* [~] **Confirm top transcript pairs and disease-linked events using long-read brain
  RNA-seq.** Tier-1 done (`longread_switch_confirm.py` +
  `real_data/brainseq/_h/29.longread_switch_confirm.sh`): ONT DLPFC BA9/46, Zenodo
  8180677 Bambu quants, confirms 60.5% of GTEx cortical switch genes.
  **GAP: the 8.4% "switch-like" rate has no matched null, so it is not yet
  interpretable.** Build the null before the number goes in the paper.
* [~] **Frozen-module projection between cohorts.** `replication.py` does bidirectional
  module preservation (best-Jaccard + size-preserving label-permutation null) and
  `scz_age_projection.py` freezes out-of-cohort aging modules and projects them onto
  the independent SCZD cohort. **GAP: no BrainSEQ->GTEx / GTEx->BrainSEQ eigengene
  projection specifically** — preservation-by-overlap is not the same test.
* [~] **Workflow-level comparators.** Done: conventional abundance WGCNA
  (`wgcna_gene`), matched-feature WGCNA (`wgcna_switch_only`, `wgcna_multiplex`,
  `run_matched_wgcna.py`), CLR switch coordinates + Leiden without the VAE
  (`isograph_spearman_leiden`). **GAP: no junction/PSI network comparator, and no
  transcript-level (as opposed to gene-level-on-switch-features) WGCNA.**
* [~] **Select Leiden resolution and edge thresholds using phenotype-blind stability
  criteria.** Resolution 5.0 is canonical and `sweep_leiden.py` +
  `05.sweep_leiden.sh` / `07.sweep_leiden_with_abundance.sh` compute the sweep; the
  giant-module size criterion (>=900 genes) is itself phenotype-blind.
  **GAP: the published justification is phenotype-AWARE** —
  `real_data/gwas/_m/GWAS_RESOLUTION_SUMMARY.md` argues 5.0 by showing it removes a
  GWAS enrichment artifact. A reviewer asked for selection on stability BEFORE any
  trait is consulted. Needed: a phenotype-blind stability-vs-resolution curve with
  the choice made on it, and the GWAS result demoted to post-hoc confirmation.
  Edge thresholds have no documented selection criterion at all.

## Open

* [ ] **Test schizophrenia findings for medication, toxicology, smoking and related
  confounding where available.** Nothing implemented. `qc_covariate_test.py` is
  RNA-quality covariates (RIN / 3' bias / exonic rate), not clinical confounders.
* [ ] **Remaining sensitivity analyses** — no harness exists for any of these:
  * pseudocount choice
  * transcript-expression filtering
  * minor-isoform thresholds
  * transcript number and identifiability
  * quantification pipeline

## Mechanistic validation (was "Additional analysis for Nature Communications")

At least one decisive mechanistic validation was required; the RBP arm has been
carried furthest.

* [x] **Demonstrate direct genetic regulation of a co-switching module.**
  `module_genetic_anchoring.py` (see above).
* [~] **Validate a recurrent predicted RBP regulon using binding or perturbation
  data.** Binding: `rbp_binding.py` (eCLIP, length-normalised density + Haldane OR),
  `rbp_regulon.py` (covariate-adjusted binomial GLM with an estimability gate),
  neuronal-CLIP suite (`neuronal_clip_*.py`, stages 21b-27), `nova2_ctag_clip.py`,
  `nova2_perturbation.py`. Verdict is an **informative sparse null** driven by the
  cross-species mm10/hg38 reciprocal-mapping bottleneck (~3-5% of windows retained).
  Perturbation is designed but not performed:
  `real_data/_m/rbp_target_panel/WETLAB_PERTURBATION_DESIGN.md` with the assayable
  pair list from `rbp_pair_assayability.py` + `real_data/_h/34.rbp_pair_assayability.sh`
  (NONO 11 / ELAVL1 6 / KHDRBS1 8 measurable-and-responsive pairs).
  **GAP: a wet-lab result, which is out of scope for this submission — say so
  explicitly in the limitations rather than leaving it implied.**
  Caveat still to be written into `real_data/_m/rbp/RBP_REGULON_SUMMARY.md`: the
  ENCODE eCLIP panels are HepG2/K562, **not brain**.
* [~] **Validate the shared SNCA alternative-first-exon mechanism across independent
  data types or cohorts.** SNCA carries an LBD splicing colocalization
  (`coloc_isoform_events.py`, `coloc_direction.py`; written up in
  `real_data/_m/GENETIC_ANCHORING_RESULTS.md` and `PER_GENE_DEEP_DIVE_PLAN.md`).
  **GAP: single data type — no orthogonal confirmation of the first-exon switch.**

## Manuscript mechanics (Cell Genomics re-target)

* [ ] Abstract re-lead for a hybrid resource+discovery framing.
* [ ] STAR Methods conversion.
* [ ] Key Resources Table.
* [ ] Fill 5 DOI placeholders.
