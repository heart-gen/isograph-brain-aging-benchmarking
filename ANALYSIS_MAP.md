# Analysis map

One row per analysis: what it asks, the CLI that implements it, the wrapper that runs
it, where its outputs land, and which display item it backs. Stage READMEs carry the
interpretation and the caveats; this file is the index.

Every analysis is a committed, parametrized CLI under `isograph_benchmark/` with a
deterministic seed and outputs written to disk — no ad-hoc inline scripts for anything
that reaches the paper.

## 01 — Synthetic benchmark

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| Synthetic grid | `benchmark/run_synthetic.py`, `run_one.py` | `01_synthetic_benchmark/01_synthetic/_h/` | `01_synthetic/_m/synthetic_results.parquet` | Fig 1 |
| Module interpretation accuracy | `benchmark/interpret_modules.py` | `02_interpret/_h/` | `02_interpret/_m/` | S11 |
| Metrics + pairwise tests | `stats/summarize.py`, `sweep_breakdown.py` | `03_metrics/_h/` | `03_metrics/_m/` | Fig 1, S1–S10, `tableS_benchmark_summary` |
| Confound ablation | `benchmark/residual_cost.py` | `01_synthetic/_h/run_residual_cost.sh` | `01_synthetic/_m/datasets_residual_cost/` | S8, S9 |
| **GPU VAE reproducibility probe** | `benchmark/gpu_reproducibility_probe.py` (`select`, `compare`) | `01_synthetic/_h/run_gpu_repro_probe.sh` | `03_metrics/_m/gpu_repro_probe/`; isolated re-runs in `01_synthetic/_o/gpu_repro_probe/` | *(GPU arm deterministic; stored rows not reproducible from current code — archive is the record; telemetry now logs git commits)* |
| **Method run counts** | `stats/summarize.py` | `03_metrics/_h/step_2_summarize.sh` | `03_metrics/_m/tableS_method_run_counts.csv` | Supplement (ablation-arm power) |

## 02 — Module discovery

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| IsoGraph fits | `real_data/run_models.py` | `02_module_discovery/_h/01–03` | `<cohort>/<region>/_m/isograph_vae/` | all *(switching transcript filter since 2026-09-14; the with-abundance arm `_h/04`/`_h/06` is retired — outputs only at tag `legacy_expression_filter`)* |
| Leiden resolution sweep | `real_data/sweep_leiden.py` | `_h/05–06` (BrainSEQ), `_h/15` (GTEx) | `<store>/isograph_vae/leiden_sweep_results.parquet` | S-real-2 *(disclosed sensitivity; 5.0 stays canonical on the ≥ 900-gene criterion, PI decision 2026-09-12)* |
| Classical WGCNA baseline | `_h/07–09.wgcna_gene_*.R` | `_h/07–09` | `<store>/wgcna_gene/` | S-real-1 |
| Matched-feature WGCNA baselines | `real_data/run_matched_wgcna.py` | `_h/10–11` (aging), `_h/16` (`caudate_sczd`, added 2026-09-12) | `<store>/wgcna_{switch_only,multiplex}/` | Fig 3 (internal control); S-real-1 / Table S1 |
| Refit / reprojection QC | `real_data/qc_covariate_test.py`, `tier_checks.py` | `_h/12–14` | — | — *(**retired 2026-09-14**: model-design QC, cited nowhere and never full-coverage; not re-run on the switching filter, outputs only at tag `legacy_expression_filter`)* |
| **Module sizes, all methods** | `real_data/module_sizes.py` | *(login; after the fits)* | `_m/module_sizes/{module_sizes.tsv,module_size_summary.tsv}` | *(giant modules occur in both IsoGraph GTEx and WGCNA; reported and compared, not engineered away — PI 2026-09-14)* |
| **Partition provenance guard** | `real_data/partition_provenance.py` (library: `partition_fingerprint`, `check_partition`, `load_enrichment`) | *(called by `module_enrichment` and ten consumers)* | sha256 fingerprint in each `<store>/module_enrichment/*_modules.parquet` file metadata | *(gate: a module-id join is only valid against the fit that wrote it)* |

## 03 — Module trust

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| Split-half stability | `real_data/stability.py` | `03_module_trust/_h/01–03` | `_m/stability/{partitions/,stability_summary.parquet}` | Fig 2 |
| Per-module trust funnel | `real_data/module_trust.py` | `_h/04, 09` | `_m/stability/{module_trust,modules_meta}/` | Fig 2, S7 |
| **Within-cohort Q2 drivers / Q3 sign** | `real_data/module_trust.py within` | `_h/17` (after `_h/04`) | `_m/stability/module_trust/within_cohort__<cohort>__<region>__<method>.parquet` | Fig 2B |
| **Age-model curvature test** | `real_data/age_model_curvature.py` | `_h/14` | `_m/stability/module_trust/age_model_curvature__isograph.{parquet,json}` | *(licenses the linear age model behind 23/130)* |
| **Cross-cohort eigengene projection** | `real_data/eigengene_projection.py` (`run`, `aggregate`) | `_h/18` (array: {isograph, wgcna} × 3 pairs; `PROJECTION_AGGREGATE=1`) | `_m/stability/eigengene_projection/{eigengene_projection__<pair>__<method>__<direction>.parquet,switch_axis_alignment__<pair>.parquet,EIGENGENE_PROJECTION.md}` | *(frozen-weight preservation vs a type-matched null; module-specific age z against same-weight random projections — raw projected age r is cohort-confounded)* |
| **Pooled cross-cohort replication** | `real_data/module_trust.py replication-pooled` | `_h/15` | `_m/stability/module_trust/module_aging_replication_pooled__<method>{.parquet,__stats.json}` | *(null for both methods; reported, not a figure)* |
| **Phenotype-blind resolution sweep** | `real_data/stability.py sweep` | `_h/16` (after `_h/01` with `STABILITY_SAVE_EDGES=1`), then `_h/03` | `_m/stability/stability_summary.parquet` (rows with method isograph_resXpY) | *(disclosed sensitivity; cannot select a resolution)* |
| LR / software robustness | `real_data/stability.py` | `_h/05–07` | — | — *(**retired 2026-09-14**: settled the single-LR promotion; not re-run on the switching filter, outputs only at tag `legacy_expression_filter`)* |
| Giant-cap ablation | `real_data/stability.py` | `_h/08` | `_m/stability_gcap_ab/` | — *(**retired 2026-09-12**: a resolution-2.0 A/B of a cap never promoted; superseded by resolution 5.0, not cited)* |
| Cross-cohort replication | `real_data/replication.py`, `replication_go.py` | `_h/10–11` | `_m/replication/` | Fig 2 |
| Replication null + function | `real_data/replication_permutation.py`, `replication_functional.py` | `_h/12–13` | `_m/replication/FUNCTIONAL_PRESERVATION*` | Fig 2 |

## 04 — Module characterization

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| Module interpretation | `real_data/interpret_modules.py` | `04_module_characterization/_h/01–02` | `<store>/isograph_vae/module_interpret/` | S-real-5 |
| GO:BP enrichment | `real_data/module_enrichment.py`, `go_enrichment.py` | `_h/03–04` | `<store>/module_enrichment/` | S-real-1 |
| GO-invisible gate | `real_data/go_invisible_gate.py` | `_h/05` | `<caudate_sczd store>/go_invisible_gate.parquet` | S-real-3, S6 |
| Incremental association | `real_data/incremental_association.py` | `_h/06–07` | `<store>/isograph_vae/incremental_association/` | S-real-4 |
| Abundance/switch separation | `real_data/abundance_structure_separation.py` | `_h/08` | `_m/incremental_effect_sizes.parquet` | S-real-4 |
| Composition-unique genes | `real_data/characterize_composition_unique.py` | `_h/09` | `_m/composition_adjustment.parquet` | S-real-4 |
| Cell-type composition | `real_data/celltype_composition.py` + MuSiC R | `_h/10–11` | `_m/composition_adjustment.parquet`, `02_module_discovery/gtex/_m/composition/` | **Fig 5** `figCompositionRobustness`, Table S13 |
| Tier projection | `real_data/project_tiers.py` | `_h/12–13` | — | — *(**retired 2026-09-14**: channel-tier ablation, cited nowhere, 6/17 coverage; not re-run on the switching filter, outputs only at tag `legacy_expression_filter`)* |
| Three-baseline comparison | `real_data/baseline_comparison.py` | `_h/14` | `_m/baseline_comparison/` | S-real-1, S1/S2 |

## 05 — Genetic anchoring

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| xQTL anchoring (17 analyses) | `real_data/qtl_anchoring.py` | `05_genetic_anchoring/_h/01–02` | `<store>/qtl_anchoring.parquet` | Fig 3 |
| xQTL anchoring sensitivity (17 x 3 arms) | `real_data/qtl_anchoring.py --outcome/--covariate-set` | `05_genetic_anchoring/_h/17` | `<store>/qtl_anchoring_{constraint,continuous,dose}.parquet` | Table S19 |
| Contrast meta-analysis | `real_data/qtl_anchoring_meta.py` | (login) | `_m/qtl_anchoring_meta/` | **Fig 3, Table 1**, S3–S5 |
| Module-level coloc convergence (5 traits) | `real_data/module_coloc_convergence.py` | `05_genetic_anchoring/_h/18` | `_m/module_coloc_convergence/` | S-real-10, Tables S20a/b |
| Sensitivity meta-analysis | `real_data/qtl_anchoring_meta.py --outcome/--covariate-set` | (login) | `_m/qtl_anchoring_meta/sensitivity/<arm>/` | Table S19 |
| sQTL direction concordance | `real_data/sqtl_concordance.py`, `sqtl_concordance_meta.py` | `_h/03` | `_m/sqtl_concordance_meta/` | — |
| Module genetic anchoring | `real_data/module_genetic_anchoring.py` | `_h/04` | `_m/module_genetic_anchoring_meta/` | Table S16 only — **no figure by design** (pooled-gene, not per-module) |
| MAGMA module-GWAS | `gwas/prepare_magma_inputs.py` | `_h/05–07` | `_m/gwas/magma_results_combined.parquet` | S-real-2 |
| Colocalization (SuSiE + eCAVIAR) | `real_data/coloc_{prep,summary,isoform_events,direction,meta}.py` | `_h/08–11` | `_m/coloc/` | Fig 4, S8/S9 |
| **Signal-level coloc (`coloc.susie`)** | `real_data/coloc_signal_susie.py` + `_h/22.coloc_gwas_susie.R`, `_h/23.coloc_signal_susie.R` | `_h/22–23` | `_m/coloc_signal_susie/` | *(estimator hierarchy: susie > abf > CLPP)* |
| **Locus event audit** | `real_data/locus_event_audit.py` (`configs/known_splice_events.yaml`) | `_h/24` | `_m/locus_event_audit/<nominations>/` | *(tiering of UNC13A / PICALM / SNCA)* |
| **BrainSEQ switch-QTL (S_g) vs abundance-QTL (A_g)** | `real_data/brainseq_switch_qtl.py` | `_h/25` (GPU) | `_m/brainseq_switch_qtl/{all_samples,ea_only}/{caudate,dlpfc,hippocampus}/qtl/cis_qtl_switch.parquet`, `_m/brainseq_switch_qtl/{all_samples,ea_only}/modality_contrast.parquet` | *(same-tissue genetic anchoring, not replication)* |
| **BrainSEQ switch vs abundance effect sizes, selection-symmetric** | `real_data/brainseq_switch_qtl.py --stage effect_size` | `_h/31` (array over arms) | `_m/brainseq_switch_qtl/{all_samples,ea_only}/{effect_size_symmetric.parquet,lead_variant_frequency.parquet,EFFECT_SIZE_SYMMETRIC.md}` | *(every gene tested on both axes; own-lead, cross-lead and common-variant contrasts; |slope| and conservative |z| both reported, with lead-variant MAF)* |
| **BrainSEQ QTL checks: S_g reproduction / sign pin, A_g positive control** | `real_data/brainseq_qtl_checks.py` | `_h/28` | `_m/brainseq_switch_qtl/{all_samples,ea_only}/checks` | *(gate before any BrainSEQ QTL result is read)* |
| **Locus LD robustness audit** | `_h/26.locus_ld_robustness.R` | `_h/26` | `_m/locus_ld_robustness/<analysis>__<LOCUS_ID>/` | *(required before a recovered or headline locus carries a claim)* |
| **SMR + HEIDI** | `real_data/smr_heidi.py` | `_h/27` compute (`SMR_QTL_SOURCE=brainseq,SMR_ARM=ea_only` for BrainSEQ; sensitivity arms `SMR_PEQTL=1e-6` and `SMR_MULTI=1`, both with `SMR_SKIP_BESD=1`), `_h/30` meta (`--peqtl-smr 1e-6` / `--smr-multi`, pass-through) | `_m/smr_heidi/{gtex,brainseq/ea_only}`, each with `_m/smr_heidi/{gtex,brainseq/ea_only}/sensitivity/{peqtl_smr_1e-06,smr_multi}` | *(orthogonal corroboration beneath coloc; BrainSEQ gated on the `_h/28` checks and read against BrainSEQ's own coloc; two separately-corrected testing families, `smr_status` separates untested from tested-null, instrument F reported, attrition trace shows GWAS coverage — not allele QC — is the dominant filter)* |
| **BrainSEQ signal-level coloc (S_g / A_g, EA-only)** | `real_data/coloc_brainseq.py` + `_h/29.coloc_brainseq_susie.R` | `_h/29` | `_m/coloc_brainseq/ea_only/` | *(same-tissue genetic anchoring; in-sample QTL LD; gated on the `_h/28` checks)* |
| Per-gene sQTL-vs-eQTL coloc contrast | `real_data/coloc_modality_contrast.py` + `_h/20.coloc_modality_abf.R` | `_h/19–21` | `_m/coloc_modality_contrast/` | *(run 2026-09-01, 4 arms — **null**; reported as a negative control, no figure)* |
| S-LDSC partitioned heritability | `real_data/ldsc_annot_prep.py`, `ldsc_summary.py` | `_h/12–14` | `_m/ldsc/ldsc_partitioned.parquet` | Fig 4B |
| Per-gene deep dive | `real_data/gene_deep_dive.py` | `_h/15` | `_m/deep_dive/` | **Table 2**, S8–S12 |
| SCZ age projection | `real_data/scz_age_projection.py` | `_h/16` | `_m/scz_age_projection/` | **Fig 4E** (folded from `figSczConvergence`), Table S18 |

## 06 — Switch mechanism

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| Short-read junction confirmation (SNCA/CTSH) | `real_data/junction_coloc_confirm.py` | `06_switch_mechanism/_h/12` | `_m/junction_coloc_confirm/` | Fig 4A (confirms SNCA; CTSH withheld) |
| Switch consequence + meta | `real_data/switch_consequence[_meta].py` | `06_switch_mechanism/_h/01–02` | `_m/switch_consequence_meta.parquet` | S-real-5 |
| PSI / junction validation | `real_data/validate_switch_splicing.py` | `_h/03–04` | `_m/switch_validation/` | — |
| **phASER allelic feasibility screen** | `real_data/ase_switch_direction.py` | *(login)* | `_m/ase_switch_direction/<region>/` | *(gate: is an allele-specific switch test supported at all?)* |
| Orthogonal confirmation | `real_data/switch_orthogonal_confirm.py` | `_h/05` | `_m/switch_orthogonal_confirm/` | **S-real-8** `figOrthogonalConfirm`, Table S14 |
| **Orthogonal confirmation, signal-level coloc** | `real_data/coloc_isoform_events.py --layer signal` → `switch_orthogonal_confirm.py --events signal` | `_h/05 --events signal` | `05_genetic_anchoring/_m/coloc_signal_susie/all_introns/coloc_isoform_events.parquet`, `_m/switch_orthogonal_confirm/signal_coloc/` | *(Analysis 6: long-read check over the coloc.susie nominations; the CLPP arm above is unchanged)* |
| ISA / satuRn concordance | `real_data/isa_concordance.py` | `_h/06` | `_m/isa_concordance/` | **S-real-9** `figIsaConcordance`, Table S15 |
| Long-read confirmation | `real_data/longread_switch_confirm.py` | `_h/07` | `_m/longread_switch_confirm/` | — |
| Clinical consequence | `real_data/clinical_consequence[_meta].py` | `_h/08–09` | `_m/clinical_consequence_meta.parquet` | S-real-7 |
| SCZ confound sensitivity | `real_data/scz_confound_sensitivity.py` | `_h/10` | `_m/scz_confound_sensitivity/` | — |
| Feature sensitivity (5 axes, partition fixed) | `real_data/switch_feature_sensitivity.py` | `_h/11` (one region; `--aggregate` for the cross-region report + quantification axis), `_h/13` (13 GTEx regions) | `_m/switch_feature_sensitivity/{SWITCH_FEATURE_SENSITIVITY.md,sensitivity_summary_all.parquet,quantification_concordance.parquet}`; per region under `_m/switch_feature_sensitivity/brainseq/caudate` and `_m/switch_feature_sensitivity/gtex/*` | — |
| **Feature sensitivity, full refit per setting** | `real_data/switch_feature_refit.py` | `_h/14` (array, one fit per setting; `REFIT_AGGREGATE=1` to combine) | `_m/switch_feature_sensitivity/refit/brainseq/caudate/{refit_summary.parquet,REFIT_SENSITIVITY.md}` (one subdirectory per setting) | *(published-setting refit is the noise floor)* |
| **PSI / junction validation — interpretation** | *(written from `summary_*.json`)* | — | `_m/switch_validation/SWITCH_VALIDATION_SUMMARY.md`, `_m/longread_switch_confirm/LONGREAD_SWITCH_CONFIRM.md` | — |

## 07 — RBP regulation

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| Motif families | `real_data/rbp_motif_families.py` | `07_rbp_regulation/_h/01` | `_m/rbp/rbp_motif_families.parquet` | S-real-6 |
| Regulons (mature + intronic) | `real_data/rbp_scan[_intronic].py`, `rbp_regulon.py` | `_h/02–03` | `_m/rbp/rbp_regulon*.parquet` | S-real-6, S10 |
| CLIP binding evidence | `real_data/rbp_binding.py` | `_h/04–05` | `_m/rbp/rbp_binding*.parquet` | S-real-6 panel B, Table S17 |
| Target panel / assayability | `real_data/rbp_pair_assayability.py`, `rbp_target_panel.py` | `_h/06–07` | `_m/rbp_target_panel/` | — |
| Neuronal CLIP | `real_data/neuronal_clip_*.py` | `_h/08–11` | `_m/neuronal_clip/` | — |
| NOVA family / NOVA2 | `real_data/nova_family_renomination.py`, `nova2_*.py` | `_h/12–14` | `_m/neuronal_clip/nova*/` | — |

`<store>` = `02_module_discovery/<cohort>/<region>/_m/`, the shared cohort × region
artifact store.
