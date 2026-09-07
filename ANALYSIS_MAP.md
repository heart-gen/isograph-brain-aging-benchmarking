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

## 02 — Module discovery

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| IsoGraph fits | `real_data/run_models.py` | `02_module_discovery/_h/01–04` | `<cohort>/<region>/_m/isograph_vae[_with_abundance]/` | all |
| Leiden resolution sweep | `real_data/sweep_leiden.py` | `_h/05–06` | `<store>/isograph_vae/leiden_sweep_results.parquet` (BrainSEQ only) | S-real-2 |
| Classical WGCNA baseline | `_h/07–09.wgcna_gene_*.R` | `_h/07–09` | `<store>/wgcna_gene/` | S-real-1 |
| Matched-feature WGCNA baselines | `real_data/run_matched_wgcna.py` | `_h/10–11` | `<store>/wgcna_{switch_only,multiplex}/` | Fig 3 (internal control) |
| Refit / reprojection QC | `real_data/qc_covariate_test.py`, `tier_checks.py` | `_h/12–14` | `<store>/{qc_covariate_test,tier_checks}/` | — |

## 03 — Module trust

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| Split-half stability | `real_data/stability.py` | `03_module_trust/_h/01–03` | `_m/stability/{partitions/,stability_summary.parquet}` | Fig 2 |
| Per-module trust funnel | `real_data/module_trust.py` | `_h/04, 09` | `_m/stability/{module_trust,modules_meta}/` | Fig 2, S7 |
| LR / software robustness | `real_data/stability.py` | `_h/05–07` | `_m/stability/lr_validation/` | — |
| Giant-cap ablation | `real_data/stability.py` | `_h/08` | `_m/stability_gcap_ab/` | — |
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
| Tier projection | `real_data/project_tiers.py` | `_h/12–13` | `<store>/tier_checks/` | — |
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
| Orthogonal confirmation | `real_data/switch_orthogonal_confirm.py` | `_h/05` | `_m/switch_orthogonal_confirm/` | **S-real-8** `figOrthogonalConfirm`, Table S14 |
| ISA / satuRn concordance | `real_data/isa_concordance.py` | `_h/06` | `_m/isa_concordance/` | **S-real-9** `figIsaConcordance`, Table S15 |
| Long-read confirmation | `real_data/longread_switch_confirm.py` | `_h/07` | `_m/longread_switch_confirm/` | — |
| Clinical consequence | `real_data/clinical_consequence[_meta].py` | `_h/08–09` | `_m/clinical_consequence_meta.parquet` | S-real-7 |
| SCZ confound sensitivity | `real_data/scz_confound_sensitivity.py` | `_h/10` | `_m/scz_confound_sensitivity/` | — |
| Feature sensitivity (5 axes) | `real_data/switch_feature_sensitivity.py` | `_h/11` | `_m/switch_feature_sensitivity/` | — |

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
