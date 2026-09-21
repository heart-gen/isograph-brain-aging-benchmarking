# Analysis map

One row per analysis: what it asks, the CLI that implements it, the wrapper that runs
it, where its outputs land, and which display item it backs. Stage READMEs carry the
interpretation and the caveats; this file is the index.

Every analysis is a committed, parametrized CLI under `isograph_benchmark/` with a
deterministic seed and outputs written to disk — no ad-hoc inline scripts for anything
that reaches the paper.

Wrapper names are `<tier><letter>.<name>.sh` under each stage's `_h/`: the tier is the run
order within the stage, and steps sharing a tier run in parallel. `<stage>/_h/run_stage.sh`
submits a stage in that order and `run_pipeline.sh` chains the stages (README "Running").

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
| IsoGraph fits | `real_data/run_models.py` | `02_module_discovery/_h/01a–01c` (same wrappers with `--leiden-resolution 2.0` for `isograph_vae_res2`) | `<cohort>/<region>/_m/isograph_vae/` | all *(switching transcript filter since 2026-09-14; the with-abundance arm is retired to `_h/retired/` — outputs only at tag `legacy_expression_filter`)* |
| Leiden resolution sweep | `real_data/sweep_leiden.py` | `_h/02a` (BrainSEQ), `_h/02b` (GTEx) | `<store>/isograph_vae/leiden_sweep_results.parquet` | S-real-2 *(disclosed sensitivity; the sweep cannot select a resolution — production is **2.0** since the PI decision of 2026-09-16, on the ≥ 900-gene criterion plus coverage, and 5.0 is the `isograph_vae_res5` arm)* |
| Classical WGCNA baseline | `_h/01d–01f.wgcna_gene_*.R` | `_h/01d–01f` | `<store>/wgcna_gene/` | S-real-1 |
| Matched-feature WGCNA baselines | `real_data/run_matched_wgcna.py` | `_h/01g` (BrainSEQ aging), `_h/01h` (`caudate_sczd`, added 2026-09-12), `_h/01i` (GTEx) | `<store>/wgcna_{switch_only,multiplex}/` | Fig 3 (internal control); S-real-1 / Table S1 |
| Refit / reprojection QC | `real_data/qc_covariate_test.py`, `tier_checks.py` | `_h/retired/{refit_qc_brainseq,refit_qc_gtex,reproject_qc_gtex}.sh` | — | — *(**retired 2026-09-14**: model-design QC, cited nowhere and never full-coverage; outputs only at tag `legacy_expression_filter`)* |
| **Module sizes, all methods** | `real_data/module_sizes.py` | `_h/02c` (after every fit) | `_m/module_sizes/{module_sizes.tsv,module_size_summary.tsv}` | *(giant modules occur in both IsoGraph GTEx and WGCNA; reported and compared, not engineered away — PI 2026-09-14)* |
| **Partition provenance guard** | `real_data/partition_provenance.py` (library: `partition_fingerprint`, `check_partition`, `load_enrichment`) | *(called by `module_enrichment` and ten consumers)* | sha256 fingerprint in each `<store>/module_enrichment/*_modules.parquet` file metadata | *(gate: a module-id join is only valid against the fit that wrote it)* |

## 03 — Module characterization

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| Module interpretation | `real_data/interpret_modules.py` | `03_module_characterization/_h/01a–01b` | `<store>/isograph_vae/module_interpret/` | S-real-5 |
| GO:BP enrichment | `real_data/module_enrichment.py`, `go_enrichment.py` | `_h/01c–01d` | `<store>/module_enrichment/` | S-real-1 |
| Incremental association | `real_data/incremental_association.py` | `_h/01e–01f` | `<store>/isograph_vae/incremental_association/` | S-real-4 |
| Cell-type composition | `real_data/celltype_composition.py` + MuSiC R (`gtex_music_deconv.R`) | `_h/01g–01h`, then `_h/02c` (`meta` rollup) | `_m/composition_adjustment.parquet`, `02_module_discovery/gtex/_m/composition/` | **Fig 5** `figCompositionRobustness`, Table S13 |
| GO-invisible gate | `real_data/go_invisible_gate.py` | `_h/02a` | `<caudate_sczd store>/go_invisible_gate.parquet` | S-real-3, S6 |
| Composition-unique genes | `real_data/characterize_composition_unique.py` | `_h/02b` | `_m/composition_adjustment.parquet` | S-real-4 |
| Abundance/switch separation | `real_data/abundance_structure_separation.py` | `_h/03a` | `_m/incremental_effect_sizes.parquet` | S-real-4 |
| Tier projection | `real_data/project_tiers.py` | `_h/retired/{project_tiers_pilot,tiers_fanout}.sh` | — | — *(**retired 2026-09-14**: channel-tier ablation, cited nowhere, 6/17 coverage; outputs only at tag `legacy_expression_filter`)* |

## 04 — Module trust

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| Split-half stability | `real_data/stability.py` | `04_module_trust/_h/01a–01b`, then `_h/03a` (aggregate) | `_m/stability/{partitions/,stability_summary.parquet}` | Fig 2 |
| Per-module trust funnel | `real_data/module_trust.py` | `_h/02b` (meta), `_h/02c` (gate), `_h/03c` (replication) | `_m/stability/{module_trust,modules_meta}/` | Fig 2, S7 |
| **Q1 trust gate** | `real_data/module_trust.py stability` | `_h/02c` (array: 6 regions × {isograph, wgcna}; before `_h/03c–03g`) | `_m/stability/module_trust/module_stability__<cohort>__<region>__<method>.parquet` | Fig 2 |
| **Q4 complementarity** | `real_data/module_trust.py complementarity` | `_h/03g` (array: 6 regions × {isograph, wgcna}; after `_h/02c` and stage 03 `_h/01a–01b`, `_h/02b`) | `_m/stability/module_trust/module_complementarity__<cohort>__<region>__<method>.parquet` | S7 |
| **Linear-vs-spline Q3 decomposition** | `real_data/module_trust.py replication-model-contrast` | `_h/04a` (after `_h/03c` linear and `--model spline`) | `_m/stability/module_trust/replication_model_contrast__<method>.parquet` | Fig 2C |
| **Within-cohort Q2 drivers / Q3 sign** | `real_data/module_trust.py within` | `_h/03b` (after `_h/02b`) | `_m/stability/module_trust/within_cohort__<cohort>__<region>__<method>.parquet` | Fig 2B |
| **Age-model curvature test** | `real_data/age_model_curvature.py` | `_h/01d` | `_m/stability/module_trust/age_model_curvature__isograph.{parquet,json}` | *(licenses the linear age model behind 23/130)* |
| **Cross-cohort eigengene projection** | `real_data/eigengene_projection.py` (`run`, `aggregate`) | `_h/03f` (array: {isograph, wgcna} × 3 pairs), then `_h/04c` (aggregate) | `_m/stability/eigengene_projection/{eigengene_projection__<pair>__<method>__<direction>.parquet,switch_axis_alignment__<pair>.parquet,EIGENGENE_PROJECTION.md}` | *(frozen-weight preservation vs a type-matched null; module-specific age z against same-weight random projections — raw projected age r is cohort-confounded)* |
| **Module context and held-out DTU evidence (PI 12a)** | `real_data/module_dtu_added_value.py` (`saturn`, `analyze`, `summarize`) | `_h/02e` (array: 6 regions × 5 seeds, satuRn on each split half), `_h/03i` (array: 6 regions), then `_h/04d` (giant-module supplement + summary) | `_m/dtu_added_value/{<cohort>__<region>__{replicates,null}.parquet,giant_module_sensitivity.parquet,region_summary.parquet,summary.json,DTU_ADDED_VALUE.md}` (`halves/` is a gitignored cache) | *(does the DTU evidence of a gene's module neighbours in one donor subset predict its evidence in an independent subset, given its own? Complementary to gene-wise DTU, not a test of it; giant-module table is supplementary)* |
| **Pooled cross-cohort replication** | `real_data/module_trust.py replication-pooled` | `_h/03e` | `_m/stability/module_trust/module_aging_replication_pooled__<method>{.parquet,__stats.json}` | *(null for both methods; reported, not a figure)* |
| **Phenotype-blind resolution sweep** | `real_data/stability.py sweep` | `_h/02a` (after `_h/01a` with `STABILITY_SAVE_EDGES=1`), then `_h/03a` | `_m/stability/stability_summary.parquet` (rows with method isograph_resXpY) | *(disclosed sensitivity; cannot select a resolution)* |
| Cross-cohort replication | `real_data/replication.py`, `replication_go.py` | `_h/01c`, `_h/02d` | `_m/replication/` | Fig 2 |
| Replication permutation null | `real_data/replication_permutation.py` | `_h/03d` (one array per `--covariates {full,complement,none}`), then `_h/04b` (report) | `_m/replication/`, `REPLICATION_PERMUTATION.md` | Fig 2 |
| Three-baseline comparison | `real_data/baseline_comparison.py` | `_h/03h` | `_m/baseline_comparison/` | S-real-1, S1/S2 |
| LR / software robustness | `real_data/stability.py` | `_h/retired/lr_{validation,validation_launch,aggregate}.sh` | — | — *(**retired 2026-09-14**: settled the single-LR promotion; outputs only at tag `legacy_expression_filter`)* |
| Giant-cap ablation | `real_data/stability.py` | `_h/retired/gcap_ab.sh` | `_m/stability_gcap_ab/` | — *(**retired 2026-09-12**: an A/B of a cap never promoted, not cited anywhere; note that production has since returned to resolution 2.0, which does not revive it)* |

## 05 — Genetic anchoring

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| xQTL anchoring (17 analyses) | `real_data/qtl_anchoring.py` | `05_genetic_anchoring/_h/01a–01b` | `<store>/qtl_anchoring.parquet` | Fig 3 |
| xQTL anchoring sensitivity (17 x 3 arms) | `real_data/qtl_anchoring.py --outcome/--covariate-set` | `_h/01c` (IsoGraph), `_h/01b` with the same flags (matched baselines) | `<store>/qtl_anchoring_{constraint,continuous,dose}.parquet` | Table S19 |
| Contrast meta-analysis | `real_data/qtl_anchoring_meta.py` | `_h/02a` | `_m/qtl_anchoring_meta/` | **Fig 3, Table 1**, S3–S5 |
| Sensitivity meta-analysis | `real_data/qtl_anchoring_meta.py --outcome/--covariate-set` | `_h/02a` with those flags | `_m/qtl_anchoring_meta/sensitivity/<arm>/` | Table S19 |
| sQTL direction concordance | `real_data/sqtl_concordance.py`, `sqtl_concordance_meta.py` | `_h/01d`, then `_h/02b` | `_m/sqtl_concordance_meta/` | — |
| Module genetic anchoring | `real_data/module_genetic_anchoring.py` | `_h/01e`, then `_h/02c` | `_m/module_genetic_anchoring_meta/` | Table S16 only — **no figure by design** (pooled-gene, not per-module) |
| MAGMA module-GWAS | `gwas/prepare_magma_inputs.py` | `_h/01f`, `_h/02d`, `_h/03a` (each also with `MAGMA_ISOGRAPH_BACKEND=isograph_vae_res2`) | `_m/gwas/magma_results_combined{,_res2}.parquet` | S-real-2 |
| Colocalization (SuSiE + eCAVIAR) | `real_data/coloc_{prep,summary,isoform_events,direction,meta}.py` | `_h/01g` (prep), `_h/02e` (LD), `_h/03b` (CLPP), `_h/04a` (meta), `_h/04b` (direction + events) | `_m/coloc/` | Fig 4, S8/S9 |
| Module-level coloc convergence (5 traits) | `real_data/module_coloc_convergence.py` | `_h/04c` | `_m/module_coloc_convergence/` | S-real-10, Tables S20a/b |
| **Deep-dive event table** | `real_data/gene_deep_dive.py --part events` | `_h/05c` | `_m/deep_dive/deep_dive_events.{parquet,tsv}` | *(read by stage 06; the panel is stage 08)* |
| **Signal-level coloc (`coloc.susie`)** | `real_data/coloc_signal_susie.py` + `_h/03c.coloc_gwas_susie.R`, `_h/06b.coloc_signal_susie.R` | `_h/03c` (GWAS cache), `_h/05b` (prep), `_h/06b`, `_h/07b` (meta) | `_m/coloc_signal_susie/` | *(estimator hierarchy: susie > abf > CLPP)* |
| Per-gene sQTL-vs-eQTL coloc contrast | `real_data/coloc_modality_contrast.py` + `_h/05a.coloc_modality_abf.R` | `_h/04d` (prep), `_h/05a`, `_h/06a` (meta, per arm), `_h/07a` (compare) | `_m/coloc_modality_contrast/` | *(run 2026-09-01, 4 arms — **null**; reported as a negative control, no figure)* |
| **Locus event audit** | `real_data/locus_event_audit.py` (`configs/known_splice_events.yaml`) | `_h/08b` | `_m/locus_event_audit/<nominations>/` | *(tiering of UNC13A / PICALM / SNCA)* |
| **Locus LD robustness audit** | `_h/08d.locus_ld_robustness.R` | `_h/08d` (manual, per locus) | `_m/locus_ld_robustness/<analysis>__<LOCUS_ID>/` | *(required before a recovered or headline locus carries a claim)* |
| **BrainSEQ switch-QTL (S_g) vs abundance-QTL (A_g)** | `real_data/brainseq_switch_qtl.py` | `_h/01i` (GPU mapping), `_h/02h` (meta) | `_m/brainseq_switch_qtl/{all_samples,ea_only}/{caudate,dlpfc,hippocampus}/qtl/cis_qtl_switch.parquet`, `_m/brainseq_switch_qtl/{all_samples,ea_only}/modality_contrast.parquet` | *(same-tissue genetic anchoring, not replication)* |
| **BrainSEQ QTL checks: S_g reproduction / sign pin, A_g positive control** | `real_data/brainseq_qtl_checks.py` | `_h/02g` | `_m/brainseq_switch_qtl/{all_samples,ea_only}/checks` | *(gate before any BrainSEQ QTL result is read)* |
| **BrainSEQ switch vs abundance effect sizes, selection-symmetric** | `real_data/brainseq_switch_qtl.py --stage effect_size` | `_h/03e` (array over arms) | `_m/brainseq_switch_qtl/{all_samples,ea_only}/{effect_size_symmetric.parquet,lead_variant_frequency.parquet,EFFECT_SIZE_SYMMETRIC.md}` | *(every gene tested on both axes; own-lead, cross-lead and common-variant contrasts; |slope| and conservative |z| both reported, with lead-variant MAF)* |
| **BrainSEQ signal-level coloc (S_g / A_g, EA-only)** | `real_data/coloc_brainseq.py` + `_h/07c.coloc_brainseq_susie.R` | `_h/06c` (prep), `_h/07c`, `_h/08c` (meta) | `_m/coloc_brainseq/ea_only/` | *(same-tissue genetic anchoring; in-sample QTL LD; gated on the `_h/02g` checks)* |
| **SMR + HEIDI** | `real_data/smr_heidi.py` | `_h/09a` runs prep and submits `_h/10a` (compute; sensitivity arms `SMR_PEQTL=1e-6` and `SMR_MULTI=1`, both with `SMR_SKIP_BESD=1`) and `_h/11a` (meta, `--peqtl-smr 1e-6` / `--smr-multi`); `SMR_QTL_SOURCE=brainseq,SMR_ARM=ea_only` for BrainSEQ | `_m/smr_heidi/{gtex,brainseq/ea_only}`, each with `sensitivity/{peqtl_smr_1e-06,smr_multi}` | *(orthogonal corroboration beneath coloc; BrainSEQ gated on the `_h/02g` checks and read against BrainSEQ's own coloc; two separately-corrected testing families, `smr_status` separates untested from tested-null, instrument F reported, attrition trace shows GWAS coverage — not allele QC — is the dominant filter)* |
| S-LDSC partitioned heritability | `real_data/ldsc_annot_prep.py`, `ldsc_summary.py` | `_h/01h`, `_h/02f`, `_h/03d`, then `_h/04e` (summary) | `_m/ldsc/ldsc_partitioned.parquet` | Fig 4B |

## 06 — Switch mechanism

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| Switch consequence + meta | `real_data/switch_consequence[_meta].py` | `06_switch_mechanism/_h/01a`, then `_h/02a` | `_m/switch_consequence_meta.parquet` | S-real-5 |
| PSI / junction validation | `real_data/validate_switch_splicing.py` | `_h/01b–01c` | `_m/switch_validation/` | — |
| ISA / satuRn concordance | `real_data/isa_concordance.py` | `_h/01d` | `_m/isa_concordance/` | **S-real-9** `figIsaConcordance`, Table S15 |
| Long-read confirmation | `real_data/longread_switch_confirm.py` | `_h/01e` | `_m/longread_switch_confirm/` | — |
| Short-read junction confirmation, **all concordant genes** | `real_data/junction_coloc_confirm.py` | `_h/01k` | `_m/junction_coloc_confirm/{targets,junction_confirm}.parquet`, `params.json`, `JUNCTION_COLOC_CONFIRM.md` | *(generalised 2026-09-20 from two hardcoded genes to all 30; genes in tissues BrainSEQ cannot sequence are reported as `no_matched_brainseq_region` rather than dropped; BH within the paired family. SNCA and CTSH are no longer in the target set)* |
| Orthogonal confirmation | `real_data/switch_orthogonal_confirm.py` | `_h/02b` (`--mode anchored`, `--mode global-null`) | `_m/switch_orthogonal_confirm/` | **S-real-8** `figOrthogonalConfirm`, Table S14 |
| **Orthogonal confirmation, signal-level coloc** | `real_data/coloc_isoform_events.py --layer signal` → `switch_orthogonal_confirm.py --events signal` | `05_genetic_anchoring/_h/08a` → `_h/02b --events signal` | `05_genetic_anchoring/_m/coloc_signal_susie/all_introns/coloc_isoform_events.parquet`, `_m/switch_orthogonal_confirm/signal_coloc/` | *(Analysis 6: long-read check over the coloc.susie nominations; the CLPP arm above is unchanged)* |
| Clinical consequence | `real_data/clinical_consequence[_meta].py` | `_h/01f` (download, login node), `_h/02c`, then `_h/03a` (meta) | `_m/clinical_consequence_meta.parquet` | S-real-7 |
| SCZ confound sensitivity | `real_data/scz_confound_sensitivity.py` | `_h/01g` | `_m/scz_confound_sensitivity/` | — |
| Feature sensitivity (5 axes, partition fixed) | `real_data/switch_feature_sensitivity.py` | `_h/01h` (BrainSEQ caudate), `_h/01i` (13 GTEx regions), then `_h/02d` (`--aggregate`: cross-region report + quantification axis) | `_m/switch_feature_sensitivity/{SWITCH_FEATURE_SENSITIVITY.md,sensitivity_summary_all.parquet,quantification_concordance.parquet}`; per region under `_m/switch_feature_sensitivity/brainseq/caudate` and `_m/switch_feature_sensitivity/gtex/*` | — |
| **Feature sensitivity, full refit per setting** | `real_data/switch_feature_refit.py` | `_h/01j` (array, one fit per setting), then `_h/02e` (aggregate) | `_m/switch_feature_sensitivity/refit/brainseq/caudate/{refit_summary.parquet,REFIT_SENSITIVITY.md}` (one subdirectory per setting) | *(published-setting refit is the noise floor)* |
| **phASER allelic feasibility screen** | `real_data/ase_switch_direction.py` | *(login)* | `_m/ase_switch_direction/<region>/` | *(gate: is an allele-specific switch test supported at all?)* |
| **Allele-aware junction recount (PI 10a), Quest only** | `real_data/ase_junction_switch.py` (`targets`, `count`, `concordance`, `screen`), `real_data/ase_junction_allelic.py` (`test`, step 4) | `_h/01l`, `_h/02f` (array, one WASP BAM per task), `_h/03b`, `_h/04a`; chained by `_h/ase_junction_quest.sh`, manual steps in the Bridges DAG | `_m/ase_junction_switch/<region>/{junction_allelic_counts,pair_feasibility,allelic_test,allelic_donor_counts}.parquet`, `ASE_JUNCTION_SCREEN.md`, `ASE_JUNCTION_ALLELIC.md` (per-BAM `counts/` gitignored) | *(the two-sided rescue of the screen above: same pre-registered gate, on the gate family of switch-QTL genes; the arm proceeds only if >= 30 pairs pass)* |
| **Risk-allele orientation of the allelic test (PI 10a, step 5), Bridges-2** | `real_data/ase_risk_orientation.py` | `_h/05a.ase_risk_orientation.sh` | `_m/ase_junction_switch/<region>/{risk_orientation.parquet,risk_orientation_summary.json,ASE_RISK_ORIENTATION.md}` | *(re-signs `beta` from the QTL lead's ALT allele to the GWAS risk allele of the colocalizing locus, via signed LD in the BrainSEQ panel; gated at |r| >= 0.8)* |
| **Allelic test anchored on the GWAS lead (PI 10a, second arm), Bridges-2** | `real_data/ase_gwas_lead_arm.py` | `_h/05b.ase_gwas_lead_arm.sh` | `_m/ase_junction_switch/<region>/{gwas_lead_arm.parquet,gwas_lead_arm_summary.json,GWAS_LEAD_ARM.md}` | *(fits the same within-donor contrast AT the GWAS lead instead of transferring the sign by LD, so orientation is exact and the cost is power. No recount: it relabels the haplotype counts. Anchor genotypes come from phASER's own per-sample VCFs -- the TOPMed panel is phased in a different frame (agreement 0.51) and must not be used -- and the run self-checks that before fitting)* |
| **Module-level cis control of isoform choice (06a P2 item 2), Bridges-2** | `real_data/ase_junction_allelic.py` (`module`) | `_h/05c.ase_module_cis_control.sh` | `_m/ase_junction_switch/{<region>/,}module_cis_control{,_min5}.parquet`, `module_cis_control_summary{,_min5}.json`, `MODULE_CIS_CONTROL{,_MIN5}.md` | *(aggregates the fitted allelic test to the co-switching module and asks whether the modules with proportionally more cis-controlled genes carry the stronger age association. Refits nothing and maps no module eigengene. Age strength is ranked within (region, trait) because linear and spline effects are not comparable; the null permutes gene -> module labels, holding module sizes and the cis-gene count fixed. **Null**: caudate rho -0.691 (p = 0.070, opposite sign), hippocampus +0.326 -> +0.018 under the >= 5-gene floor, DLPFC not estimable; pooled -0.183 (p = 0.47))* |
| **PSI / junction validation — interpretation** | *(written from `summary_*.json`)* | — | `_m/switch_validation/SWITCH_VALIDATION_SUMMARY.md`, `_m/longread_switch_confirm/LONGREAD_SWITCH_CONFIRM.md` | — |

## 07 — RBP regulation

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| Motif families | `real_data/rbp_motif_families.py` | `07_rbp_regulation/_h/01a` | `_m/rbp/rbp_motif_families.parquet` | S-real-6 |
| Regulons (mature + intronic) | `real_data/rbp_scan[_intronic].py`, `rbp_regulon.py` | `_h/02a`, `_h/03a` | `_m/rbp/rbp_regulon*.parquet` | S-real-6, S10 |
| CLIP binding evidence | `real_data/rbp_binding.py` | `_h/03b` (fetch, login node), `_h/04a` | `_m/rbp/rbp_binding*.parquet` | S-real-6 panel B, Table S17 |
| Neuronal CLIP | `real_data/neuronal_clip_*.py` | `_h/04b`, `_h/05a`, `_h/06a`, `_h/07a` | `_m/neuronal_clip/` | — |
| NOVA family / NOVA2 | `real_data/nova_family_renomination.py`, `nova2_*.py` | `_h/07b`, `_h/08a`, `_h/09a` | `_m/neuronal_clip/nova*/` | — |

## 08 — Integration

| Analysis | CLI | Wrapper | Outputs | Display |
|---|---|---|---|---|
| Per-gene deep dive | `real_data/gene_deep_dive.py --part panel` | `08_integration/_h/01a` | `_m/deep_dive/` | **Table 2**, S8–S12 |
| SCZ age projection | `real_data/scz_age_projection.py` | `_h/01b` | `_m/scz_age_projection/` | **Fig 4E** (folded from `figSczConvergence`), Table S18 |
| Functional preservation of matched modules | `real_data/replication_functional.py` | `_h/01c` | `_m/functional_preservation/` | Fig 2 |
| Target panel / assayability | `real_data/rbp_target_panel.py`, `rbp_pair_assayability.py` | `_h/02a` (per RBP arm), `_h/03a` | `_m/rbp_target_panel/` | — |
| **Anchored gene summary (genetics + orthogonal validation)** | `real_data/anchored_gene_summary.py` | `_h/04a` | `_m/anchored_gene_summary/{anchored_gene_summary.{parquet,tsv},ANCHORED_GENE_SUMMARY.md,provenance.json}` | *(one table for the 30 concordant genes: coloc.abf as primary, SMR/HEIDI, eCAVIAR CLPP as a trailing secondary column, beside long-read / junction-recount / PSI validation. Every number generated; `provenance.json` carries the commit, dirty flag, thresholds and per-input digests. Re-run it whenever stage 05 or 06 is re-run)* |

`<store>` = `02_module_discovery/<cohort>/<region>/_m/`, the shared cohort × region
artifact store.
