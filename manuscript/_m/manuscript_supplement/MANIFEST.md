# Manuscript supplementary items

Written by `manuscript/_h/build_supplementary_data.py` on 2026-10-04.
Supplementary Tables stay in the PDF (<= 40 rows and <= 12 columns); Supplementary Data are separate Excel workbooks (README + Data sheets) with CSV copies.

| Item | Manuscript file | Source CSV | Rows | Columns | Title |
| --- | --- | --- | --- | --- | --- |
| Table S1 | `supplementary_tables/tableS01_benchmark_summary.csv` | `01_synthetic_benchmark/03_metrics/_m/tableS_benchmark_summary.csv` | 36 | 12 | Synthetic benchmark performance across scenarios and methods. |
| Table S2 | `supplementary_tables/tableS02_scale_compute_summary.csv` | `01_synthetic_benchmark/03_metrics/_m/tableS_scale_compute_summary.csv` | 15 | 7 | Runtime and memory use across network sizes and compute settings. |
| Table S3 | `supplementary_tables/tableS03_axis_orthogonality.csv` | `manuscript/_m/supp_tables/tableS13a_axis_orthogonality.csv` | 17 | 12 | Separation of the switch and abundance coordinates, per analysis. |
| Table S4 | `supplementary_tables/tableS04_baseline_pooled.csv` | `manuscript/_m/supp_tables/tableS1_baseline_pooled.csv` | 4 | 10 | Pooled per-module association rates across the 17 applications. |
| Table S5 | `supplementary_tables/tableS05_clinical_consequence.csv` | `manuscript/_m/supp_tables/tableS37_clinical_consequence.csv` | 6 | 11 | Constraint and ClinVar density of switch genes and switched exons. |
| Table S6 | `supplementary_tables/tableS06_qtl_specificity_contrast.csv` | `manuscript/_m/supp_tables/tableS3_qtl_specificity_contrast.csv` | 10 | 8 | Splicing-QTL specificity contrast, all analyses. |
| Table S7 | `supplementary_tables/tableS07_qtl_specificity_matched_baseline.csv` | `manuscript/_m/supp_tables/tableS4_qtl_specificity_matched_baseline.csv` | 10 | 8 | Splicing-QTL specificity contrast on the analyses shared by all methods. |
| Table S8 | `supplementary_tables/tableS08_module_cis_control.csv` | `manuscript/_m/supp_tables/tableS23_module_cis_control.csv` | 24 | 9 | Module-level cis control of isoform choice. |
| Table S9 | `supplementary_tables/tableS09_synthetic_scenarios.csv` | `manuscript/_m/supp_tables/tableS35_synthetic_scenarios.csv` | 16 | 8 | Synthetic benchmark scenarios. |
| Data S1 | `supplementary_data/dataS01_cohort_description.xlsx` | `manuscript/_m/supp_tables/tableS34_cohort_description.csv` | 17 | 22 | Discovery cohorts. |
| Data S2 | `supplementary_data/dataS02_composition_adjustment.xlsx` | `manuscript/_m/supp_tables/tableS13_composition_adjustment.csv` | 12 | 14 | Switch-unique genes before and after adjustment for estimated cellular composition. |
| Data S3 | `supplementary_data/dataS03_module_trust_funnel.xlsx` | `manuscript/_m/supp_tables/tableS7_module_trust_funnel.csv` | 12 | 16 | Per-region module trust funnel. |
| Data S4 | `supplementary_data/dataS04_baseline_per_region.xlsx` | `manuscript/_m/supp_tables/tableS2_baseline_per_region.csv` | 68 | 12 | Per-analysis module association rates for all four module sets. |
| Data S5 | `supplementary_data/dataS05_split_half_module_ledger.xlsx` | `manuscript/_m/supp_tables/tableS7a_split_half_module_ledger.csv` | 190 | 12 | Split-half module ledger. |
| Data S6 | `supplementary_data/dataS06_resolution_sensitivity.xlsx` | `manuscript/_m/supp_tables/tableS7e_resolution_sensitivity.csv` | 54 | 13 | Split-half agreement across the Leiden resolution sweep. |
| Data S7 | `supplementary_data/dataS07_projection_module_ledger.xlsx` | `manuscript/_m/supp_tables/tableS7b_projection_module_ledger.csv` | 155 | 20 | Cross-cohort eigengene projection ledger. |
| Data S8 | `supplementary_data/dataS08_crosscohort_permutation.xlsx` | `manuscript/_m/supp_tables/tableS7c_crosscohort_permutation.csv` | 36 | 14 | Cross-cohort matched-pair count against its permutation nulls. |
| Data S9 | `supplementary_data/dataS09_functional_preservation.xlsx` | `manuscript/_m/supp_tables/tableS7d_functional_preservation.csv` | 12 | 15 | Functional preservation of matched cross-cohort module pairs. |
| Data S10 | `supplementary_data/dataS10_saturn_concordance.xlsx` | `manuscript/_m/supp_tables/tableS15_isa_concordance.csv` | 17 | 14 | satuRn concordance of IsoGraph switch genes, per analysis. |
| Data S11 | `supplementary_data/dataS11_psi_gene_corroboration.xlsx` | `manuscript/_m/supp_tables/tableS38_psi_gene_corroboration.csv` | 19 | 18 | Gene-level corroboration of the switch layer by junction usage. |
| Data S12 | `supplementary_data/dataS12_module_context_heldout_dtu.xlsx` | `manuscript/_m/supp_tables/tableS21_module_context_heldout_dtu.csv` | 6 | 31 | Module context and held-out DTU evidence, per analysis. |
| Data S13 | `supplementary_data/dataS13_rbp_regulons.xlsx` | `manuscript/_m/supp_tables/tableS32_rbp_regulons.csv` | 11040 | 24 | RNA-binding-protein motif regulons of the switch modules. |
| Data S14 | `supplementary_data/dataS14_rbp_eclip_binding.xlsx` | `manuscript/_m/supp_tables/tableS33_rbp_eclip_binding.csv` | 35 | 19 | ENCODE eCLIP binding at switched and constitutive exons. |
| Data S15 | `supplementary_data/dataS15_allelic_imbalance_regions.xlsx` | `manuscript/_m/supp_tables/tableS22_allelic_imbalance_regions.csv` | 4 | 27 | Within-donor allelic test of isoform choice, per region. |
| Data S16 | `supplementary_data/dataS16_signal_coloc_nominations.xlsx` | `manuscript/_m/supp_tables/tableS25_signal_coloc_nominations.csv` | 127 | 16 | Signal-level sQTL colocalization nominations. |
| Data S17 | `supplementary_data/dataS17_coloc_modality_contrast.xlsx` | `manuscript/_m/supp_tables/tableS26_coloc_modality_contrast.csv` | 7 | 18 | Paired sQTL and eQTL colocalization contrast. |
| Data S18 | `supplementary_data/dataS18_brainseq_axis_coloc_contrast.xlsx` | `manuscript/_m/supp_tables/tableS31_brainseq_axis_coloc_contrast.csv` | 7 | 18 | BrainSEQ switch- and abundance-axis colocalization contrast. |
| Data S19 | `supplementary_data/dataS19_coloc_isoform_events.xlsx` | `manuscript/_m/supp_tables/tableS28_coloc_isoform_events.csv` | 410 | 24 | Colocalization events and their mapping to IsoGraph switch pairs. |
| Data S20 | `supplementary_data/dataS20_smr_heidi.xlsx` | `manuscript/_m/supp_tables/tableS30_smr_heidi.csv` | 820 | 20 | SMR and HEIDI for the colocalization nominations. |
| Data S21 | `supplementary_data/dataS21_clpp_isoform_events.xlsx` | `manuscript/_m/supp_tables/tableS24_clpp_isoform_events.csv` | 467 | 21 | CLPP-nominated isoform events and their mapping to IsoGraph switch pairs. |
| Data S22 | `supplementary_data/dataS22_ldsc_partitioned.xlsx` | `manuscript/_m/supp_tables/tableS27_ldsc_partitioned.csv` | 30 | 14 | Partitioned heritability of the switch-derived QTL annotations. |
| Data S23 | `supplementary_data/dataS23_coloc_convergence.xlsx` | `manuscript/_m/supp_tables/tableS20a_coloc_convergence_global.csv` | 10 | 14 | Module concentration of colocalizing switch genes. |
| Data S24 | `supplementary_data/dataS24_magma_module_gwas.xlsx` | `manuscript/_m/supp_tables/tableS36_magma_module_gwas.csv` | 8064 | 13 | MAGMA competitive gene-set tests of every module. |
| Data S25 | `supplementary_data/dataS25_longread_coloc_confirmation.xlsx` | `manuscript/_m/supp_tables/tableS29_longread_coloc_confirmation.csv` | 55 | 18 | Long-read confirmation of colocalization-anchored switch pairs. |
| Data S26 | `supplementary_data/dataS26_junction_pair_corroboration.xlsx` | `manuscript/_m/supp_tables/tableS40_junction_pair_corroboration.csv` | 260 | 30 | Short-read junction corroboration of colocalization-prioritized transcript pairs. |
| Data S27 | `supplementary_data/dataS27_psi_junction_confirmation.xlsx` | `manuscript/_m/supp_tables/tableS39_psi_junction_confirmation.csv` | 173 | 20 | Short-read junction confirmation of the CLPP-anchored switch pairs. |
| Data S28 | `supplementary_data/dataS28_projection_sign_scale.xlsx` | `manuscript/_m/supp_tables/tableS7f_projection_sign_scale.csv` | 12 | 17 | Raw against null-standardized projected-age sign agreement, per region. |
