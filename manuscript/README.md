# Manuscript

Every display item in the paper, its builder, and the analysis outputs it reads.
Numbers are never hand-edited here: each item regenerates from committed ledgers.

- `_h/` — figure and table builders (one per display item).
- `_m/figures/` — all real-data figure panels (PDF + PNG).
- `_m/main_tables/`, `_m/supp_tables/` — CSV + Markdown tables.
- [`FIGURE_ORDERING.md`](FIGURE_ORDERING.md) — the ordering, with the one honest
  claim each item carries. **Authoritative for figure numbering.**
- [`MANUSCRIPT_PLAN.md`](MANUSCRIPT_PLAN.md) — evidence inventory, claims hierarchy,
  Results/Methods outlines, journal requirements.
- [`GENETIC_ANCHORING_RESULTS.md`](GENETIC_ANCHORING_RESULTS.md) — drafted Results
  prose for the anchoring section, with Manubot citekeys.
- [`CELL_GENOMICS_IMPACT_PROGRESS.md`](CELL_GENOMICS_IMPACT_PROGRESS.md),
  [`PER_GENE_DEEP_DIVE_PLAN.md`](PER_GENE_DEEP_DIVE_PLAN.md) — planning notes.

## Main display items

| Item | Builder | Reads from |
| --- | --- | --- |
| Fig 1 `figConceptOverview` | `_h/concept_overview_figure.R` | `01_synthetic_benchmark/01_synthetic/_m/synthetic_results.parquet` (core scenarios for D–E; the `genetic_anchoring` scenario for F–G) |
| Fig 2 `figTrustFunnel` | `_h/trust_funnel_figure.R` | `04_module_trust/_m/stability/module_trust/`, `04_module_trust/_m/stability/eigengene_projection/eigengene_projection_summary.parquet` (panel A) |
| Fig 3 `figQtlSpecificity` | `_h/qtl_specificity_figure.R` | `05_genetic_anchoring/_m/qtl_anchoring_meta/` |
| Fig 4 `figGeneticAnchoring` | `_h/genetic_anchoring_figure.R` | `05_genetic_anchoring/_m/{deep_dive,ldsc,scz_age_projection}/` |
| Fig 5 `figCompositionRobustness` | `_h/composition_robustness_figure.R` | `03_module_characterization/_m/`, `02_module_discovery/gtex/_m/composition/` |
| Table 1 `table2_qtl_specificity_contrast` | `_h/assemble_main_tables.py` | `05_genetic_anchoring/_m/qtl_anchoring_meta/` |
| Table 2 `table3_splicing_led_genes` | `_h/assemble_main_tables.py` | `08_integration/_m/deep_dive/` |

The synthetic supplements (S1–S12) are built **in place** inside
`01_synthetic_benchmark/03_metrics/figures/` by
`isograph_benchmark/figures/synthetic_benchmark.R` and are not copied here. That script
also still builds `fig1_benchmark_overview`, the 24-panel grid Fig 1 replaced; it is kept
as supplementary material, not as Fig 1.

## Supplementary real-data figures

| Item | Builder | Reads from |
| --- | --- | --- |
| S-real-1 `figBaselineRates` | `_h/baseline_rates_figure.R` | `04_module_trust/_m/baseline_comparison/` |
| S-real-2 `figGwasResolution` | `_h/gwas_resolution_figure.R` | `05_genetic_anchoring/_m/gwas/` |
| S-real-3 `figGoInvisible` | `_h/go_invisible_figure.R` | `02_module_discovery/brainseq/caudate_sczd/_m/` |
| S-real-4 `figSeparation` | `_h/abundance_structure_figure.R` | `03_module_characterization/_m/axis_orthogonality_{all,summary}.parquet` (A, all 17 analyses); `02_module_discovery/brainseq/caudate_sczd/_m/isograph_vae/abundance_structure/` (B-C) |
| S-real-5 `figSwitchConsequence` | `_h/switch_consequence_figure.R` | `06_switch_mechanism/_m/` |
| S-real-6 `figRbpRegulon` | `_h/rbp_regulon_figure.R` | `07_rbp_regulation/_m/rbp/` |
| S-real-7 `figClinicalConsequence` | `_h/clinical_consequence_figure.R` | `06_switch_mechanism/_m/`, `08_integration/_m/deep_dive/` |
| S-real-8 `figOrthogonalConfirm` | `_h/orthogonal_confirm_figure.R` | `06_switch_mechanism/_m/switch_orthogonal_confirm/` (`anchored_summary.json`, `global_null_summary.json`, `anchored_gene_confirmation.parquet`) |
| S-real-9 `figIsaConcordance` | `_h/isa_concordance_figure.R` | `06_switch_mechanism/_m/isa_concordance/` |

`_h/scz_convergence_figure.R` builds `figSczConvergence`, which is **folded into Fig 4E**
and is no longer a numbered display item on its own. The builder is kept because it is
the reproducible standalone form of the same result; its per-module candidate RBP
regulators, which do not fit at half width, are in Table S18.

Supplementary tables S1–S18 are built by `_h/assemble_supp_tables.py`; see
`_m/supp_tables/SUPPLEMENTARY_TABLES.md`.

## Rebuilding

Login node, no SLURM. Figures use the `rnaseq` R environment; tables use the
`isograph` python environment.

```bash
/ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript manuscript/_h/qtl_specificity_figure.R
/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python manuscript/_h/assemble_main_tables.py
/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python manuscript/_h/assemble_supp_tables.py
```
