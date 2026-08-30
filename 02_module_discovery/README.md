# 02 — Module discovery

Fits IsoGraph and the matched WGCNA baselines on both cohorts and chooses module
granularity. Everything downstream reads the artifacts this stage writes.

**Question:** what co-switching modules exist in postmortem human brain?

## Cohorts

| Cohort | Regions | Trait | Notes |
|---|---|---|---|
| **BrainSEQ** | `caudate` (Phase 3), `hippocampus`, `dlpfc` (Phase 2) | Age | Controls only, adults (Age ≥ 18) |
| **BrainSEQ SCZD** | `caudate_sczd` | Dx | Control + schizophrenia |
| **GTEx v11** | 13 brain regions (`amygdala` … `substantia_nigra`) | AGE | Exact age (v8-preferred) |

Covariates come from `configs/real_data.yaml`: BrainSEQ adjusts for Sex, MoD, RIN,
mapping rate, mito rate and SNP PC1–PC5; GTEx for SEX, SMRIN, SMTSISCH, SMMAPRT. Age is
modelled linearly (eigengene–age correlation) and with a natural spline (df = 4), the
spline compared against the linear fit. Canonical Leiden resolution **5.0**, giant-cap
**off**, seed 13. Residualization is discovery-only — see
`04_module_characterization/README.md`.

## The artifact store

`{brainseq,gtex}/<region>/_m/` is the cohort × region artifact store and is shared by
every later stage. It is keyed by cohort × region rather than by stage, so it holds both
the fits written here (`isograph_vae/`, `wgcna_gene/`, `wgcna_switch_only/`,
`wgcna_multiplex/`, `isograph_vae_with_abundance/`) and the per-region tables the
downstream stages write into it (`module_enrichment/`, `qtl_anchoring.parquet`,
`tier_checks/`, …). Address it in code with `region_store()` / `region_artifact_dir()`
from `isograph_benchmark.paths` — never by literal path.

Four heavy artifacts per fit (`feature_scores`, `feature_reconstruction`,
`high_vs_low_table`, `edges`) are gitignored and distributed via Zenodo; see `zenodo/`.

## Order

| Step | Wrapper | Produces |
|---|---|---|
| 01–03 | `run_isograph_{brainseq_aging,brainseq_sczd,gtex}` | IsoGraph fits: `modules`, `edges`, `traits`, `age_{linear,spline}`, `feature_scores`, `module_gene_roles`, `calibration` |
| 04 | `run_isograph_brainseq_with_abundance` | Abundance-channel variant fit |
| 05–06 | `sweep_leiden_brainseq[_with_abundance]` | Resolution sweep (re-clusters saved edges, no refit) |
| 07–09 | `wgcna_gene_{brainseq_aging,brainseq_sczd,gtex}` | Classical gene-level WGCNA baseline |
| 10–11 | `wgcna_matched_features_{brainseq,gtex}` | `wgcna_switch_only` + `wgcna_multiplex` baselines on **identical** features — the primary internal control |
| 12–14 | `refit_qc_{brainseq,gtex}`, `reproject_qc_gtex` | Refit QC and reprojection checks |

Three WGCNA baselines exist on purpose: `wgcna_gene` (classical abundance),
`wgcna_switch_only` and `wgcna_multiplex` hold the switch features fixed and vary only
the network inference, which is what isolates a method effect from a feature effect.

**CLIs:** `isograph_benchmark/real_data/{run_models,run_matched_wgcna,sweep_leiden}.py`.

## Display items

Feeds every figure; directly backs S-real-4 (`figSeparation`) via the
`caudate_sczd` abundance-structure outputs.
