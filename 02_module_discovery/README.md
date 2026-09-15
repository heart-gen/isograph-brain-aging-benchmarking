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
`wgcna_multiplex/`, `isograph_vae_res2/`) and the per-region tables the
downstream stages write into it (`module_enrichment/`, `qtl_anchoring.parquet`, …).
Address it in code with `region_store()` / `region_artifact_dir()`
from `isograph_benchmark.paths` — never by literal path.

**Transcript filter (2026-09-14).** Every IsoGraph fit and matched baseline uses the switching
filter (`run_models.filter_production_transcripts`): gene count ≥ 10 in ≥ 70% of samples;
transcript count ≥ 10 and share of its gene ≥ 0.10, each in ≥ 10% of samples. The pipeline and
results before the switch — BrainSEQ aging on count > 10 in ≥ 70%, SCZD and GTEx unfiltered — are
frozen at tag `legacy_expression_filter`. Retired with that switch and absent from the live tree:
the with-abundance arm (steps 04, 06), refit/reprojection QC and tier projection (steps 12–14,
stage 04 steps 12–13), and the transcript-filter arms (step 17), which chose this filter.

Four heavy artifacts per fit (`feature_scores`, `feature_reconstruction`,
`high_vs_low_table`, `edges`) are gitignored and distributed via Zenodo; see `zenodo/`.

## Order

| Step | Wrapper | Produces |
|---|---|---|
| 01–03 | `run_isograph_{brainseq_aging,brainseq_sczd,gtex}` | IsoGraph fits: `modules`, `edges`, `traits`, `age_{linear,spline}`, `feature_scores`, `module_gene_roles`, `calibration` |
| 04 | `run_isograph_brainseq_with_abundance` | *Retired 2026-09-14* — abundance-channel variant at the pre-5.0 resolutions |
| 05–06 | `sweep_leiden_brainseq[_with_abundance]` | Resolution sweep (re-clusters saved edges, no refit); step 06 retired with step 04 |
| 07–09 | `wgcna_gene_{brainseq_aging,brainseq_sczd,gtex}` | Classical gene-level WGCNA baseline |
| 10–11 | `wgcna_matched_features_{brainseq,gtex}` | `wgcna_switch_only` + `wgcna_multiplex` baselines on **identical** features — the primary internal control |
| 12–14 | `refit_qc_{brainseq,gtex}`, `reproject_qc_gtex` | *Retired 2026-09-14* — refit QC and reprojection checks |
| 15 | `sweep_leiden_gtex` | Resolution sweep on the 13 GTEx fits — a disclosed sensitivity; the CLI refuses `--write-best` for GTEx |

Three WGCNA baselines exist on purpose: `wgcna_gene` (classical abundance),
`wgcna_switch_only` and `wgcna_multiplex` hold the switch features fixed and vary only
the network inference, which is what isolates a method effect from a feature effect.

## Granularity: resolution and edge thresholds

**Resolution.** Canonical Leiden resolution 5.0, stated on the ≥ 900-gene giant-module
criterion: at 2.0, 17 of 28 significant module–GWAS hits were modules of ≥ 900 genes; at 5.0
none are (`05_genetic_anchoring/_m/gwas/GWAS_RESOLUTION_SUMMARY.md`). The size criterion is
phenotype-blind; the GWAS result is its confirmation. **Decided 2026-09-12 (PI): keep this as
the stated basis.** The phenotype-blind split-half sweep (`03_module_trust/_h/16`) was run and
cannot select a resolution — ARI is U-shaped and NMI monotone, both tracking granularity — so
it is reported as a disclosed sensitivity with that reason, alongside the BrainSEQ (`_h/05`)
and GTEx (`_h/15`) module-count sweeps.

**Edge thresholds.** The gene graph Leiden clusters is built from VAE feature similarity with
two thresholds, and neither was tuned against a trait:

| Threshold | Edges it admits | How it is set | Value in production |
|---|---|---|---|
| `alpha_switch` | switch–switch | Fixed default; no grid is passed, so `select_alpha_switch` never runs | **0.50** in all 17 fits |
| `alpha_abundance` | abundance–abundance | Selected per fit from the grid 0.70–0.95 by `isograph.models.multiplex.select_alpha_abundance`: the *smallest* value whose abundance edges merge no two switch-defined components of the graph built without them. Falls back to 0.95 when every candidate merges | 0.95 in all 4 BrainSEQ fits and 6 GTEx; 0.90 in 5 GTEx; 0.85 in frontal cortex BA9 and putamen |

So the abundance threshold has a stated, structural, trait-blind criterion — abundance may
connect genes but may not fuse switch modules — recorded per fit as
`selected_alpha_abundance` in `calibration.parquet`. `alpha_switch` is a fixed default, and
**its sensitivity has not been measured**; the resulting giant component is 1.8–5.9% of
assigned genes in every fit, so it is not producing degenerate partitions, but that is a
sanity bound, not an insensitivity result.

**CLIs:** `isograph_benchmark/real_data/{run_models,run_matched_wgcna,sweep_leiden}.py`.

## Display items

Feeds every figure; directly backs S-real-4 (`figSeparation`) via the
`caudate_sczd` abundance-structure outputs.
