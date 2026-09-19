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
`03_module_characterization/README.md`.

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
frozen at tag `legacy_expression_filter`. Retired with that switch, absent from the live tree and
now under `_h/retired/`: the with-abundance arm, refit/reprojection QC, the transcript-filter arms
that chose this filter, and (in stage 03) tier projection.

Four heavy artifacts per fit (`feature_scores`, `feature_reconstruction`,
`high_vs_low_table`, `edges`) are gitignored and distributed via Zenodo; see `zenodo/`.

## Order

Run the whole stage with `bash 02_module_discovery/_h/run_stage.sh` (add `--dry-run` to print
the plan). The leading number of a wrapper is its tier; steps in one tier run in parallel.

| Step | Wrapper | Waits on | Produces |
|---|---|---|---|
| 00a | `production_gene_universe` | bundles | `production_gene_universe.{parquet,json}` per region — the genes surviving the production transcript filter. The gene-level WGCNA baseline is fit in R and cannot call that filter, so it reads this |
| 01a–01c | `run_isograph_{brainseq_aging,brainseq_sczd,gtex}` | bundles | IsoGraph fits: `modules`, `edges`, `traits`, `age_{linear,spline}`, `feature_scores`, `module_gene_roles`, `calibration`. The same wrappers with `--leiden-resolution 2.0` write `isograph_vae_res2`, the disclosed resolution comparison |
| 01d–01f | `wgcna_gene_{brainseq_aging,brainseq_sczd,gtex}` | 00a | Classical gene-level WGCNA baseline, on the 00a gene universe |
| 01g–01i | `wgcna_matched_features_{brainseq,sczd,gtex}` | bundles | `wgcna_switch_only` + `wgcna_multiplex` baselines on **identical** features — the primary internal control |
| 02a–02b | `sweep_leiden_{brainseq,gtex}` | 01a–01c | Resolution sweep (re-clusters saved edges, no refit); GTEx is a disclosed sensitivity — the CLI refuses `--write-best` for GTEx |
| 02c | `module_sizes` | 01a–01i | Module-size tables for every method in every store (`_m/module_sizes/`) |

Retired (`_h/retired/`; outputs only at tag `legacy_expression_filter`):
`run_isograph_brainseq_with_abundance` and `sweep_leiden_brainseq_with_abundance` (the
abundance-channel variant at the pre-5.0 resolutions), `refit_qc_{brainseq,gtex}` and
`reproject_qc_gtex` (refit and reprojection checks), `transcript_filter_arms`, and the
unreferenced `aging_models.R`.

Three WGCNA baselines exist on purpose: `wgcna_gene` (classical abundance),
`wgcna_switch_only` and `wgcna_multiplex` hold the switch features fixed and vary only
the network inference, which is what isolates a method effect from a feature effect.

All three are fit on the same gene universe as IsoGraph. The matched baselines get it for
free — `run_matched_wgcna` applies `filter_production_transcripts` itself. `wgcna_gene` is
fit in R, so it reads the universe 00a writes; before 2026-09-15 it took the bundle's full
gene list instead, which is wider (BrainSEQ caudate 20,365 vs 19,309) and made every
IsoGraph-vs-WGCNA head-to-head — MAGMA GSA, module trust, replication — compare module sets
drawn from different gene pools. A missing universe file stops the fit; it never falls back.

## Granularity: resolution and edge thresholds

**Resolution.** Canonical Leiden resolution **2.0**. **Decided 2026-09-16 (PI), reversing the
2026-09-12 decision to keep 5.0.** Resolution 5.0 had been stated on the ≥ 900-gene
giant-module criterion — at 2.0, 17 of 28 significant module–GWAS hits were modules of
≥ 900 genes, and at 5.0 none were. That result predates both the PGC3 SCZ swap and the
switching transcript filter, and it reverses under the current data: of the FDR-significant
MAGMA hits, **20 of 26 are ≥ 900-gene modules at 5.0 versus 24 of 37 at 2.0**
(`05_genetic_anchoring/_m/gwas/magma_results_combined{,_res5}.parquet`), so 5.0 is now the
more giant-dominated of the two. Resolution 2.0 also assigns **38 % more genes to modules**
(88,318 → 121,521 gene-region assignments across the 17 regions); the cost is a larger
largest module (2,215 → 3,212 genes). The retired 5.0 fits remain beside the canonical ones
as `isograph_vae_res5` and are carried as the disclosed resolution comparison. The
phenotype-blind split-half sweep (`04_module_trust/_h/02a`) was run and cannot select a
resolution — ARI is U-shaped and NMI monotone, both tracking granularity — so it is reported
as a disclosed sensitivity with that reason, alongside the BrainSEQ (`_h/02a`) and GTEx
(`_h/02b`) module-count sweeps.

Note that resolution is **not** what sets module coverage: only 25–38 % of the gene universe
is assigned at 5.0, and the dominant cause is the `alpha_switch = 0.5` edge threshold, which
leaves 39–68 % of genes with no edge at all (`node_diagnostics.parquet`, `fate` column).
`min_module_size = 20` is the only channel through which resolution affects coverage.

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

**CLIs:** `isograph_benchmark/real_data/{run_models,run_matched_wgcna,sweep_leiden,module_sizes}.py`.

## Display items

Feeds every figure; directly backs S-real-4 (`figSeparation`) via the
`caudate_sczd` abundance-structure outputs.
