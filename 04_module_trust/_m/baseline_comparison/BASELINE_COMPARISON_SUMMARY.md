# Three-baseline module comparison (features vs method; honest scope)

Modular analysis summary for Manubot integration. Generated from
`04_module_characterization/_m/baseline_comparison/baseline_comparison{,_pooled}.parquet`. Every numeric
claim is reproduced from those tables; do not edit the numbers by hand — regenerate.

## Purpose

State, honestly, what IsoGraph does and does not beat WGCNA at on module-level metrics, and
separate two confounded explanations for any difference: the *input representation* (does
the method see isoform switching?) from the *network inference* (VAE + Leiden vs classical
WGCNA). This is the scope-setting analysis that keeps every other real-data claim honest:
it establishes that IsoGraph is **not** globally superior on module enrichment/phenotype
rates, so the defensible value must be DTU-without-DGE *content*, not better rates.

## Inputs

- **Four module sets per analysis**, all from the same samples:
  - `isograph` — VAE + Leiden on switch+abundance features (the method).
  - `wgcna_switch_only` — classical WGCNA on the SAME switch-only features.
  - `wgcna_multiplex` — classical WGCNA on the SAME switch+abundance features.
  - `wgcna_gene` — classical WGCNA on abundance only (the input control).
- **Module enrichment outputs** per method: per-module phenotype association (`pheno_fdr`),
  GO enrichment (`n_go_terms`), module size and count.
- 17 analyses (1 SCZD + 3 BrainSEQ aging + 13 GTEx aging); 16 have all four methods.

The two matched WGCNA baselines consume the identical feature matrix as IsoGraph, so the
isograph-vs-matched contrast isolates inference from representation, while
isograph/wgcna_switch_only-vs-wgcna_gene isolates representation.

## Methods text

For each analysis and method we computed per-module rates rather than totals, because module
totals scale with module count and IsoGraph runs at a finer resolution (median 35 modules,
median size 71) than the WGCNA baselines (median 8–18.5 modules, size 127–387), making raw
counts non-comparable. A module was phenotype-significant at `pheno_fdr ≤ 0.1`, GO-enriched
at `n_go_terms > 0`, and "both" if it was phenotype-significant and GO-enriched. Per-module
fractions were averaged across analyses to give pooled rates per method. The phenotype rate
contrasts switch-fed (isograph, wgcna_switch_only) against abundance-fed (wgcna_gene,
wgcna_multiplex) methods to test where phenotype sensitivity comes from; the
isograph-vs-wgcna_multiplex contrast on identical multiplex features isolates the inference
effect. Cross-cohort GO replication (preserved aging modules with GO overlap vs a
permutation null) was additionally computed for isograph vs classical wgcna_gene. Analyses
used the project Python 3.12 environment with deterministic seeds.

## Results text

**Phenotype signal lives in the switch features, not the method.** Pooled per-module
phenotype-significant rate: `wgcna_switch_only` **0.336** > `isograph` **0.268** >
`wgcna_multiplex` **0.189** ≈ `wgcna_gene` **0.180**. Both switch-fed methods beat both
abundance-fed methods — representing isoform switching is what buys phenotype sensitivity,
independent of the network-inference method.

**IsoGraph is NOT globally superior on module-level metrics.** Classical
`wgcna_switch_only` matches or exceeds IsoGraph's phenotype-significant rate, and IsoGraph
has the **lowest** "both" (phenotype-sig AND GO) rate (0.074 vs 0.136–0.196) because its
GO-enrichment is low by construction. GO-enriched rate is abundance-dominated:
`wgcna_gene` 0.885 > `wgcna_multiplex` 0.758 ≫ `wgcna_switch_only` 0.381 > `isograph`
0.217. This is the expected picture — abundance dominates gene-level GO enrichment — and is
honest, not a failure.

**One clean method effect survives.** On the identical switch+abundance multiplex features,
`isograph` (0.268) exceeds `wgcna_multiplex` (0.189): VAE + Leiden extracts more
phenotype-linked switch structure from the full multiplex than classical WGCNA, which
dilutes the switch signal back toward the abundance baseline when abundance is added. The
representation+inference combination, not either alone, is where IsoGraph adds module-level
value.

**Cross-cohort GO replication favors abundance.** For preserved aging modules with GO
overlap, classical `wgcna_gene` carries more cross-cohort GO Jaccard (median 0.077,
perm p=1e-3) than `isograph` (median 0.0, mean 0.048, perm p=0.013) — again the abundance/GO
advantage — though both beat their permutation null. (Replication covers isograph vs
wgcna_gene only; the matched baselines were not run through replication_go.)

**Headline:** *On per-module rates IsoGraph is not globally superior to WGCNA — phenotype
sensitivity comes from the switch features (both switch-fed methods win), GO enrichment is
abundance-dominated, and the only clean method effect is that VAE + Leiden beats classical
multiplex WGCNA on identical features. IsoGraph's defensible value is therefore the
DTU-without-DGE content, not better module-level enrichment.*

## Figure and table notes

- **Supplementary figure — three-baseline rates
  (`manuscript/_m/figures/figBaselineRates.{pdf,png}`, built by
  `manuscript/_h/baseline_rates_figure.R`).** Faceted grouped bars: phenotype-sig rate /
  both rate / GO-enriched rate per method, fill encoded by feature class
  (switch+abundance / switch-only / abundance) so the switch-vs-abundance story reads off
  the colour. No in-panel titles; facet strips carry the metric, caption carries the read.
  Key message: switch features win phenotype rate; abundance wins GO; IsoGraph is not
  globally superior.
- **Supplementary table:** `baseline_comparison.parquet` (per cohort × region × method:
  n_modules, median size, frac_pheno_sig, frac_both, frac_go_enriched) and
  `baseline_comparison_pooled.parquet` (pooled per method).

## Reproducibility information

- Analysis directory: `04_module_characterization/_m/baseline_comparison/`.
- Primary script: `isograph_benchmark/real_data/baseline_comparison.py`
  (`python -m isograph_benchmark.real_data.baseline_comparison`; login-node aggregation,
  no SLURM — reads saved module_enrichment + replication_go outputs).
- Inputs: per-region `module_enrichment` for all four methods; `replication_go` (isograph,
  wgcna_gene).
- Outputs: `baseline_comparison.parquet`, `baseline_comparison_pooled.parquet`,
  `BASELINE_COMPARISON.md`.
- Key parameters: phenotype FDR ≤ 0.1; GO-enriched = `n_go_terms > 0`; per-module rates
  averaged across analyses; cross-cohort GO replication vs permutation null.
- Compute environment: project Python 3.12
  (`/ocean/projects/bio260021p/shared/opt/envs/isograph`).
- Missing reproducibility information: per-package versions from the live environment, not a
  per-run lockfile; matched baselines absent from the replication_go arm.

## Limitations and integration notes

- **Read rates, not totals.** Raw phenotype-sig totals (isograph 158 > wgcna_switch_only 75
  > wgcna_multiplex 53 > wgcna_gene 34) scale with module count and are NOT a superiority
  claim; IsoGraph's finer partition inflates totals. Only per-module rates are comparable.
- The matched WGCNA baselines were not run through cross-cohort GO replication, so that
  contrast is isograph vs classical abundance WGCNA only.
- This analysis deliberately **bounds** the claim: it shows IsoGraph does not win
  module-level enrichment/phenotype rates outright. Its purpose is to force the manuscript to
  rest the value claim on DTU-without-DGE content. Integrate with the GO-invisible gate
  (the content is real disease isoform switching), the QTL splicing-specificity contrast
  (the same content is genetically anchored, and there the matched-baseline comparison DOES
  show a clean IsoGraph-only effect), and the de-confounded gene-level result (abundance
  dominates the bulk signal). Together: IsoGraph is a complementary layer, not a superior
  one.
- This is a **scope-setting primary** analysis (it defines the honest boundary of the
  central claim), not a sensitivity check.
