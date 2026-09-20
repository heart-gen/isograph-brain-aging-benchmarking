# Three-baseline module comparison (features vs method; honest scope)

Modular analysis summary for Manubot integration. **Hand-written** — it is NOT emitted by
`baseline_comparison.py`, so it does not update when the analysis re-runs. Re-quote it by hand
against `04_module_trust/_m/baseline_comparison/baseline_comparison{,_pooled}.parquet` and the
generated `BASELINE_COMPARISON.md` whenever those change. *(Last re-quoted 2026-09-19 against
the switching-filter re-run at Leiden 2.0; it had carried the legacy expression-filter rates.)*

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
totals scale with module count and the four methods partition at different granularity
(IsoGraph median 25 modules of median size 108; WGCNA baselines median 8–25 modules, size
159–468), making raw counts non-comparable. A module was phenotype-significant at `pheno_fdr ≤ 0.1`, GO-enriched
at `n_go_terms > 0`, and "both" if it was phenotype-significant and GO-enriched. Per-module
fractions were averaged across analyses to give pooled rates per method. The phenotype rate
contrasts switch-fed (isograph, wgcna_switch_only) against abundance-fed (wgcna_gene,
wgcna_multiplex) methods to test where phenotype sensitivity comes from; the
isograph-vs-wgcna_multiplex contrast on identical multiplex features isolates the inference
effect. Cross-cohort GO replication (preserved aging modules with GO overlap vs a
permutation null) was additionally computed for isograph vs classical wgcna_gene. Analyses
used the project Python 3.12 environment with deterministic seeds.

## Results text

**The switch-only representation carries the highest phenotype rate — in GTEx.** Pooled
per-module phenotype-significant rate: `wgcna_switch_only` **0.214** > `wgcna_multiplex`
**0.189** > `wgcna_gene` **0.178** > `isograph` **0.146**. This pooled figure is the unweighted
mean of per-analysis fractions over 17 analyses, 13 of which are GTEx, and **the two cohorts
rank the methods in opposite orders**:

| method | BrainSEQ (k=4) | GTEx (k=13) | pooled |
|---|---|---|---|
| isograph | **0.176** | 0.137 | 0.146 |
| wgcna_multiplex | 0.163 | 0.196 | 0.189 |
| wgcna_gene | 0.094 | 0.203 | 0.178 |
| wgcna_switch_only | 0.046 | **0.265** | 0.214 |

Quote the split, not the pooled ranking. Two further cautions: 4–8 of the 17 analyses return
zero phenotype-significant modules for a given method (many GTEx regions are null for every
method and drag all four means down equally), and module counts per analysis range from 3 to
35, so one module moves a small-denominator fraction by tens of points.

**IsoGraph is NOT superior on any module-level rate.** It is last on pooled phenotype rate
(though first in BrainSEQ — see above) and has the **lowest** "both" (phenotype-sig AND GO) rate
(0.052 vs 0.117–0.159), because its GO enrichment is low by construction. GO-enriched rate is abundance-dominated: `wgcna_gene`
**0.862** > `wgcna_multiplex` 0.649 ≫ `isograph` 0.253 > `wgcna_switch_only` 0.184. This is
the expected picture — abundance dominates gene-level GO enrichment — and is honest, not a
failure.

**The legacy "one clean method effect" did not survive the re-run.** On the identical
switch+abundance multiplex features the legacy expression-filter production had `isograph`
(0.268) above `wgcna_multiplex` (0.189); on the switching filter at Leiden 2.0 the contrast
inverts (0.146 vs 0.189). There is now no per-module rate on which IsoGraph leads, and the
manuscript must not claim an inference-level module-rate advantage.

**Cross-cohort GO replication favors abundance.** For preserved aging modules with GO
overlap, classical `wgcna_gene` carries more cross-cohort GO Jaccard (median 0.077,
perm p=1e-3) than `isograph` (median 0.0, mean 0.048, perm p=0.013) — again the abundance/GO
advantage — though both beat their permutation null. (Replication covers isograph vs
wgcna_gene only; the matched baselines were not run through replication_go.)

**Headline:** *On per-module rates IsoGraph is not superior to WGCNA on any metric —
GO enrichment is abundance-dominated, and no inference-level advantage survives on identical
features in the pooled mean. Phenotype rate is cohort-dependent and should not be quoted as a
single ranking: IsoGraph leads in BrainSEQ and the switch-only baseline leads in GTEx.
IsoGraph's defensible value is therefore the DTU-without-DGE content and the genetic
anchoring, not module-level enrichment rates.*

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

- Analysis directory: `04_module_trust/_m/baseline_comparison/`.
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

- **Read rates, not totals.** Raw phenotype-sig totals (wgcna_multiplex 75 >
  wgcna_switch_only 64 > isograph 54 > wgcna_gene 28) scale with module count and are NOT a
  superiority claim for anyone. Only per-module rates are comparable.
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
