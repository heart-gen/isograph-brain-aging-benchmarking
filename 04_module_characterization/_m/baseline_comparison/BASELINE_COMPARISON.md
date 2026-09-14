# Three-baseline comparison — IsoGraph vs matched WGCNA baselines

Per-region module_enrichment across 17 analyses (17 with all four methods). Phenotype-significant = `pheno_fdr <= 0.1`; GO-enriched = `n_go_terms > 0`; BOTH = phenotype-significant AND GO-enriched. Rates are per-module fractions averaged across regions; raw totals scale with module count (IsoGraph runs at finer resolution) and are NOT directly comparable.

Reproduce: `python -m isograph_benchmark.real_data.baseline_comparison`.

## Pooled across regions (per method)

| method | features | n_regions | med_n_mod | med_size | pheno_sig_rate | both_rate | go_rate | tot_pheno_sig |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | switch+abundance | 17 | 39.0 | 70.0 | 0.274 | 0.082 | 0.235 | 180 |
| wgcna_gene | abundance | 17 | 8.0 | 387.0 | 0.18 | 0.171 | 0.885 | 34 |
| wgcna_switch_only | switch-only | 17 | 13.0 | 144.0 | 0.389 | 0.197 | 0.401 | 77 |
| wgcna_multiplex | switch+abundance | 17 | 20.0 | 293.0 | 0.187 | 0.137 | 0.755 | 61 |

## Caudate (representative single region)

| method | features | n_modules | frac_pheno_sig | frac_both | frac_go_enriched |
| --- | --- | --- | --- | --- | --- |
| isograph | switch+abundance | 46 | 0.304 | 0.0 | 0.109 |
| wgcna_gene | abundance | 21 | 0.048 | 0.048 | 0.81 |
| wgcna_switch_only | switch-only | 10 | 0.3 | 0.1 | 0.3 |
| wgcna_multiplex | switch+abundance | 27 | 0.111 | 0.111 | 0.815 |

## The features-vs-method question

`wgcna_switch_only` and `wgcna_multiplex` are fed the SAME switch / switch+abundance feature matrix as IsoGraph but inferred with classical WGCNA, so the comparison isolates network inference (VAE + Leiden) from the input representation. `wgcna_gene` is the abundance-only input control.

## Honest read (per-module rates, not totals)

- **The phenotype signal lives in the switch features, not the method.** Both switch-fed methods (isograph, wgcna_switch_only) carry a higher phenotype-significant rate than the abundance-fed ones (wgcna_gene, wgcna_multiplex). Representing isoform switching is what buys phenotype sensitivity.
- **IsoGraph is NOT globally superior on module-level metrics.** Classical `wgcna_switch_only` matches or exceeds IsoGraph's phenotype-significant rate, and IsoGraph has the LOWEST BOTH rate (few modules are both phenotype-sig and GO-enriched, because its GO-enrichment is low by construction). This is the expected picture: abundance dominates module-level enrichment.
- **One clean method effect survives:** on IDENTICAL switch+abundance features, IsoGraph's phenotype-significant rate exceeds `wgcna_multiplex` — VAE + Leiden extracts more phenotype-linked structure from the full multiplex than classical WGCNA, which dilutes the switch signal back toward the abundance baseline when abundance is added.
- **Classical / multiplex WGCNA win GO-enrichment** (gene-level, abundance-biased GO). IsoGraph's low GO fraction is expected, not a failure — see the GO-invisible biology gate.
- **Bottom line:** IsoGraph's defensible value is the DTU-without-DGE *content* (genes/modules structurally invisible to any abundance pipeline — biology gate, incremental association, sQTL-vs-eQTL specificity), not better module-level enrichment or phenotype rates than every WGCNA baseline. Frame the method as a complementary layer, consistent with the de-confounded gene-level result.

## Cross-cohort GO replication consistency

| method | n_preserved_aging_pairs | n_pairs_with_go | mean_go_jaccard | median_go_jaccard | perm_p | null_mean_go_jaccard |
| --- | --- | --- | --- | --- | --- | --- |
| isograph | 32 | 19 | 0.0529 | 0.0 | 2.00e-03 | 0.003 |
| wgcna_gene | 26 | 26 | 0.2005 | 0.0767 | 9.99e-04 | 0.0395 |

Replication covers isograph vs classical `wgcna_gene` only (the matched baselines were not run through replication_go). Classical WGCNA's preserved aging modules carry more cross-cohort GO overlap — again the abundance/GO advantage — while both beat their permutation null.
