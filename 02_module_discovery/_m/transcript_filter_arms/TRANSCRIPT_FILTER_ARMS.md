# Transcript-filter arms (BrainSEQ)

`transcript_filter_arms.py`: the production IsoGraph fit re-run with only the transcript filter changed, written outside the production stores. `production` is the committed fit. Aging production uses the abundance filter (count > 10 in ≥ 70%), so its `abundance` refit is the noise floor; SCZD production is unfiltered, so its `none` refit is the floor and `abundance` is the BrainSEQ method applied to SCZD.

Switching filter: gene count ≥ 10 in ≥ 70% of samples; transcript count ≥ 10 and share of its gene ≥ 0.10, each in ≥ 10% of samples.

Trait tables: aging `age_linear` (linear age), SCZD `diagnosis_assoc`; significant = FDR < 0.05.

## Fits

| region | arm | seed | n_transcripts | n_genes_input | n_multi_isoform_genes | n_modules | n_genes_assigned | median_module_size | max_module_size | largest_module_frac | n_modules_ge_900 | n_trait_tested | n_trait_sig | n_genes_in_trait_sig | selected_alpha_abundance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| caudate | production | 13 |  |  |  | 46 | 5795 | 54 | 687 | 0.119 | 0 | 46 | 20 | 3115 | 0.95 |
| caudate | abundance | 13 | 5.91e+04 | 1.8e+04 | 1.16e+04 | 46 | 5795 | 54 | 687 | 0.119 | 0 | 46 | 20 | 3115 | 0.95 |
| caudate | switching | 13 | 5.04e+04 | 1.93e+04 | 1.32e+04 | 28 | 5742 | 53.5 | 797 | 0.139 | 0 | 28 | 4 | 1461 | 0.95 |
| hippocampus | production | 13 |  |  |  | 50 | 3941 | 46 | 446 | 0.113 | 0 | 50 | 47 | 3734 | 0.95 |
| hippocampus | abundance | 13 | 5.56e+04 | 1.75e+04 | 1.11e+04 | 50 | 3941 | 46 | 446 | 0.113 | 0 | 50 | 47 | 3734 | 0.95 |
| hippocampus | switching | 13 | 5.1e+04 | 1.89e+04 | 1.31e+04 | 32 | 4112 | 81 | 544 | 0.132 | 0 | 32 | 29 | 3888 | 0.9 |
| dlpfc | production | 13 |  |  |  | 39 | 4458 | 68 | 410 | 0.092 | 0 | 39 | 6 | 659 | 0.95 |
| dlpfc | abundance | 13 | 6.04e+04 | 1.81e+04 | 1.16e+04 | 39 | 4458 | 68 | 410 | 0.092 | 0 | 39 | 6 | 659 | 0.95 |
| dlpfc | switching | 13 | 4.93e+04 | 1.95e+04 | 1.29e+04 | 38 | 5349 | 84.5 | 658 | 0.123 | 0 | 38 | 4 | 650 | 0.95 |
| caudate_sczd | production | 13 |  |  |  | 30 | 3838 | 40.5 | 719 | 0.187 | 0 | 30 | 7 | 315 | 0.95 |
| caudate_sczd | none | 13 | 2.33e+05 | 2.06e+04 | 1.72e+04 | 30 | 3838 | 40.5 | 719 | 0.187 | 0 | 30 | 7 | 315 | 0.95 |
| caudate_sczd | switching | 13 | 5.08e+04 | 1.94e+04 | 1.32e+04 | 33 | 6169 | 63 | 795 | 0.129 | 0 | 33 | 7 | 461 | 0.95 |
| caudate_sczd | abundance | 13 | 5.86e+04 | 1.79e+04 | 1.16e+04 | 42 | 5649 | 77.5 | 697 | 0.123 | 0 | 42 | 0 | 0 | 0.95 |

## Pairwise

Read every production-vs-arm row against the production-vs-floor row for the same region. A floor at ARI = 1 means the refit is deterministic at the production seed (identical edges), NOT that seed or input-perturbation variation is zero; the floor bounds refit reproducibility only. `n_ref_sig_retained_in_qry`: significant `ref` modules whose best-Jaccard `qry` counterpart is significant with the same sign (and the reverse column); it is base-rate blind, so read it with `median_best_jaccard_*` and the gene-level columns. `frac_ref_sig_genes_kept_same_sign` against `qry_sig_gene_rate` (its chance level) and `sig_gene_odds_ratio` are module-free; `gene_effect_spearman` correlates each gene's module trait effect across the two fits. `trait_gene_jaccard`: overlap of the genes in significant modules. ARI/NMI and the gene-level columns are over genes assigned in both fits.

| region | ref | qry | ref_seed | qry_seed | n_genes_both | assigned_gene_jaccard | ari | nmi | n_trait_sig_ref | n_trait_sig_qry | n_ref_sig_retained_in_qry | median_best_jaccard_ref_sig | n_qry_sig_retained_in_ref | median_best_jaccard_qry_sig | trait_gene_jaccard | gene_effect_spearman | n_ref_sig_genes | frac_ref_sig_genes_kept_same_sign | qry_sig_gene_rate | sig_gene_odds_ratio | transcript_jaccard |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| caudate | production | abundance | 13 | 13 | 5795 | 1 | 1 | 1 | 20 | 20 | 20 | 1 | 20 | 1 | 1 | 1 | 3115 | 1 | 0.538 | inf |  |
| caudate | production | switching | 13 | 13 | 3560 | 0.446 | 0.337 | 0.38 | 20 | 4 | 15 | 0.0346 | 4 | 0.058 | 0.252 | 0.624 | 2006 | 0.456 | 0.285 | 13 |  |
| caudate | abundance | switching | 13 | 13 | 3560 | 0.446 | 0.337 | 0.38 | 20 | 4 | 15 | 0.0346 | 4 | 0.058 | 0.252 | 0.624 | 2006 | 0.456 | 0.285 | 13 | 0.507 |
| hippocampus | production | abundance | 13 | 13 | 3941 | 1 | 1 | 1 | 47 | 47 | 47 | 1 | 47 | 1 | 1 | 1 | 3734 | 1 | 0.947 | inf |  |
| hippocampus | production | switching | 13 | 13 | 2339 | 0.409 | 0.306 | 0.378 | 47 | 29 | 46 | 0.0588 | 28 | 0.0872 | 0.388 | 0.434 | 2223 | 0.958 | 0.951 | 5.36 |  |
| hippocampus | abundance | switching | 13 | 13 | 2339 | 0.409 | 0.306 | 0.378 | 47 | 29 | 46 | 0.0588 | 28 | 0.0872 | 0.388 | 0.434 | 2223 | 0.958 | 0.951 | 5.36 | 0.509 |
| dlpfc | production | abundance | 13 | 13 | 4458 | 1 | 1 | 1 | 6 | 6 | 6 | 1 | 6 | 1 | 1 | 1 | 659 | 1 | 0.148 | inf |  |
| dlpfc | production | switching | 13 | 13 | 3045 | 0.45 | 0.263 | 0.362 | 6 | 4 | 3 | 0.131 | 3 | 0.133 | 0.089 | 0.319 | 470 | 0.213 | 0.123 | 2.53 |  |
| dlpfc | abundance | switching | 13 | 13 | 3045 | 0.45 | 0.263 | 0.362 | 6 | 4 | 3 | 0.131 | 3 | 0.133 | 0.089 | 0.319 | 470 | 0.213 | 0.123 | 2.53 | 0.499 |
| caudate_sczd | production | none | 13 | 13 | 3838 | 1 | 1 | 1 | 7 | 7 | 7 | 1 | 7 | 1 | 1 | 1 | 315 | 1 | 0.0821 | inf |  |
| caudate_sczd | production | switching | 13 | 13 | 2434 | 0.321 | 0.389 | 0.373 | 7 | 7 | 2 | 0.024 | 1 | 0.0323 | 0.0197 | 0.664 | 99 | 0.131 | 0.0468 | 4.03 |  |
| caudate_sczd | production | abundance | 13 | 13 | 2122 | 0.288 | 0.382 | 0.412 | 7 | 0 | 0 | 0.025 | 0 |  | 0 | 0.665 | 93 | 0 | 0 |  |  |
| caudate_sczd | none | switching | 13 | 13 | 2434 | 0.321 | 0.389 | 0.373 | 7 | 7 | 2 | 0.024 | 1 | 0.0323 | 0.0197 | 0.664 | 99 | 0.131 | 0.0468 | 4.03 | 0.218 |
| caudate_sczd | none | abundance | 13 | 13 | 2122 | 0.288 | 0.382 | 0.412 | 7 | 0 | 0 | 0.025 | 0 |  | 0 | 0.665 | 93 | 0 | 0 |  | 0.251 |
| caudate_sczd | switching | abundance | 13 | 13 | 3757 | 0.466 | 0.37 | 0.4 | 7 | 0 | 0 | 0.0638 | 0 |  | 0 | 0.663 | 256 | 0 | 0 |  | 0.508 |

## Seed floor

Every arm refit at seeds 13, 14, 15 (`random_state` drives VAE initialisation, minibatch order and Leiden). `within_arm` pairs are seed-to-seed variation for one filter: the floor. `between_arm` pairs change the filter (and, off the diagonal, the seed). A filter effect is a `between_arm` value outside the `within_arm` range of both arms.

Per arm, values in seed order:

| region | arm | seeds | n_modules | n_trait_sig | n_genes_in_trait_sig | n_modules_ge_900 |
|---|---|---|---|---|---|---|
| caudate | abundance | 13 / 14 / 15 | 46 / 44 / 47 | 20 / 23 / 30 | 3115 / 3273 / 3559 | 0 |
| caudate | switching | 13 / 14 / 15 | 28 / 33 / 30 | 4 / 5 / 7 | 1461 / 1408 / 2185 | 0 |
| hippocampus | abundance | 13 / 14 / 15 | 50 / 44 / 47 | 47 / 40 / 42 | 3734 / 3080 / 3506 | 0 |
| hippocampus | switching | 13 / 14 / 15 | 32 / 33 / 32 | 29 / 30 / 30 | 3888 / 3501 / 3554 | 0 |
| dlpfc | abundance | 13 / 14 / 15 | 39 / 40 / 37 | 6 / 1 / 3 | 659 / 290 / 263 | 0 |
| dlpfc | switching | 13 / 14 / 15 | 38 / 34 / 38 | 4 / 6 / 7 | 650 / 835 / 617 | 0 |
| caudate_sczd | switching | 13 / 14 / 15 | 33 / 33 / 27 | 7 / 6 / 7 | 461 / 319 / 438 | 1 |
| caudate_sczd | none | 13 / 14 / 15 | 30 / 30 / 27 | 7 / 2 / 2 | 315 / 106 / 117 | 0 |
| caudate_sczd | abundance | 13 / 14 / 15 | 42 / 43 / 43 | 0 / 5 / 4 | 0 / 505 / 377 | 1 |

| region | kind | ref | qry | n_pairs | ari_median | ari_min | ari_max | gene_effect_spearman_median | frac_ref_sig_genes_kept_same_sign_median | qry_sig_gene_rate_median | trait_gene_jaccard_median |
|---|---|---|---|---|---|---|---|---|---|---|---|
| caudate | within_arm | abundance | abundance | 3 | 0.691 | 0.67 | 0.696 | 0.882 | 0.966 | 0.636 | 0.604 |
| caudate | between_arm | abundance | switching | 9 | 0.321 | 0.297 | 0.337 | 0.616 | 0.445 | 0.281 | 0.241 |
| caudate | within_arm | switching | switching | 3 | 0.722 | 0.719 | 0.746 | 0.928 | 0.939 | 0.378 | 0.535 |
| hippocampus | within_arm | abundance | abundance | 3 | 0.611 | 0.592 | 0.632 | 0.707 | 0.964 | 0.932 | 0.595 |
| hippocampus | between_arm | abundance | switching | 9 | 0.306 | 0.298 | 0.323 | 0.431 | 0.963 | 0.95 | 0.363 |
| hippocampus | within_arm | switching | switching | 3 | 0.638 | 0.63 | 0.701 | 0.74 | 0.988 | 0.948 | 0.67 |
| dlpfc | within_arm | abundance | abundance | 3 | 0.372 | 0.355 | 0.599 | 0.706 | 0.245 | 0.0706 | 0.155 |
| dlpfc | between_arm | abundance | switching | 9 | 0.224 | 0.173 | 0.263 | 0.319 | 0.302 | 0.134 | 0.0745 |
| dlpfc | within_arm | switching | switching | 3 | 0.497 | 0.484 | 0.69 | 0.555 | 0.522 | 0.142 | 0.204 |
| caudate_sczd | within_arm | none | none | 3 | 0.63 | 0.615 | 0.685 | 0.903 | 0.41 | 0.0273 | 0.268 |
| caudate_sczd | between_arm | none | switching | 9 | 0.308 | 0.282 | 0.396 | 0.652 | 0 | 0.0457 | 0.00184 |
| caudate_sczd | between_arm | none | abundance | 9 | 0.319 | 0.31 | 0.382 | 0.648 | 0 | 0.0441 | 0 |
| caudate_sczd | within_arm | switching | switching | 3 | 0.709 | 0.667 | 0.738 | 0.914 | 0.667 | 0.0617 | 0.264 |
| caudate_sczd | between_arm | switching | abundance | 9 | 0.347 | 0.285 | 0.37 | 0.647 | 0.18 | 0.042 | 0.0551 |
| caudate_sczd | within_arm | abundance | abundance | 3 | 0.603 | 0.595 | 0.739 | 0.842 | 0.382 | 0.0587 | 0 |
