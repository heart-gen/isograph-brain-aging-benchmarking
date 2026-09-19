# Module-level colocalization convergence

Do GWAS colocalizing switch genes concentrate in particular IsoGraph co-switch modules, beyond what module sizes alone predict? Generalises the SCZ-only layer in `scz_age_projection.py` to every trait, so the four aging traits are reported on the same partitions and with the same statistics.

Reproduce: `python -m isograph_benchmark.real_data.module_coloc_convergence`

## Reading

- **Per-module counts are small.** With single-digit colocalizing genes per trait, the per-module hypergeometric is underpowered; the global concentration permutation is the better-powered test and should lead.
- **`n_coloc_genes` vs `n_coloc_loci`.** Several genes can sit under one GWAS peak. A module whose genes collapse to one locus is a single-locus result, not a converging program. `loo_worst_p` drops each contributing locus in turn and reports the worst case — read it before believing any module row.
- Module assignment is re-derived from each source's own `modules.parquet`; the `module_id` column in the coloc outputs is an arbitrary contributing analysis's label for bundle gene sources and is deliberately not used.
- The permutation holds module sizes fixed by construction, so a large module collecting hits in proportion to its size is not evidence.

## Global concentration per (trait, source)

| trait | source | n_pool | n_coloc_genes | n_modules_with_coloc | concentration_obs | concentration_null_mean | concentration_p | n_anchored_modules | frac_coloc_anchored | frac_pool_anchored | anchored_hyperg_p | frac_allgenes_anchored | anchored_hyperg_p_allgenes_denom |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AD | gtex_caudate_bg | 500 | 13 | 5.0 | 35.0 | 33.694 | 0.412 | 4.0 | 0.154 | 0.218 | 0.815 | 0.20655983975963946 | 0.7836686172748338 |
| AD | brainseq_caudate | 391 | 7 | 4.0 | 13.0 | 12.136 | 0.43 | 2.0 | 0.571 | 0.322 | 0.155 | 0.29840116279069767 | 0.12385952279236992 |
| PD | gtex_caudate_bg | 302 | 10 | 5.0 | 28.0 | 22.805 | 0.202 | 3.0 | 0.4 | 0.242 | 0.202 | 0.2398597896845268 | 0.20081723548817945 |
| PD | brainseq_caudate | 224 | 7 | 4.0 | 15.0 | 12.062 | 0.235 | 1.0 | 0.0 | 0.134 | 1 | 0.13066860465116278 | 1.0 |
| LBD | gtex_caudate_bg | 147 | 3 | 2.0 | 5.0 | 3.821 | 0.358 | 3.0 | 0.0 | 0.014 | 1 | 0.03242363545317977 | 1.0 |
| LBD | brainseq_caudate | 109 | 0 |  |  |  | NA |  |  |  | NA |  |  |
| ALS | gtex_caudate_bg | 389 | 24 | 8.0 | 96.0 | 101.561 | 0.579 | 6.0 | 0.458 | 0.398 | 0.34 | 0.37205808713069605 | 0.25039837171750395 |
| ALS | brainseq_caudate | 297 | 20 | 9.0 | 64.0 | 69.747 | 0.659 | 1.0 | 0.3 | 0.212 | 0.231 | 0.1677325581395349 | 0.1040290601076745 |
| SCZ | gtex_caudate_bg | 807 | 43 | 9.0 | 337.0 | 313.8 | 0.305 | 5.0 | 0.698 | 0.713 | 0.66 | 0.6402103154732098 | 0.26826619042878963 |
| SCZ | brainseq_caudate | 641 | 37 | 10.0 | 173.0 | 200.501 | 0.858 | 4.0 | 0.459 | 0.488 | 0.702 | 0.5197674418604651 | 0.8162540456030424 |

## Modules carrying colocalizing switch genes

| trait | source | module_id | testable | n_module_genes_in_pool | n_coloc_genes | n_coloc_loci | n_go_invisible | hyperg_p | hyperg_fdr | loo_worst_p | magma_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AD | gtex_caudate_bg | M005 | True | 34 | 3 | 3 | 1 | 0.0514 | 0.257 | 0.193 | 0.977 |
| AD | gtex_caudate_bg | M004 | True | 46 | 3 | 3 | 0 | 0.108 | 0.271 | 0.304 | 0.252 |
| AD | gtex_caudate_bg | M000 | True | 121 | 3 | 3 | 2 | 0.646 | 0.724 | 0.83 | 0.48 |
| AD | gtex_caudate_bg | M003 | True | 56 | 2 | 2 | 2 | 0.438 | 0.724 | 0.764 | 0.528 |
| AD | gtex_caudate_bg | M001 | True | 92 | 2 | 2 | 1 | 0.724 | 0.724 | 0.915 | 0.0452 |
| AD | brainseq_caudate | M003 | True | 55 | 2 | 2 | 0 | 0.257 | 0.372 | 0.6 | 0.042 |
| AD | brainseq_caudate | M000 | True | 66 | 2 | 2 | 1 | 0.336 | 0.372 | 0.673 | 0.581 |
| AD | brainseq_caudate | M001 | True | 71 | 2 | 2 | 2 | 0.372 | 0.372 | 0.702 | 0.00695 |
| AD | brainseq_caudate | M008 | True | 20 | 1 | 1 | 1 | 0.31 | 0.372 | 1 | 0.713 |
| PD | gtex_caudate_bg | M000 | True | 69 | 4 | 3 | 1 | 0.173 | 0.806 | 0.58 | 0.0394 |
| PD | gtex_caudate_bg | M001 | True | 63 | 3 | 3 | 0 | 0.348 | 0.806 | 0.593 | 0.343 |
| PD | gtex_caudate_bg | M003 | True | 28 | 1 | 1 | 0 | 0.628 | 0.806 | 1 | 0.443 |
| PD | gtex_caudate_bg | M005 | True | 32 | 1 | 1 | 0 | 0.68 | 0.806 | 1 | 0.445 |
| PD | gtex_caudate_bg | M002 | True | 45 | 1 | 1 | 1 | 0.806 | 0.806 | 1 | 0.49 |
| PD | brainseq_caudate | M008 | True | 8 | 3 | 3 | 0 | 0.00099 | 0.00396 | 0.0156 | 0.0776 |
| PD | brainseq_caudate | M004 | True | 24 | 2 | 2 | 0 | 0.166 | 0.331 | 0.497 | 0.644 |
| PD | brainseq_caudate | M002 | True | 30 | 1 | 1 | 1 | 0.64 | 0.701 | 1 | 0.27 |
| PD | brainseq_caudate | M001 | True | 35 | 1 | 1 | 0 | 0.701 | 0.701 | 1 | 0.0995 |
| LBD | gtex_caudate_bg | M001 | True | 28 | 2 | 2 | 0 | 0.093 | 0.186 | 0.346 | 0.742 |
| LBD | gtex_caudate_bg | M002 | True | 19 | 1 | 1 | 0 | 0.342 | 0.342 | 1 | 0.0809 |
| ALS | gtex_caudate_bg | M000 | True | 98 | 6 | 6 | 1 | 0.591 | 0.795 | 0.731 | 0.00167 |
| ALS | gtex_caudate_bg | M002 | True | 52 | 5 | 5 | 4 | 0.204 | 0.795 | 0.37 | 0.0311 |
| ALS | gtex_caudate_bg | M005 | True | 31 | 4 | 3 | 2 | 0.113 | 0.795 | 0.429 | 0.771 |
| ALS | gtex_caudate_bg | M001 | True | 72 | 3 | 3 | 0 | 0.856 | 0.856 | 0.949 | 0.114 |
| ALS | gtex_caudate_bg | M004 | True | 25 | 2 | 2 | 1 | 0.467 | 0.795 | 0.793 | 0.271 |
| ALS | gtex_caudate_bg | M003 | True | 46 | 2 | 2 | 0 | 0.803 | 0.856 | 0.942 | 0.491 |
| ALS | gtex_caudate_bg | M008 | True | 6 | 1 | 1 | 0 | 0.319 | 0.795 | 1 | 0.364 |
| ALS | gtex_caudate_bg | M007 | True | 14 | 1 | 1 | 0 | 0.596 | 0.795 | 1 | 0.0986 |
| ALS | brainseq_caudate | M001 | True | 63 | 6 | 5 | 0 | 0.231 | 0.542 | 0.555 | 0.00795 |
| ALS | brainseq_caudate | M000 | True | 54 | 3 | 3 | 1 | 0.741 | 0.834 | 0.893 | 0.0625 |
| ALS | brainseq_caudate | M010 | True | 6 | 2 | 2 | 0 | 0.055 | 0.495 | 0.33 | 0.647 |
| ALS | brainseq_caudate | M008 | True | 11 | 2 | 2 | 1 | 0.164 | 0.542 | 0.523 | 0.291 |
| ALS | brainseq_caudate | M006 | True | 14 | 2 | 2 | 1 | 0.241 | 0.542 | 0.612 | 0.0857 |
| ALS | brainseq_caudate | M003 | True | 34 | 2 | 2 | 2 | 0.695 | 0.834 | 0.908 | 0.167 |
| ALS | brainseq_caudate | M009 | True | 9 | 1 | 1 | 0 | 0.471 | 0.834 | 1 | 0.332 |
| ALS | brainseq_caudate | M007 | True | 15 | 1 | 1 | 1 | 0.658 | 0.834 | 1 | 0.327 |
| ALS | brainseq_caudate | M002 | True | 43 | 1 | 1 | 0 | 0.961 | 0.961 | 1 | 0.571 |
| SCZ | gtex_caudate_bg | M000 | True | 216 | 12 | 12 | 3 | 0.491 | 0.689 | 0.595 | 4.73e-05 |
| SCZ | gtex_caudate_bg | M003 | True | 77 | 9 | 8 | 5 | 0.0157 | 0.141 | 0.0768 | 0.271 |
| SCZ | gtex_caudate_bg | M004 | True | 63 | 6 | 5 | 3 | 0.11 | 0.494 | 0.4 | 0.0476 |
| SCZ | gtex_caudate_bg | M002 | True | 109 | 6 | 6 | 4 | 0.536 | 0.689 | 0.693 | 0.0257 |
| SCZ | gtex_caudate_bg | M001 | True | 162 | 6 | 6 | 3 | 0.894 | 0.937 | 0.948 | 0.00435 |
| SCZ | gtex_caudate_bg | M012 | True | 5 | 1 | 1 | 0 | 0.24 | 0.689 | 1 | 0.0739 |
| SCZ | gtex_caudate_bg | M010 | True | 7 | 1 | 1 | 0 | 0.319 | 0.689 | 1 | 0.165 |
| SCZ | gtex_caudate_bg | M006 | True | 13 | 1 | 1 | 0 | 0.512 | 0.689 | 1 | 0.451 |
| SCZ | gtex_caudate_bg | M005 | True | 49 | 1 | 1 | 0 | 0.937 | 0.937 | 1 | 0.0643 |
| SCZ | brainseq_caudate | M001 | True | 123 | 7 | 7 | 0 | 0.587 | 0.838 | 0.721 | 0.00165 |
| SCZ | brainseq_caudate | M004 | True | 79 | 6 | 6 | 3 | 0.299 | 0.685 | 0.463 | 0.166 |
| SCZ | brainseq_caudate | M003 | True | 82 | 6 | 5 | 4 | 0.332 | 0.685 | 0.678 | 1.84e-05 |
| SCZ | brainseq_caudate | M000 | True | 98 | 4 | 4 | 1 | 0.846 | 0.933 | 0.934 | 0.00324 |
| SCZ | brainseq_caudate | M006 | True | 30 | 3 | 3 | 1 | 0.247 | 0.685 | 0.513 | 0.0521 |
| SCZ | brainseq_caudate | M005 | True | 40 | 3 | 3 | 3 | 0.411 | 0.685 | 0.674 | 0.211 |
| SCZ | brainseq_caudate | M002 | True | 95 | 3 | 3 | 1 | 0.933 | 0.933 | 0.98 | 0.216 |
| SCZ | brainseq_caudate | M009 | True | 18 | 2 | 2 | 1 | 0.279 | 0.685 | 0.652 | 0.494 |
| SCZ | brainseq_caudate | M007 | True | 24 | 2 | 2 | 0 | 0.409 | 0.685 | 0.757 | 0.776 |
| SCZ | brainseq_caudate | M008 | True | 25 | 1 | 1 | 1 | 0.78 | 0.933 | 1 | 0.109 |

