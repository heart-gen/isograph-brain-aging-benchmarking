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
| AD | gtex_caudate_bg | 212 | 7 | 3.0 | 19.0 | 21.494 | 0.686 | 1.0 | 0.429 | 0.231 | 0.203 | 0.19966348850252383 | 0.1473246969623006 |
| AD | brainseq_caudate | 148 | 4 | 4.0 | 4.0 | 5.19 | 1 | 2.0 | 0.0 | 0.243 | 1 | 0.13375130616509928 | 1.0 |
| PD | gtex_caudate_bg | 118 | 11 | 4.0 | 37.0 | 49.121 | 0.861 | 1.0 | 0.364 | 0.508 | 0.908 | 0.3380071041316134 | 0.5408688038534893 |
| PD | brainseq_caudate | 72 | 4 | 3.0 | 6.0 | 5.451 | 0.53 | 3.0 | 0.25 | 0.278 | 0.737 | 0.14350400557297108 | 0.46194672782800716 |
| LBD | gtex_caudate_bg | 57 | 2 | 2.0 | 2.0 | 2.637 | 1 | 3.0 | 0.0 | 0.018 | 1 | 0.06169377453729669 | 1.0 |
| LBD | brainseq_caudate | 41 | 0 |  |  |  | NA |  |  |  | NA |  |  |
| ALS | gtex_caudate_bg | 144 | 9 | 3.0 | 41.0 | 36.258 | 0.371 | 2.0 | 0.667 | 0.556 | 0.37 | 0.4099831744251262 | 0.11070810967671942 |
| ALS | brainseq_caudate | 98 | 6 | 5.0 | 8.0 | 9.796 | 0.849 | 4.0 | 0.5 | 0.398 | 0.451 | 0.19592476489028213 | 0.09382987911555346 |
| SCZ | gtex_caudate_bg | 496 | 26 | 4.0 | 414.0 | 268.981 | 0.0103 | 3.0 | 0.962 | 0.879 | 0.153 | 0.5978687605159843 | 2.779959073695038e-05 |
| SCZ | brainseq_caudate | 365 | 21 | 9.0 | 93.0 | 70.658 | 0.113 | 8.0 | 0.429 | 0.463 | 0.708 | 0.3793103448275862 | 0.39853655090447354 |

## Modules carrying colocalizing switch genes

| trait | source | module_id | testable | n_module_genes_in_pool | n_coloc_genes | n_coloc_loci | n_go_invisible | hyperg_p | hyperg_fdr | loo_worst_p | magma_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AD | gtex_caudate_bg | M001 | True | 49 | 3 | 3 | 1 | 0.203 | 0.304 | 0.422 | 0.0367 |
| AD | gtex_caudate_bg | M000 | True | 112 | 3 | 3 | 2 | 0.821 | 0.821 | 0.918 | 0.424 |
| AD | gtex_caudate_bg | M009 | True | 3 | 1 | 1 | 1 | 0.0963 | 0.289 | 1 | 0.43 |
| AD | brainseq_caudate | M011 | True | 4 | 1 | 1 | 1 | 0.105 | 0.314 | 1 | 0.184 |
| AD | brainseq_caudate | M006 | True | 13 | 1 | 1 | 1 | 0.31 | 0.461 | 1 | 0.673 |
| AD | brainseq_caudate | M001 | True | 21 | 1 | 1 | 1 | 0.461 | 0.461 | 1 | 0.0671 |
| AD | brainseq_caudate | M021 | False | 1 | 1 | 1 | 0 | 0.027 | NA | 1 | 0.809 |
| PD | gtex_caudate_bg | M001 | True | 35 | 4 | 4 | 3 | 0.42 | 0.56 | 0.616 | 0.377 |
| PD | gtex_caudate_bg | M000 | True | 60 | 4 | 3 | 0 | 0.908 | 0.908 | 0.986 | 0.0383 |
| PD | gtex_caudate_bg | M004 | True | 7 | 2 | 2 | 0 | 0.128 | 0.513 | 0.471 | 0.415 |
| PD | gtex_caudate_bg | M003 | True | 4 | 1 | 1 | 1 | 0.328 | 0.56 | 1 | 0.429 |
| PD | brainseq_caudate | M003 | True | 7 | 2 | 2 | 1 | 0.0447 | 0.134 | 0.268 | 0.717 |
| PD | brainseq_caudate | M001 | True | 10 | 1 | 1 | 0 | 0.458 | 0.687 | 1 | 0.17 |
| PD | brainseq_caudate | M002 | True | 20 | 1 | 1 | 0 | 0.737 | 0.737 | 1 | 0.012 |
| LBD | gtex_caudate_bg | M003 | True | 5 | 1 | 1 | 1 | 0.169 | 0.338 | 1 | 0.4 |
| LBD | gtex_caudate_bg | M001 | True | 14 | 1 | 1 | 1 | 0.434 | 0.434 | 1 | 0.512 |
| ALS | gtex_caudate_bg | M000 | True | 77 | 6 | 6 | 0 | 0.321 | 0.481 | 0.439 | 0.000547 |
| ALS | gtex_caudate_bg | M001 | True | 44 | 2 | 2 | 2 | 0.823 | 0.823 | 0.951 | 0.171 |
| ALS | gtex_caudate_bg | M002 | True | 3 | 1 | 1 | 1 | 0.177 | 0.481 | 1 | 0.846 |
| ALS | brainseq_caudate | M002 | True | 28 | 2 | 2 | 0 | 0.553 | 0.553 | 0.822 | 0.0109 |
| ALS | brainseq_caudate | M007 | True | 8 | 1 | 1 | 1 | 0.409 | 0.553 | 1 | 0.0581 |
| ALS | brainseq_caudate | M013 | False | 1 | 1 | 1 | 1 | 0.0612 | NA | 1 | 0.133 |
| ALS | brainseq_caudate | M017 | False | 1 | 1 | 1 | 0 | 0.0612 | NA | 1 | 0.886 |
| ALS | brainseq_caudate | M019 | False | 2 | 1 | 1 | 0 | 0.119 | NA | 1 | 0.0299 |
| SCZ | gtex_caudate_bg | M000 | True | 269 | 20 | 19 | 7 | 0.0129 | 0.0514 | 0.0278 | 2.45e-05 |
| SCZ | gtex_caudate_bg | M001 | True | 135 | 3 | 3 | 2 | 0.987 | 0.987 | 0.997 | 0.00234 |
| SCZ | gtex_caudate_bg | M004 | True | 32 | 2 | 2 | 2 | 0.512 | 0.863 | 0.819 | 0.0294 |
| SCZ | gtex_caudate_bg | M005 | True | 19 | 1 | 1 | 1 | 0.647 | 0.863 | 1 | 0.322 |
| SCZ | brainseq_caudate | M002 | True | 91 | 8 | 7 | 4 | 0.122 | 0.506 | 0.327 | 0.00383 |
| SCZ | brainseq_caudate | M003 | True | 43 | 4 | 4 | 2 | 0.225 | 0.506 | 0.426 | 0.24 |
| SCZ | brainseq_caudate | M000 | True | 8 | 2 | 2 | 1 | 0.0716 | 0.506 | 0.366 | 0.139 |
| SCZ | brainseq_caudate | M001 | True | 56 | 2 | 2 | 2 | 0.862 | 0.862 | 0.968 | 0.102 |
| SCZ | brainseq_caudate | M025 | True | 4 | 1 | 1 | 0 | 0.212 | 0.506 | 1 | 0.167 |
| SCZ | brainseq_caudate | M009 | True | 8 | 1 | 1 | 1 | 0.38 | 0.685 | 1 | 0.591 |
| SCZ | brainseq_caudate | M008 | True | 15 | 1 | 1 | 0 | 0.596 | 0.862 | 1 | 0.0578 |
| SCZ | brainseq_caudate | M007 | True | 26 | 1 | 1 | 0 | 0.798 | 0.862 | 1 | 0.382 |
| SCZ | brainseq_caudate | M004 | True | 30 | 1 | 1 | 1 | 0.844 | 0.862 | 1 | 0.00493 |

