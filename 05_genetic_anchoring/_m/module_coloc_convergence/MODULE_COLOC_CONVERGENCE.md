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
| AD | gtex_caudate_bg | 177 | 8 | 5 | 16.0 | 13.029 | 0.241 | 5 | 0.625 | 0.401 | 0.17 | 0.27854806044555225 | 0.042755672357479714 |
| AD | brainseq_caudate | 106 | 4 | 3 | 6.0 | 4.883 | 0.358 | 4 | 0.5 | 0.245 | 0.252 | 0.20345125107851597 | 0.18607164957593011 |
| PD | gtex_caudate_bg | 131 | 9 | 6 | 17.0 | 14.814 | 0.309 | 3 | 0.222 | 0.168 | 0.466 | 0.0945630160461131 | 0.206493764708039 |
| PD | brainseq_caudate | 71 | 5 | 4 | 7.0 | 6.776 | 0.602 | 3 | 0.0 | 0.014 | 1 | 0.04227782571182054 | 1.0 |
| LBD | gtex_caudate_bg | 45 | 2 | 2 | 2.0 | 2.198 | 1 | 1 | 0.0 | 0.0 | NA | 0.011060912914784234 | 1.0 |
| LBD | brainseq_caudate | 30 | 1 | 1 | 1.0 | 1.0 | 1 | 2 | 0.0 | 0.0 | NA | 0.02985332182916307 | 1.0 |
| ALS | gtex_caudate_bg | 128 | 11 | 9 | 15.0 | 21.905 | 0.962 | 4 | 0.091 | 0.055 | 0.475 | 0.09300514098769279 | 0.6585940662942692 |
| ALS | brainseq_caudate | 71 | 3 | 3 | 3.0 | 3.451 | 1 | 4 | 0.333 | 0.155 | 0.401 | 0.12148403796376187 | 0.3220184577036663 |
| SCZ | gtex_caudate_bg | 482 | 27 | 14 | 87.0 | 96.118 | 0.647 | 5 | 0.444 | 0.38 | 0.302 | 0.22589188347094563 | 0.009526550908027911 |
| SCZ | brainseq_caudate | 303 | 18 | 11 | 52.0 | 42.323 | 0.188 | 5 | 0.333 | 0.343 | 0.627 | 0.2286453839516825 | 0.21227344361798486 |

## Modules carrying colocalizing switch genes

| trait | source | module_id | testable | n_module_genes_in_pool | n_coloc_genes | n_coloc_loci | n_go_invisible | hyperg_p | hyperg_fdr | loo_worst_p | magma_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AD | gtex_caudate_bg | M000 | True | 34 | 3 | 3 | 3 | 0.182 | 0.261 | 0.403 | 0.0292 |
| AD | gtex_caudate_bg | M010 | True | 13 | 2 | 2 | 1 | 0.109 | 0.261 | 0.419 | 0.0201 |
| AD | gtex_caudate_bg | M005 | True | 3 | 1 | 1 | 1 | 0.13 | 0.261 | 1 | 0.814 |
| AD | gtex_caudate_bg | M007 | True | 5 | 1 | 1 | 1 | 0.209 | 0.261 | 1 | 0.305 |
| AD | gtex_caudate_bg | M006 | True | 7 | 1 | 1 | 0 | 0.281 | 0.281 | 1 | 0.719 |
| AD | brainseq_caudate | M001 | True | 18 | 2 | 2 | 2 | 0.133 | 0.266 | 0.431 | 0.0411 |
| AD | brainseq_caudate | M003 | True | 18 | 1 | 1 | 1 | 0.531 | 0.531 | 1 | 0.188 |
| AD | brainseq_caudate | M021 | False | 2 | 1 | 1 | 1 | 0.0744 | NA | 1 | 0.078 |
| PD | gtex_caudate_bg | M000 | True | 27 | 3 | 3 | 3 | 0.274 | 0.411 | 0.52 | 0.664 |
| PD | gtex_caudate_bg | M004 | True | 12 | 2 | 2 | 2 | 0.193 | 0.411 | 0.547 | 0.306 |
| PD | gtex_caudate_bg | M005 | True | 4 | 1 | 1 | 1 | 0.25 | 0.411 | 1 | 0.201 |
| PD | gtex_caudate_bg | M029 | True | 4 | 1 | 1 | 1 | 0.25 | 0.411 | 1 | 1.41e-05 |
| PD | gtex_caudate_bg | M008 | True | 7 | 1 | 1 | 0 | 0.4 | 0.48 | 1 | 0.172 |
| PD | gtex_caudate_bg | M002 | True | 15 | 1 | 1 | 1 | 0.678 | 0.678 | 1 | 0.00191 |
| PD | brainseq_caudate | M001 | True | 14 | 2 | 2 | 1 | 0.254 | 0.366 | 0.593 | 0.0966 |
| PD | brainseq_caudate | M012 | True | 3 | 1 | 1 | 1 | 0.199 | 0.366 | 1 | 0.359 |
| PD | brainseq_caudate | M000 | True | 6 | 1 | 1 | 1 | 0.366 | 0.366 | 1 | 0.441 |
| PD | brainseq_caudate | M016 | False | 1 | 1 | 1 | 1 | 0.0704 | NA | 1 | 0.522 |
| LBD | gtex_caudate_bg | M029 | False | 1 | 1 | 1 | 1 | 0.0444 | NA | 1 | 0.256 |
| LBD | gtex_caudate_bg | M007 | False | 2 | 1 | 1 | 1 | 0.0879 | NA | 1 | 0.299 |
| LBD | brainseq_caudate | M012 | False | 1 | 1 | 1 | 1 | 0.0333 | NA | 0 | 0.144 |
| ALS | gtex_caudate_bg | M004 | True | 16 | 2 | 2 | 2 | 0.412 | 0.593 | 0.751 | 0.145 |
| ALS | gtex_caudate_bg | M000 | True | 29 | 2 | 2 | 2 | 0.763 | 0.784 | 0.931 | 0.54 |
| ALS | gtex_caudate_bg | M005 | True | 3 | 1 | 1 | 1 | 0.238 | 0.593 | 1 | 0.937 |
| ALS | gtex_caudate_bg | M001 | True | 5 | 1 | 1 | 1 | 0.367 | 0.593 | 1 | 0.0814 |
| ALS | gtex_caudate_bg | M007 | True | 5 | 1 | 1 | 0 | 0.367 | 0.593 | 1 | 0.412 |
| ALS | gtex_caudate_bg | M006 | True | 6 | 1 | 1 | 0 | 0.423 | 0.593 | 1 | 0.0131 |
| ALS | gtex_caudate_bg | M002 | True | 16 | 1 | 1 | 1 | 0.784 | 0.784 | 1 | 0.0543 |
| ALS | gtex_caudate_bg | M017 | False | 2 | 1 | 1 | 1 | 0.165 | NA | 1 | 0.796 |
| ALS | gtex_caudate_bg | M019 | False | 2 | 1 | 1 | 1 | 0.165 | NA | 1 | 0.173 |
| ALS | brainseq_caudate | M003 | True | 10 | 1 | 1 | 1 | 0.37 | 0.37 | 1 | 0.0228 |
| ALS | brainseq_caudate | M017 | False | 2 | 1 | 1 | 1 | 0.0833 | NA | 1 | 0.0579 |
| ALS | brainseq_caudate | M031 | False | 2 | 1 | 1 | 1 | 0.0833 | NA | 1 | 0.545 |
| SCZ | gtex_caudate_bg | M002 | True | 58 | 6 | 6 | 5 | 0.0915 | 0.672 | 0.191 | 0.0109 |
| SCZ | gtex_caudate_bg | M000 | True | 105 | 5 | 5 | 5 | 0.738 | 0.853 | 0.856 | 0.134 |
| SCZ | gtex_caudate_bg | M008 | True | 47 | 3 | 3 | 1 | 0.502 | 0.818 | 0.744 | 0.000748 |
| SCZ | gtex_caudate_bg | M009 | True | 20 | 2 | 2 | 1 | 0.31 | 0.672 | 0.678 | 0.437 |
| SCZ | gtex_caudate_bg | M004 | True | 57 | 2 | 2 | 2 | 0.853 | 0.853 | 0.965 | 0.0478 |
| SCZ | gtex_caudate_bg | M016 | True | 3 | 1 | 1 | 1 | 0.159 | 0.672 | 1 | 0.3 |
| SCZ | gtex_caudate_bg | M029 | True | 4 | 1 | 1 | 1 | 0.207 | 0.672 | 1 | 0.113 |
| SCZ | gtex_caudate_bg | M005 | True | 6 | 1 | 1 | 1 | 0.294 | 0.672 | 1 | 0.185 |
| SCZ | gtex_caudate_bg | M018 | True | 6 | 1 | 1 | 1 | 0.294 | 0.672 | 1 | 0.346 |
| SCZ | gtex_caudate_bg | M012 | True | 12 | 1 | 1 | 1 | 0.503 | 0.818 | 1 | 0.768 |
| SCZ | gtex_caudate_bg | M006 | True | 17 | 1 | 1 | 0 | 0.631 | 0.853 | 1 | 0.035 |
| SCZ | gtex_caudate_bg | M007 | True | 21 | 1 | 1 | 1 | 0.71 | 0.853 | 1 | 0.659 |
| SCZ | gtex_caudate_bg | M001 | True | 32 | 1 | 1 | 1 | 0.852 | 0.853 | 1 | 0.174 |
| SCZ | gtex_caudate_bg | M026 | False | 2 | 1 | 1 | 1 | 0.109 | NA | 1 | 0.146 |
| SCZ | brainseq_caudate | M003 | True | 42 | 6 | 6 | 6 | 0.0255 | 0.229 | 0.0694 | 0.447 |
| SCZ | brainseq_caudate | M004 | True | 25 | 2 | 2 | 2 | 0.449 | 0.673 | 0.778 | 0.0449 |
| SCZ | brainseq_caudate | M001 | True | 54 | 2 | 2 | 0 | 0.865 | 0.904 | 0.968 | 0.0168 |
| SCZ | brainseq_caudate | M022 | True | 3 | 1 | 1 | 1 | 0.168 | 0.673 | 1 | 0.704 |
| SCZ | brainseq_caudate | M011 | True | 7 | 1 | 1 | 1 | 0.352 | 0.673 | 1 | 0.0933 |
| SCZ | brainseq_caudate | M014 | True | 7 | 1 | 1 | 1 | 0.352 | 0.673 | 1 | 0.687 |
| SCZ | brainseq_caudate | M012 | True | 8 | 1 | 1 | 1 | 0.391 | 0.673 | 1 | 0.0125 |
| SCZ | brainseq_caudate | M009 | True | 14 | 1 | 1 | 0 | 0.584 | 0.751 | 1 | 5.8e-05 |
| SCZ | brainseq_caudate | M000 | True | 36 | 1 | 1 | 1 | 0.904 | 0.904 | 1 | 0.338 |
| SCZ | brainseq_caudate | M016 | False | 2 | 1 | 1 | 1 | 0.115 | NA | 1 | 0.958 |
| SCZ | brainseq_caudate | M028 | False | 2 | 1 | 1 | 1 | 0.115 | NA | 1 | 0.186 |

