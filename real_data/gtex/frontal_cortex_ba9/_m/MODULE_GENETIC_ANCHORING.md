# Module-level genetic anchoring — gtex-aging/frontal_cortex_ba9 vs GTEx Brain_Frontal_Cortex_BA9 xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 768 | False | 0.953 | 0.996 | -0.044 | -0.005 | 0.636 |
| M001 | 681 | False | 0.532 | 1.021 | -0.652 | 0.001 | 1.000 |
| M003 | 538 | True | 1.319 | 1.303 | 0.012 | -0.008 | 0.456 |
| M004 | 532 | False | 0.578 | 0.631 | -0.088 | 0.001 | 0.746 |
| M005 | 504 | False | 1.676 | 1.25 | 0.293 | 0.003 | 0.020 |
| M006 | 480 | True | 1.038 | 0.899 | 0.144 | -0.002 | 0.148 |
| M007 | 426 | False | 0.749 | 0.91 | -0.195 | -0.007 | 0.895 |
| M008 | 307 | False | 1.064 | 1.007 | 0.055 | -0.004 | 0.394 |
| M009 | 290 | True | 1.824 | 1.408 | 0.259 | -0.013 | 0.065 |
| M010 | 276 | True | 1.087 | 1.389 | -0.246 | -0.025 | 0.868 |
| M011 | 176 | True | 0.72 | 0.715 | 0.007 | -0.004 | 0.495 |
| M012 | 168 | False | 0.695 | 0.899 | -0.258 | -0.024 | 0.834 |
| M013 | 147 | True | 0.832 | 1.057 | -0.239 | -0.017 | 0.796 |
| M014 | 146 | True | 1.828 | 1.753 | 0.042 | -0.016 | 0.426 |
| M015 | 142 | True | 0.99 | 1.137 | -0.139 | 0.0 | 0.702 |
| M016 | 141 | False | 0.881 | 0.913 | -0.036 | -0.006 | 0.552 |
| M017 | 135 | True | 1.006 | 0.845 | 0.174 | -0.014 | 0.239 |
| M018 | 125 | True | 0.89 | 1.228 | -0.322 | -0.006 | 0.865 |
| M019 | 125 | True | 0.828 | 0.954 | -0.141 | -0.018 | 0.667 |
| M021 | 101 | True | 1.534 | 1.508 | 0.017 | -0.023 | 0.455 |
| M023 | 90 | False | 1.017 | 0.91 | 0.112 | -0.028 | 0.360 |
| M024 | 87 | False | 0.514 | 1.015 | -0.68 | -0.003 | 0.959 |
| M025 | 83 | False | 0.923 | 0.617 | 0.403 | -0.031 | 0.098 |
| M027 | 78 | True | 1.463 | 1.014 | 0.366 | -0.022 | 0.134 |
| M028 | 77 | True | 0.681 | 1.1 | -0.48 | -0.026 | 0.896 |
| M029 | 70 | True | 0.812 | 0.913 | -0.117 | -0.007 | 0.629 |
| M030 | 68 | True | 1.089 | 1.029 | 0.056 | -0.025 | 0.429 |
| M031 | 64 | True | 1.36 | 1.851 | -0.308 | -0.052 | 0.748 |
| M032 | 60 | True | 0.925 | 1.108 | -0.181 | -0.02 | 0.664 |
| M033 | 58 | False | 0.167 | 0.635 | -1.334 | -0.011 | 0.998 |
| M034 | 51 | True | 2.118 | 1.525 | 0.328 | -0.1 | 0.197 |
| M035 | 47 | True | 1.69 | 1.357 | 0.219 | -0.057 | 0.293 |
| M036 | 45 | True | 0.673 | 0.564 | 0.176 | -0.065 | 0.305 |
| M037 | 39 | True | 0.915 | 0.938 | -0.025 | -0.047 | 0.504 |
| M038 | 33 | True | 1.178 | 1.25 | -0.059 | -0.103 | 0.529 |
| M039 | 33 | True | 0.737 | 1.299 | -0.567 | -0.107 | 0.825 |
| M040 | 28 | True | 0.999 | 1.488 | -0.399 | -0.059 | 0.714 |
| M041 | 27 | True | 0.772 | 1.455 | -0.634 | -0.236 | 0.799 |
| M042 | 26 | True | 0.796 | 1.227 | -0.432 | -0.16 | 0.688 |
| M043 | 25 | True | 0.771 | 0.872 | -0.123 | -0.214 | 0.528 |
| M044 | 25 | True | 0.599 | 0.635 | -0.059 | -0.1 | 0.495 |
| M045 | 23 | True | 0.252 | 0.831 | -1.194 | -0.586 | 0.933 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._