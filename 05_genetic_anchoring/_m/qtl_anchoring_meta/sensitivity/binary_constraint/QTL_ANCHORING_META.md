# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph, wgcna_switch_only, wgcna_multiplex.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 101117.0 | 0.92 | 0.9 | 0.94 | 5.68e-19 | 0.92 | 2.75e-15 | 0.21 |
| isograph | all_modules | sQTL | 17 | 91823.0 | 0.95 | 0.92 | 0.97 | 5.60e-05 | 0.94 | 3.25e-03 | 0.52 |
| isograph | pheno_sig_modules | eQTL | 10 | 18088.0 | 0.93 | 0.89 | 0.96 | 2.03e-05 | 0.91 | 8.14e-03 | 0.68 |
| isograph | pheno_sig_modules | sQTL | 9 | 16210.0 | 0.97 | 0.93 | 1.02 | 2.64e-01 | 0.95 | 3.14e-01 | 0.77 |
| isograph | go_invisible_modules | eQTL | 7 | 4574.0 | 0.99 | 0.93 | 1.05 | 6.38e-01 | 0.99 | 6.38e-01 | 0.0 |
| isograph | go_invisible_modules | sQTL | 7 | 4307.0 | 1.04 | 0.96 | 1.12 | 3.75e-01 | 1.03 | 4.81e-01 | 0.12 |
| isograph | go_visible_modules | eQTL | 9 | 13514.0 | 0.91 | 0.87 | 0.95 | 4.99e-06 | 0.88 | 6.66e-03 | 0.73 |
| isograph | go_visible_modules | sQTL | 8 | 11903.0 | 0.95 | 0.9 | 1.0 | 6.13e-02 | 0.97 | 7.21e-01 | 0.82 |
| isograph | all_modules | eQTL | 17 | 101117.0 | 0.89 | 0.88 | 0.91 | 9.36e-38 | 0.89 | 1.25e-08 | 0.8 |
| isograph | all_modules | sQTL | 17 | 91823.0 | 0.97 | 0.95 | 1.0 | 4.58e-02 | 0.97 | 1.24e-01 | 0.58 |
| isograph | pheno_sig_modules | eQTL | 10 | 18088.0 | 0.89 | 0.86 | 0.92 | 2.01e-11 | 0.87 | 2.06e-03 | 0.81 |
| isograph | pheno_sig_modules | sQTL | 10 | 16256.0 | 0.99 | 0.94 | 1.03 | 5.78e-01 | 0.95 | 4.14e-01 | 0.78 |
| isograph | go_invisible_modules | eQTL | 7 | 4574.0 | 0.99 | 0.93 | 1.06 | 8.13e-01 | 0.99 | 7.63e-01 | 0.3 |
| isograph | go_invisible_modules | sQTL | 7 | 4307.0 | 1.06 | 0.98 | 1.14 | 1.64e-01 | 1.04 | 4.96e-01 | 0.4 |
| isograph | go_visible_modules | eQTL | 9 | 13514.0 | 0.86 | 0.83 | 0.9 | 6.32e-14 | 0.84 | 5.31e-04 | 0.79 |
| isograph | go_visible_modules | sQTL | 9 | 11949.0 | 0.96 | 0.91 | 1.01 | 1.21e-01 | 0.98 | 7.43e-01 | 0.81 |
| wgcna_switch_only | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | pheno_sig_modules | eQTL | 6 | 19315.0 | 1.06 | 1.0 | 1.12 | 3.37e-02 | 1.06 | 3.37e-02 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | sQTL | 6 | 18235.0 | 0.99 | 0.92 | 1.06 | 6.75e-01 | 0.98 | 6.91e-01 | 0.59 |
| wgcna_switch_only | go_invisible_modules | eQTL | 6 | 17935.0 | 1.06 | 1.01 | 1.12 | 3.02e-02 | 1.06 | 3.02e-02 | 0.0 |
| wgcna_switch_only | go_invisible_modules | sQTL | 6 | 16934.0 | 0.99 | 0.92 | 1.06 | 7.57e-01 | 0.98 | 7.44e-01 | 0.54 |
| wgcna_switch_only | go_visible_modules | eQTL | 2 | 1380.0 | 1.0 | 0.87 | 1.14 | 9.92e-01 | 1.0 | 9.92e-01 | 0.0 |
| wgcna_switch_only | go_visible_modules | sQTL | 2 | 1301.0 | 0.97 | 0.81 | 1.18 | 7.81e-01 | 0.97 | 7.81e-01 | 0.0 |
| wgcna_switch_only | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | pheno_sig_modules | eQTL | 6 | 19315.0 | 1.04 | 0.99 | 1.1 | 1.18e-01 | 1.04 | 1.20e-01 | 0.01 |
| wgcna_switch_only | pheno_sig_modules | sQTL | 6 | 18235.0 | 0.95 | 0.89 | 1.02 | 1.49e-01 | 0.94 | 4.13e-01 | 0.78 |
| wgcna_switch_only | go_invisible_modules | eQTL | 6 | 17935.0 | 1.04 | 0.99 | 1.1 | 1.33e-01 | 1.04 | 1.39e-01 | 0.03 |
| wgcna_switch_only | go_invisible_modules | sQTL | 6 | 16934.0 | 0.95 | 0.89 | 1.02 | 1.53e-01 | 0.94 | 4.06e-01 | 0.76 |
| wgcna_switch_only | go_visible_modules | eQTL | 2 | 1380.0 | 1.02 | 0.89 | 1.16 | 8.23e-01 | 1.02 | 8.23e-01 | 0.0 |
| wgcna_switch_only | go_visible_modules | sQTL | 2 | 1301.0 | 1.0 | 0.83 | 1.2 | 9.65e-01 | 1.0 | 9.65e-01 | 0.0 |
| wgcna_multiplex | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | pheno_sig_modules | eQTL | 8 | 28771.0 | 0.92 | 0.88 | 0.95 | 2.72e-06 | 0.87 | 1.10e-01 | 0.95 |
| wgcna_multiplex | pheno_sig_modules | sQTL | 8 | 24500.0 | 1.02 | 0.96 | 1.07 | 5.50e-01 | 1.01 | 8.80e-01 | 0.79 |
| wgcna_multiplex | go_invisible_modules | eQTL | 4 | 4938.0 | 0.96 | 0.9 | 1.03 | 2.45e-01 | 0.86 | 1.05e-01 | 0.79 |
| wgcna_multiplex | go_invisible_modules | sQTL | 4 | 4555.0 | 0.94 | 0.86 | 1.02 | 1.55e-01 | 0.97 | 7.03e-01 | 0.6 |
| wgcna_multiplex | go_visible_modules | eQTL | 8 | 26650.0 | 0.93 | 0.89 | 0.96 | 4.22e-05 | 0.88 | 1.38e-01 | 0.95 |
| wgcna_multiplex | go_visible_modules | sQTL | 8 | 22610.0 | 1.04 | 0.99 | 1.1 | 1.21e-01 | 1.02 | 7.89e-01 | 0.84 |
| wgcna_multiplex | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | pheno_sig_modules | eQTL | 8 | 28771.0 | 0.9 | 0.87 | 0.94 | 2.09e-08 | 0.86 | 7.50e-02 | 0.95 |
| wgcna_multiplex | pheno_sig_modules | sQTL | 8 | 24500.0 | 1.02 | 0.97 | 1.08 | 3.47e-01 | 1.01 | 8.81e-01 | 0.87 |
| wgcna_multiplex | go_invisible_modules | eQTL | 4 | 4938.0 | 0.92 | 0.86 | 0.98 | 9.14e-03 | 0.83 | 4.66e-02 | 0.78 |
| wgcna_multiplex | go_invisible_modules | sQTL | 4 | 4555.0 | 0.94 | 0.86 | 1.02 | 1.55e-01 | 0.97 | 8.11e-01 | 0.71 |
| wgcna_multiplex | go_visible_modules | eQTL | 8 | 26650.0 | 0.92 | 0.89 | 0.95 | 2.49e-06 | 0.87 | 1.15e-01 | 0.95 |
| wgcna_multiplex | go_visible_modules | sQTL | 8 | 22610.0 | 1.05 | 1.0 | 1.11 | 3.87e-02 | 1.02 | 8.10e-01 | 0.91 |

## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 17 | 1.028 | 0.995 | 1.061 | 9.66e-02 | 1.027 | 0.459 |
| isograph | pheno_sig_modules | 9 | 1.039 | 0.978 | 1.103 | 2.19e-01 | 1.038 | 0.243 |
| isograph | go_invisible_modules | 7 | 1.043 | 0.943 | 1.154 | 4.12e-01 | 1.043 | 0.0 |
| isograph | go_visible_modules | 8 | 1.03 | 0.961 | 1.105 | 4.04e-01 | 1.069 | 0.609 |
| isograph | all_modules | 17 | 1.088 | 1.055 | 1.122 | 6.25e-08 | 1.085 | 0.619 |
| isograph | pheno_sig_modules | 10 | 1.09 | 1.028 | 1.155 | 3.87e-03 | 1.09 | 0.0 |
| isograph | go_invisible_modules | 7 | 1.051 | 0.951 | 1.16 | 3.31e-01 | 1.051 | 0.0 |
| isograph | go_visible_modules | 9 | 1.095 | 1.023 | 1.171 | 8.41e-03 | 1.115 | 0.417 |
| wgcna_switch_only | pheno_sig_modules | 6 | 0.931 | 0.853 | 1.016 | 1.11e-01 | 0.921 | 0.57 |
| wgcna_switch_only | go_invisible_modules | 6 | 0.933 | 0.853 | 1.02 | 1.25e-01 | 0.924 | 0.536 |
| wgcna_switch_only | go_visible_modules | 2 | 0.973 | 0.772 | 1.226 | 8.16e-01 | 0.973 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | 6 | 0.918 | 0.842 | 1.0 | 5.00e-02 | 0.908 | 0.599 |
| wgcna_switch_only | go_invisible_modules | 6 | 0.918 | 0.841 | 1.001 | 5.36e-02 | 0.909 | 0.557 |
| wgcna_switch_only | go_visible_modules | 2 | 0.979 | 0.78 | 1.229 | 8.55e-01 | 0.979 | 0.0 |
| wgcna_multiplex | pheno_sig_modules | 8 | 1.087 | 1.019 | 1.159 | 1.09e-02 | 1.114 | 0.677 |
| wgcna_multiplex | go_invisible_modules | 4 | 0.973 | 0.871 | 1.087 | 6.29e-01 | 1.111 | 0.73 |
| wgcna_multiplex | go_visible_modules | 8 | 1.098 | 1.031 | 1.169 | 3.83e-03 | 1.107 | 0.46 |
| wgcna_multiplex | pheno_sig_modules | 8 | 1.109 | 1.042 | 1.18 | 1.08e-03 | 1.124 | 0.538 |
| wgcna_multiplex | go_invisible_modules | 4 | 1.02 | 0.916 | 1.136 | 7.16e-01 | 1.14 | 0.689 |
| wgcna_multiplex | go_visible_modules | 8 | 1.119 | 1.053 | 1.19 | 3.22e-04 | 1.121 | 0.371 |

## Cross-method contrast on the 13 tissues all methods share

The matched WGCNA baselines (`wgcna_switch_only`, `wgcna_multiplex`) consume the SAME switch / switch+abundance features as IsoGraph, so comparing their splicing-specificity on the SAME tissues isolates whether the signal lives in the switch features or in IsoGraph's VAE + Leiden inference.

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 5 | 1.019 | 0.961 | 1.08 | 5.32e-01 | 1.016 | 0.624 |
| isograph | all_modules | 5 | 1.077 | 1.019 | 1.138 | 8.42e-03 | 1.075 | 0.634 |
| isograph | pheno_sig_modules | 5 | 1.018 | 0.95 | 1.091 | 6.16e-01 | 1.017 | 0.147 |
| isograph | pheno_sig_modules | 5 | 1.077 | 1.008 | 1.151 | 2.91e-02 | 1.077 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | 5 | 0.936 | 0.851 | 1.029 | 1.71e-01 | 0.922 | 0.654 |
| wgcna_switch_only | pheno_sig_modules | 5 | 0.924 | 0.842 | 1.014 | 9.75e-02 | 0.911 | 0.675 |
| wgcna_multiplex | pheno_sig_modules | 5 | 1.059 | 0.988 | 1.134 | 1.07e-01 | 1.059 | 0.744 |
| wgcna_multiplex | pheno_sig_modules | 5 | 1.084 | 1.014 | 1.159 | 1.82e-02 | 1.082 | 0.461 |
| isograph | go_invisible_modules | 4 | 1.057 | 0.921 | 1.213 | 4.30e-01 | 1.057 | 0.0 |
| isograph | go_invisible_modules | 4 | 1.046 | 0.914 | 1.197 | 5.10e-01 | 1.046 | 0.0 |
| wgcna_switch_only | go_invisible_modules | 5 | 0.931 | 0.847 | 1.025 | 1.45e-01 | 0.919 | 0.629 |
| wgcna_switch_only | go_invisible_modules | 5 | 0.919 | 0.837 | 1.009 | 7.69e-02 | 0.908 | 0.645 |
| wgcna_multiplex | go_invisible_modules | 3 | 0.958 | 0.858 | 1.071 | 4.55e-01 | 1.013 | 0.659 |
| wgcna_multiplex | go_invisible_modules | 3 | 1.006 | 0.903 | 1.122 | 9.09e-01 | 1.054 | 0.611 |
| isograph | go_visible_modules | 5 | 1.002 | 0.93 | 1.08 | 9.56e-01 | 1.003 | 0.068 |
| isograph | go_visible_modules | 5 | 1.072 | 0.998 | 1.151 | 5.55e-02 | 1.072 | 0.0 |
| wgcna_switch_only | go_visible_modules | 1 | 1.088 | 0.696 | 1.702 | 7.11e-01 | 1.088 | 0.0 |
| wgcna_switch_only | go_visible_modules | 1 | 1.106 | 0.712 | 1.717 | 6.53e-01 | 1.106 | 0.0 |
| wgcna_multiplex | go_visible_modules | 5 | 1.078 | 1.007 | 1.154 | 2.99e-02 | 1.074 | 0.55 |
| wgcna_multiplex | go_visible_modules | 5 | 1.102 | 1.032 | 1.177 | 3.91e-03 | 1.097 | 0.284 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, carried by the phenotype-associated modules — splicing-QTL is spared relative to expression-QTL in the modules that track the trait.
- **Do NOT read a GO-invisible localisation off this table.** On the 2026-08-29 refresh the `go_invisible_modules` and `go_visible_modules` arms are statistically indistinguishable, so neither is the site of the effect and neither is a null. The DTU-without-DGE claim rests on the GO-invisible content gate, not on these genetics. An earlier framing made `go_invisible` the headline; it came from a stale `module_enrichment` join that relabelled the partition and must not be restored.
- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL signal is a property of the switch features (which both methods share), not of IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. IsoGraph's value is the switch-feature *representation* and its finer modules, not a unique network-inference effect.
- **Primary internal control = the matched WGCNA baselines**, not `go_visible_modules`. The baselines hold the switch features fixed and vary only the inference, so a null there localises the effect to IsoGraph's inference. `go_visible_modules` is a secondary, exploratory stratification on module CONTENT and is only ever a relative contrast — never an on/off null.
- High I2 flags between-tissue heterogeneity; prefer RE there. A nominally significant ratio carrying high I2 is driven by a few tissues, not by a consistent effect, and is weaker evidence than a smaller ratio at I2 near 0 — compare module sets on consistency as well as magnitude. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
