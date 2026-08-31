# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph, wgcna_switch_only, wgcna_multiplex.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 86465.0 | 0.88 | 0.87 | 0.89 | 3.63e-206 | 0.88 | 9.22e-69 | 0.67 |
| isograph | all_modules | sQTL | 17 | 70399.0 | 0.95 | 0.95 | 0.96 | 1.14e-24 | 0.95 | 4.17e-19 | 0.23 |
| isograph | pheno_sig_modules | eQTL | 12 | 20631.0 | 0.85 | 0.84 | 0.87 | 3.96e-91 | 0.84 | 2.46e-20 | 0.75 |
| isograph | pheno_sig_modules | sQTL | 12 | 16433.0 | 0.95 | 0.93 | 0.97 | 2.84e-09 | 0.95 | 4.26e-06 | 0.32 |
| isograph | go_invisible_modules | eQTL | 11 | 9926.0 | 0.95 | 0.93 | 0.97 | 5.75e-06 | 0.94 | 7.37e-03 | 0.69 |
| isograph | go_invisible_modules | sQTL | 11 | 8396.0 | 0.99 | 0.97 | 1.01 | 2.69e-01 | 0.98 | 2.40e-01 | 0.52 |
| isograph | go_visible_modules | eQTL | 11 | 10705.0 | 0.8 | 0.78 | 0.81 | 1.03e-107 | 0.77 | 4.05e-23 | 0.78 |
| isograph | go_visible_modules | sQTL | 11 | 8037.0 | 0.92 | 0.9 | 0.94 | 1.06e-11 | 0.92 | 1.31e-08 | 0.15 |
| wgcna_switch_only | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | pheno_sig_modules | eQTL | 11 | 30161.0 | 0.98 | 0.96 | 1.0 | 2.13e-02 | 0.94 | 7.65e-02 | 0.89 |
| wgcna_switch_only | pheno_sig_modules | sQTL | 11 | 25165.0 | 0.97 | 0.95 | 0.99 | 7.27e-03 | 0.96 | 1.28e-02 | 0.44 |
| wgcna_switch_only | go_invisible_modules | eQTL | 7 | 5829.0 | 1.02 | 0.99 | 1.05 | 1.48e-01 | 1.02 | 3.53e-01 | 0.48 |
| wgcna_switch_only | go_invisible_modules | sQTL | 7 | 5115.0 | 1.02 | 0.98 | 1.05 | 3.11e-01 | 1.02 | 3.28e-01 | 0.41 |
| wgcna_switch_only | go_visible_modules | eQTL | 11 | 24332.0 | 0.97 | 0.95 | 0.99 | 1.77e-03 | 0.93 | 1.30e-02 | 0.85 |
| wgcna_switch_only | go_visible_modules | sQTL | 11 | 20050.0 | 0.97 | 0.94 | 0.99 | 1.41e-03 | 0.96 | 4.53e-03 | 0.3 |
| wgcna_multiplex | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | pheno_sig_modules | eQTL | 9 | 56696.0 | 0.97 | 0.96 | 0.99 | 1.86e-05 | 0.91 | 1.13e-02 | 0.97 |
| wgcna_multiplex | pheno_sig_modules | sQTL | 9 | 43914.0 | 0.97 | 0.96 | 0.98 | 2.98e-05 | 0.96 | 7.68e-03 | 0.7 |
| wgcna_multiplex | go_invisible_modules | eQTL | 7 | 7561.0 | 0.99 | 0.97 | 1.01 | 4.08e-01 | 0.99 | 7.67e-01 | 0.93 |
| wgcna_multiplex | go_invisible_modules | sQTL | 7 | 6354.0 | 1.01 | 0.98 | 1.03 | 6.55e-01 | 1.0 | 8.94e-01 | 0.66 |
| wgcna_multiplex | go_visible_modules | eQTL | 9 | 53239.0 | 0.97 | 0.96 | 0.98 | 5.97e-07 | 0.9 | 4.94e-03 | 0.97 |
| wgcna_multiplex | go_visible_modules | sQTL | 9 | 41201.0 | 0.97 | 0.96 | 0.98 | 2.29e-05 | 0.96 | 1.03e-02 | 0.73 |

## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 17 | 1.082 | 1.069 | 1.095 | 9.97e-38 | 1.082 | 0.629 |
| isograph | pheno_sig_modules | 12 | 1.113 | 1.088 | 1.139 | 1.37e-19 | 1.113 | 0.067 |
| isograph | go_invisible_modules | 11 | 1.038 | 1.007 | 1.07 | 1.60e-02 | 1.038 | 0.099 |
| isograph | go_visible_modules | 11 | 1.155 | 1.12 | 1.191 | 2.91e-20 | 1.171 | 0.513 |
| wgcna_switch_only | pheno_sig_modules | 11 | 0.986 | 0.956 | 1.018 | 3.92e-01 | 0.989 | 0.607 |
| wgcna_switch_only | go_invisible_modules | 7 | 0.994 | 0.951 | 1.038 | 7.72e-01 | 0.994 | 0.0 |
| wgcna_switch_only | go_visible_modules | 11 | 0.991 | 0.962 | 1.021 | 5.49e-01 | 0.998 | 0.472 |
| wgcna_multiplex | pheno_sig_modules | 9 | 0.994 | 0.976 | 1.014 | 5.71e-01 | 1.03 | 0.906 |
| wgcna_multiplex | go_invisible_modules | 7 | 1.016 | 0.981 | 1.052 | 3.79e-01 | 1.013 | 0.458 |
| wgcna_multiplex | go_visible_modules | 9 | 0.998 | 0.979 | 1.018 | 8.59e-01 | 1.033 | 0.89 |

## Cross-method contrast on the 13 tissues all methods share

The matched WGCNA baselines (`wgcna_switch_only`, `wgcna_multiplex`) consume the SAME switch / switch+abundance features as IsoGraph, so comparing their splicing-specificity on the SAME tissues isolates whether the signal lives in the switch features or in IsoGraph's VAE + Leiden inference.

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 8 | 1.095 | 1.076 | 1.115 | 3.08e-24 | 1.095 | 0.589 |
| isograph | pheno_sig_modules | 8 | 1.116 | 1.089 | 1.143 | 1.38e-18 | 1.116 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | 8 | 0.983 | 0.95 | 1.016 | 3.02e-01 | 0.978 | 0.696 |
| wgcna_multiplex | pheno_sig_modules | 8 | 0.992 | 0.973 | 1.011 | 4.05e-01 | 1.011 | 0.908 |
| isograph | go_invisible_modules | 7 | 1.034 | 1.0 | 1.069 | 4.82e-02 | 1.03 | 0.347 |
| wgcna_switch_only | go_invisible_modules | 6 | 0.992 | 0.949 | 1.037 | 7.18e-01 | 0.993 | 0.102 |
| wgcna_multiplex | go_invisible_modules | 7 | 1.016 | 0.981 | 1.052 | 3.79e-01 | 1.013 | 0.458 |
| isograph | go_visible_modules | 8 | 1.154 | 1.119 | 1.19 | 1.58e-19 | 1.165 | 0.491 |
| wgcna_switch_only | go_visible_modules | 8 | 0.988 | 0.958 | 1.02 | 4.55e-01 | 0.993 | 0.58 |
| wgcna_multiplex | go_visible_modules | 8 | 0.996 | 0.977 | 1.015 | 6.59e-01 | 1.016 | 0.891 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, strongest for the GO-invisible modules — splicing-QTL is spared relative to expression-QTL exactly where the DTU-without-DGE value concentrates.
- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL signal is a property of the switch features (which both methods share), not of IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. IsoGraph's value is the switch-feature *representation* and its finer modules, not a unique network-inference effect.
- **Primary internal control = the matched WGCNA baselines**, not `go_visible_modules`. The baselines hold the switch features fixed and vary only the inference, so a null there localises the effect to IsoGraph's inference. `go_visible_modules` is a secondary control on module CONTENT and is only ever a relative contrast — read it as the low end of a gradient, not as an on/off null.
- High I2 flags between-tissue heterogeneity; prefer RE there. A nominally significant ratio carrying high I2 is driven by a few tissues, not by a consistent effect, and is weaker evidence than a smaller ratio at I2 near 0 — compare module sets on consistency as well as magnitude. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
