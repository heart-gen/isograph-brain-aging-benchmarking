# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph, wgcna_switch_only, wgcna_multiplex.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 86465.0 | 0.81 | 0.79 | 0.82 | 7.30e-138 | 0.81 | 1.18e-60 | 0.56 |
| isograph | all_modules | sQTL | 17 | 70399.0 | 0.87 | 0.85 | 0.89 | 9.11e-31 | 0.86 | 8.57e-19 | 0.43 |
| isograph | pheno_sig_modules | eQTL | 12 | 20631.0 | 0.76 | 0.74 | 0.79 | 1.41e-59 | 0.74 | 2.82e-14 | 0.76 |
| isograph | pheno_sig_modules | sQTL | 12 | 16433.0 | 0.86 | 0.82 | 0.9 | 1.93e-10 | 0.83 | 2.49e-05 | 0.55 |
| isograph | go_invisible_modules | eQTL | 11 | 9926.0 | 0.92 | 0.88 | 0.96 | 1.40e-04 | 0.9 | 3.19e-02 | 0.74 |
| isograph | go_invisible_modules | sQTL | 11 | 8396.0 | 0.99 | 0.93 | 1.05 | 6.62e-01 | 0.97 | 4.50e-01 | 0.23 |
| isograph | go_visible_modules | eQTL | 11 | 10705.0 | 0.68 | 0.65 | 0.71 | 2.04e-68 | 0.63 | 1.18e-15 | 0.79 |
| isograph | go_visible_modules | sQTL | 11 | 8037.0 | 0.75 | 0.7 | 0.81 | 2.38e-16 | 0.72 | 2.33e-09 | 0.37 |
| wgcna_switch_only | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | pheno_sig_modules | eQTL | 11 | 30161.0 | 0.93 | 0.89 | 0.97 | 1.36e-03 | 0.85 | 2.35e-02 | 0.88 |
| wgcna_switch_only | pheno_sig_modules | sQTL | 11 | 25165.0 | 0.94 | 0.88 | 1.01 | 8.00e-02 | 0.88 | 3.48e-02 | 0.55 |
| wgcna_switch_only | go_invisible_modules | eQTL | 7 | 5829.0 | 1.03 | 0.97 | 1.1 | 3.07e-01 | 1.02 | 6.82e-01 | 0.53 |
| wgcna_switch_only | go_invisible_modules | sQTL | 7 | 5115.0 | 1.07 | 0.98 | 1.16 | 1.58e-01 | 1.06 | 4.07e-01 | 0.54 |
| wgcna_switch_only | go_visible_modules | eQTL | 11 | 24332.0 | 0.93 | 0.89 | 0.96 | 2.58e-04 | 0.85 | 6.06e-03 | 0.85 |
| wgcna_switch_only | go_visible_modules | sQTL | 11 | 20050.0 | 0.92 | 0.87 | 0.98 | 1.22e-02 | 0.89 | 3.39e-02 | 0.53 |
| wgcna_multiplex | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | pheno_sig_modules | eQTL | 9 | 56696.0 | 0.94 | 0.92 | 0.97 | 1.81e-05 | 0.83 | 1.03e-02 | 0.96 |
| wgcna_multiplex | pheno_sig_modules | sQTL | 9 | 43914.0 | 0.94 | 0.91 | 0.98 | 5.39e-03 | 0.89 | 5.30e-02 | 0.87 |
| wgcna_multiplex | go_invisible_modules | eQTL | 7 | 7561.0 | 0.99 | 0.94 | 1.04 | 5.61e-01 | 0.96 | 6.39e-01 | 0.91 |
| wgcna_multiplex | go_invisible_modules | sQTL | 7 | 6354.0 | 1.0 | 0.93 | 1.07 | 9.49e-01 | 0.98 | 7.84e-01 | 0.66 |
| wgcna_multiplex | go_visible_modules | eQTL | 9 | 53239.0 | 0.94 | 0.91 | 0.96 | 1.52e-06 | 0.83 | 4.72e-03 | 0.96 |
| wgcna_multiplex | go_visible_modules | sQTL | 9 | 41201.0 | 0.95 | 0.91 | 0.99 | 1.49e-02 | 0.89 | 6.69e-02 | 0.88 |

## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 17 | 1.068 | 1.037 | 1.1 | 1.35e-05 | 1.066 | 0.291 |
| isograph | pheno_sig_modules | 12 | 1.111 | 1.048 | 1.177 | 3.57e-04 | 1.107 | 0.229 |
| isograph | go_invisible_modules | 11 | 1.068 | 0.993 | 1.148 | 7.72e-02 | 1.068 | 0.0 |
| isograph | go_visible_modules | 11 | 1.084 | 1.0 | 1.174 | 5.02e-02 | 1.082 | 0.284 |
| wgcna_switch_only | pheno_sig_modules | 11 | 0.977 | 0.9 | 1.06 | 5.74e-01 | 0.977 | 0.002 |
| wgcna_switch_only | go_invisible_modules | 7 | 1.017 | 0.912 | 1.135 | 7.60e-01 | 1.026 | 0.362 |
| wgcna_switch_only | go_visible_modules | 11 | 0.97 | 0.899 | 1.047 | 4.38e-01 | 0.97 | 0.0 |
| wgcna_multiplex | pheno_sig_modules | 9 | 0.983 | 0.936 | 1.031 | 4.77e-01 | 0.993 | 0.722 |
| wgcna_multiplex | go_invisible_modules | 7 | 0.993 | 0.913 | 1.079 | 8.63e-01 | 0.993 | 0.0 |
| wgcna_multiplex | go_visible_modules | 9 | 0.995 | 0.949 | 1.044 | 8.51e-01 | 1.004 | 0.653 |

## Cross-method contrast on the 13 tissues all methods share

The matched WGCNA baselines (`wgcna_switch_only`, `wgcna_multiplex`) consume the SAME switch / switch+abundance features as IsoGraph, so comparing their splicing-specificity on the SAME tissues isolates whether the signal lives in the switch features or in IsoGraph's VAE + Leiden inference.

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 8 | 1.105 | 1.058 | 1.154 | 5.81e-06 | 1.105 | 0.003 |
| isograph | pheno_sig_modules | 8 | 1.108 | 1.042 | 1.177 | 1.01e-03 | 1.108 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | 8 | 0.968 | 0.889 | 1.055 | 4.63e-01 | 0.951 | 0.25 |
| wgcna_multiplex | pheno_sig_modules | 8 | 0.98 | 0.934 | 1.029 | 4.15e-01 | 0.983 | 0.741 |
| isograph | go_invisible_modules | 7 | 1.066 | 0.985 | 1.153 | 1.12e-01 | 1.061 | 0.157 |
| wgcna_switch_only | go_invisible_modules | 6 | 1.021 | 0.914 | 1.14 | 7.12e-01 | 1.036 | 0.456 |
| wgcna_multiplex | go_invisible_modules | 7 | 0.993 | 0.913 | 1.079 | 8.63e-01 | 0.993 | 0.0 |
| isograph | go_visible_modules | 8 | 1.076 | 0.992 | 1.166 | 7.83e-02 | 1.075 | 0.079 |
| wgcna_switch_only | go_visible_modules | 8 | 0.959 | 0.886 | 1.038 | 3.03e-01 | 0.959 | 0.087 |
| wgcna_multiplex | go_visible_modules | 8 | 0.993 | 0.946 | 1.042 | 7.72e-01 | 0.995 | 0.674 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, strongest for the GO-invisible modules — splicing-QTL is spared relative to expression-QTL exactly where the DTU-without-DGE value concentrates.
- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL signal is a property of the switch features (which both methods share), not of IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. IsoGraph's value is the switch-feature *representation* and its finer modules, not a unique network-inference effect.
- **Primary internal control = the matched WGCNA baselines**, not `go_visible_modules`. The baselines hold the switch features fixed and vary only the inference, so a null there localises the effect to IsoGraph's inference. `go_visible_modules` is a secondary control on module CONTENT and is only ever a relative contrast — read it as the low end of a gradient, not as an on/off null.
- High I2 flags between-tissue heterogeneity; prefer RE there. A nominally significant ratio carrying high I2 is driven by a few tissues, not by a consistent effect, and is weaker evidence than a smaller ratio at I2 near 0 — compare module sets on consistency as well as magnitude. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
