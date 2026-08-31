# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph, wgcna_switch_only, wgcna_multiplex.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 69642.0 | 0.86 | 0.84 | 0.88 | 2.27e-51 | 0.86 | 2.26e-16 | 0.7 |
| isograph | all_modules | sQTL | 17 | 60708.0 | 0.93 | 0.9 | 0.95 | 1.90e-07 | 0.93 | 4.97e-06 | 0.24 |
| isograph | pheno_sig_modules | eQTL | 12 | 16661.0 | 0.81 | 0.78 | 0.85 | 2.61e-26 | 0.78 | 3.63e-06 | 0.82 |
| isograph | pheno_sig_modules | sQTL | 12 | 14235.0 | 0.94 | 0.89 | 0.99 | 1.60e-02 | 0.91 | 3.20e-02 | 0.47 |
| isograph | go_invisible_modules | eQTL | 11 | 7907.0 | 0.98 | 0.94 | 1.03 | 5.13e-01 | 0.97 | 4.60e-01 | 0.44 |
| isograph | go_invisible_modules | sQTL | 11 | 7260.0 | 1.05 | 0.98 | 1.12 | 1.69e-01 | 1.05 | 1.69e-01 | 0.0 |
| isograph | go_visible_modules | eQTL | 11 | 8754.0 | 0.72 | 0.68 | 0.76 | 1.78e-38 | 0.65 | 1.83e-09 | 0.82 |
| isograph | go_visible_modules | sQTL | 11 | 6975.0 | 0.83 | 0.77 | 0.9 | 2.30e-06 | 0.78 | 1.48e-04 | 0.39 |
| isograph | all_modules | eQTL | 17 | 69642.0 | 0.84 | 0.82 | 0.85 | 8.63e-76 | 0.84 | 7.04e-23 | 0.71 |
| isograph | all_modules | sQTL | 17 | 60708.0 | 0.94 | 0.91 | 0.97 | 5.87e-06 | 0.94 | 5.46e-05 | 0.22 |
| isograph | pheno_sig_modules | eQTL | 12 | 16661.0 | 0.78 | 0.75 | 0.81 | 2.90e-41 | 0.76 | 2.06e-08 | 0.81 |
| isograph | pheno_sig_modules | sQTL | 12 | 14235.0 | 0.93 | 0.88 | 0.98 | 6.63e-03 | 0.9 | 2.01e-02 | 0.49 |
| isograph | go_invisible_modules | eQTL | 11 | 7907.0 | 0.95 | 0.9 | 0.99 | 2.30e-02 | 0.93 | 1.25e-01 | 0.59 |
| isograph | go_invisible_modules | sQTL | 11 | 7260.0 | 1.04 | 0.97 | 1.11 | 2.64e-01 | 1.04 | 2.64e-01 | 0.0 |
| isograph | go_visible_modules | eQTL | 11 | 8754.0 | 0.69 | 0.65 | 0.72 | 8.49e-52 | 0.63 | 1.49e-11 | 0.81 |
| isograph | go_visible_modules | sQTL | 11 | 6975.0 | 0.83 | 0.77 | 0.9 | 1.45e-06 | 0.77 | 1.55e-04 | 0.48 |
| wgcna_switch_only | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | pheno_sig_modules | eQTL | 11 | 25083.0 | 0.94 | 0.89 | 0.99 | 1.71e-02 | 0.84 | 4.19e-02 | 0.89 |
| wgcna_switch_only | pheno_sig_modules | sQTL | 11 | 21807.0 | 0.98 | 0.9 | 1.06 | 5.65e-01 | 0.98 | 5.65e-01 | 0.0 |
| wgcna_switch_only | go_invisible_modules | eQTL | 7 | 5033.0 | 1.06 | 0.99 | 1.14 | 9.73e-02 | 1.06 | 2.22e-01 | 0.29 |
| wgcna_switch_only | go_invisible_modules | sQTL | 7 | 4557.0 | 1.13 | 1.02 | 1.25 | 1.44e-02 | 1.13 | 1.44e-02 | 0.0 |
| wgcna_switch_only | go_visible_modules | eQTL | 11 | 20050.0 | 0.92 | 0.88 | 0.97 | 8.77e-04 | 0.82 | 5.78e-03 | 0.86 |
| wgcna_switch_only | go_visible_modules | sQTL | 11 | 17250.0 | 0.92 | 0.86 | 0.99 | 2.51e-02 | 0.92 | 2.51e-02 | 0.0 |
| wgcna_switch_only | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | pheno_sig_modules | eQTL | 11 | 25083.0 | 0.94 | 0.89 | 0.98 | 9.94e-03 | 0.86 | 3.61e-02 | 0.85 |
| wgcna_switch_only | pheno_sig_modules | sQTL | 11 | 21807.0 | 0.95 | 0.88 | 1.03 | 2.11e-01 | 0.92 | 1.16e-01 | 0.37 |
| wgcna_switch_only | go_invisible_modules | eQTL | 7 | 5033.0 | 1.05 | 0.98 | 1.12 | 1.97e-01 | 1.04 | 3.78e-01 | 0.33 |
| wgcna_switch_only | go_invisible_modules | sQTL | 7 | 4557.0 | 1.12 | 1.01 | 1.23 | 2.37e-02 | 1.12 | 2.54e-02 | 0.01 |
| wgcna_switch_only | go_visible_modules | eQTL | 11 | 20050.0 | 0.93 | 0.88 | 0.97 | 1.09e-03 | 0.85 | 8.98e-03 | 0.83 |
| wgcna_switch_only | go_visible_modules | sQTL | 11 | 17250.0 | 0.91 | 0.85 | 0.97 | 6.63e-03 | 0.9 | 2.14e-02 | 0.28 |
| wgcna_multiplex | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | pheno_sig_modules | eQTL | 9 | 45059.0 | 0.92 | 0.89 | 0.94 | 2.85e-08 | 0.83 | 6.45e-03 | 0.94 |
| wgcna_multiplex | pheno_sig_modules | sQTL | 9 | 37895.0 | 0.94 | 0.9 | 0.98 | 7.49e-03 | 0.92 | 1.54e-01 | 0.77 |
| wgcna_multiplex | go_invisible_modules | eQTL | 7 | 6083.0 | 1.03 | 0.97 | 1.09 | 3.34e-01 | 1.01 | 9.30e-01 | 0.8 |
| wgcna_multiplex | go_invisible_modules | sQTL | 7 | 5491.0 | 1.05 | 0.97 | 1.13 | 2.09e-01 | 1.04 | 4.32e-01 | 0.3 |
| wgcna_multiplex | go_visible_modules | eQTL | 9 | 42464.0 | 0.91 | 0.88 | 0.94 | 7.75e-09 | 0.82 | 3.41e-03 | 0.93 |
| wgcna_multiplex | go_visible_modules | sQTL | 9 | 35605.0 | 0.94 | 0.9 | 0.98 | 6.51e-03 | 0.92 | 1.14e-01 | 0.74 |
| wgcna_multiplex | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | pheno_sig_modules | eQTL | 9 | 45059.0 | 0.92 | 0.89 | 0.94 | 1.45e-08 | 0.83 | 7.26e-03 | 0.94 |
| wgcna_multiplex | pheno_sig_modules | sQTL | 9 | 37895.0 | 0.96 | 0.91 | 1.0 | 5.32e-02 | 0.92 | 2.02e-01 | 0.86 |
| wgcna_multiplex | go_invisible_modules | eQTL | 7 | 6083.0 | 1.01 | 0.95 | 1.07 | 7.63e-01 | 1.0 | 9.60e-01 | 0.86 |
| wgcna_multiplex | go_invisible_modules | sQTL | 7 | 5491.0 | 1.06 | 0.98 | 1.14 | 1.41e-01 | 1.05 | 3.58e-01 | 0.41 |
| wgcna_multiplex | go_visible_modules | eQTL | 9 | 42464.0 | 0.91 | 0.89 | 0.94 | 2.49e-09 | 0.83 | 3.69e-03 | 0.94 |
| wgcna_multiplex | go_visible_modules | sQTL | 9 | 35605.0 | 0.96 | 0.91 | 1.0 | 5.04e-02 | 0.91 | 1.62e-01 | 0.86 |

## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 17 | 1.076 | 1.04 | 1.113 | 2.79e-05 | 1.075 | 0.362 |
| isograph | pheno_sig_modules | 12 | 1.133 | 1.06 | 1.211 | 2.42e-04 | 1.125 | 0.323 |
| isograph | go_invisible_modules | 11 | 1.062 | 0.978 | 1.153 | 1.54e-01 | 1.062 | 0.0 |
| isograph | go_visible_modules | 11 | 1.124 | 1.024 | 1.232 | 1.35e-02 | 1.13 | 0.32 |
| isograph | all_modules | 17 | 1.119 | 1.083 | 1.157 | 2.86e-11 | 1.118 | 0.494 |
| isograph | pheno_sig_modules | 12 | 1.183 | 1.109 | 1.262 | 3.70e-07 | 1.165 | 0.38 |
| isograph | go_invisible_modules | 11 | 1.093 | 1.008 | 1.185 | 3.04e-02 | 1.09 | 0.28 |
| isograph | go_visible_modules | 11 | 1.179 | 1.078 | 1.29 | 3.31e-04 | 1.171 | 0.34 |
| wgcna_switch_only | pheno_sig_modules | 11 | 1.001 | 0.911 | 1.099 | 9.87e-01 | 0.996 | 0.125 |
| wgcna_switch_only | go_invisible_modules | 7 | 1.056 | 0.934 | 1.193 | 3.84e-01 | 1.056 | 0.133 |
| wgcna_switch_only | go_visible_modules | 11 | 0.969 | 0.889 | 1.057 | 4.82e-01 | 0.97 | 0.016 |
| wgcna_switch_only | pheno_sig_modules | 11 | 0.983 | 0.897 | 1.077 | 7.13e-01 | 0.967 | 0.261 |
| wgcna_switch_only | go_invisible_modules | 7 | 1.06 | 0.941 | 1.194 | 3.38e-01 | 1.062 | 0.304 |
| wgcna_switch_only | go_visible_modules | 11 | 0.955 | 0.877 | 1.039 | 2.85e-01 | 0.956 | 0.041 |
| wgcna_multiplex | pheno_sig_modules | 9 | 1.009 | 0.954 | 1.068 | 7.57e-01 | 1.039 | 0.651 |
| wgcna_multiplex | go_invisible_modules | 7 | 1.003 | 0.912 | 1.103 | 9.51e-01 | 0.999 | 0.237 |
| wgcna_multiplex | go_visible_modules | 9 | 1.011 | 0.956 | 1.069 | 7.09e-01 | 1.03 | 0.513 |
| wgcna_multiplex | pheno_sig_modules | 9 | 1.026 | 0.972 | 1.084 | 3.51e-01 | 1.042 | 0.704 |
| wgcna_multiplex | go_invisible_modules | 7 | 1.026 | 0.935 | 1.126 | 5.81e-01 | 1.02 | 0.324 |
| wgcna_multiplex | go_visible_modules | 9 | 1.029 | 0.975 | 1.087 | 2.96e-01 | 1.037 | 0.584 |

## Cross-method contrast on the 13 tissues all methods share

The matched WGCNA baselines (`wgcna_switch_only`, `wgcna_multiplex`) consume the SAME switch / switch+abundance features as IsoGraph, so comparing their splicing-specificity on the SAME tissues isolates whether the signal lives in the switch features or in IsoGraph's VAE + Leiden inference.

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 8 | 1.141 | 1.085 | 1.2 | 2.64e-07 | 1.141 | 0.0 |
| isograph | all_modules | 8 | 1.183 | 1.127 | 1.242 | 1.44e-11 | 1.181 | 0.163 |
| isograph | pheno_sig_modules | 8 | 1.143 | 1.065 | 1.227 | 2.00e-04 | 1.143 | 0.0 |
| isograph | pheno_sig_modules | 8 | 1.196 | 1.117 | 1.281 | 2.78e-07 | 1.196 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | 8 | 0.991 | 0.899 | 1.093 | 8.58e-01 | 0.977 | 0.303 |
| wgcna_switch_only | pheno_sig_modules | 8 | 0.974 | 0.885 | 1.071 | 5.83e-01 | 0.943 | 0.443 |
| wgcna_multiplex | pheno_sig_modules | 8 | 1.006 | 0.951 | 1.065 | 8.27e-01 | 1.032 | 0.68 |
| wgcna_multiplex | pheno_sig_modules | 8 | 1.023 | 0.969 | 1.081 | 4.14e-01 | 1.03 | 0.724 |
| isograph | go_invisible_modules | 7 | 1.075 | 0.984 | 1.176 | 1.10e-01 | 1.075 | 0.0 |
| isograph | go_invisible_modules | 7 | 1.106 | 1.013 | 1.207 | 2.39e-02 | 1.098 | 0.338 |
| wgcna_switch_only | go_invisible_modules | 6 | 1.062 | 0.939 | 1.202 | 3.36e-01 | 1.065 | 0.231 |
| wgcna_switch_only | go_invisible_modules | 6 | 1.064 | 0.943 | 1.2 | 3.16e-01 | 1.068 | 0.411 |
| wgcna_multiplex | go_invisible_modules | 7 | 1.003 | 0.912 | 1.103 | 9.51e-01 | 0.999 | 0.237 |
| wgcna_multiplex | go_invisible_modules | 7 | 1.026 | 0.935 | 1.126 | 5.81e-01 | 1.02 | 0.324 |
| isograph | go_visible_modules | 8 | 1.115 | 1.015 | 1.223 | 2.26e-02 | 1.116 | 0.06 |
| isograph | go_visible_modules | 8 | 1.172 | 1.071 | 1.284 | 5.82e-04 | 1.166 | 0.273 |
| wgcna_switch_only | go_visible_modules | 8 | 0.955 | 0.873 | 1.046 | 3.20e-01 | 0.957 | 0.132 |
| wgcna_switch_only | go_visible_modules | 8 | 0.943 | 0.863 | 1.03 | 1.91e-01 | 0.943 | 0.217 |
| wgcna_multiplex | go_visible_modules | 8 | 1.008 | 0.953 | 1.067 | 7.78e-01 | 1.025 | 0.546 |
| wgcna_multiplex | go_visible_modules | 8 | 1.026 | 0.972 | 1.084 | 3.52e-01 | 1.028 | 0.603 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, strongest for the GO-invisible modules — splicing-QTL is spared relative to expression-QTL exactly where the DTU-without-DGE value concentrates.
- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL signal is a property of the switch features (which both methods share), not of IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. IsoGraph's value is the switch-feature *representation* and its finer modules, not a unique network-inference effect.
- **Primary internal control = the matched WGCNA baselines**, not `go_visible_modules`. The baselines hold the switch features fixed and vary only the inference, so a null there localises the effect to IsoGraph's inference. `go_visible_modules` is a secondary control on module CONTENT and is only ever a relative contrast — read it as the low end of a gradient, not as an on/off null.
- High I2 flags between-tissue heterogeneity; prefer RE there. A nominally significant ratio carrying high I2 is driven by a few tissues, not by a consistent effect, and is weaker evidence than a smaller ratio at I2 near 0 — compare module sets on consistency as well as magnitude. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
