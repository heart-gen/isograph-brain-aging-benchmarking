# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph, wgcna_switch_only, wgcna_multiplex.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 80925.0 | 0.82 | 0.8 | 0.83 | 5.16e-122 | 0.82 | 2.74e-74 | 0.39 |
| isograph | all_modules | sQTL | 17 | 65578.0 | 0.88 | 0.86 | 0.9 | 1.63e-25 | 0.88 | 1.66e-12 | 0.55 |
| isograph | pheno_sig_modules | eQTL | 12 | 16598.0 | 0.76 | 0.73 | 0.79 | 1.21e-55 | 0.75 | 1.49e-18 | 0.64 |
| isograph | pheno_sig_modules | sQTL | 12 | 12867.0 | 0.86 | 0.82 | 0.91 | 2.18e-08 | 0.86 | 6.22e-03 | 0.69 |
| isograph | go_invisible_modules | eQTL | 10 | 8646.0 | 0.88 | 0.84 | 0.92 | 1.12e-08 | 0.89 | 9.12e-04 | 0.54 |
| isograph | go_invisible_modules | sQTL | 10 | 7197.0 | 0.99 | 0.93 | 1.05 | 7.92e-01 | 0.99 | 7.98e-01 | 0.41 |
| isograph | go_visible_modules | eQTL | 11 | 7952.0 | 0.67 | 0.64 | 0.71 | 1.22e-56 | 0.6 | 1.52e-12 | 0.84 |
| isograph | go_visible_modules | sQTL | 11 | 5670.0 | 0.72 | 0.66 | 0.78 | 1.06e-15 | 0.68 | 1.73e-07 | 0.51 |
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
| isograph | all_modules | 17 | 1.071 | 1.039 | 1.104 | 7.78e-06 | 1.071 | 0.349 |
| isograph | pheno_sig_modules | 12 | 1.127 | 1.059 | 1.198 | 1.50e-04 | 1.131 | 0.222 |
| isograph | go_invisible_modules | 10 | 1.125 | 1.042 | 1.215 | 2.61e-03 | 1.12 | 0.151 |
| isograph | go_visible_modules | 11 | 1.037 | 0.942 | 1.14 | 4.62e-01 | 1.035 | 0.486 |
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
| isograph | all_modules | 8 | 1.107 | 1.059 | 1.156 | 6.49e-06 | 1.106 | 0.518 |
| isograph | pheno_sig_modules | 8 | 1.108 | 1.038 | 1.184 | 2.12e-03 | 1.105 | 0.32 |
| wgcna_switch_only | pheno_sig_modules | 8 | 0.968 | 0.889 | 1.055 | 4.63e-01 | 0.951 | 0.25 |
| wgcna_multiplex | pheno_sig_modules | 8 | 0.98 | 0.934 | 1.029 | 4.15e-01 | 0.983 | 0.741 |
| isograph | go_invisible_modules | 7 | 1.114 | 1.025 | 1.211 | 1.13e-02 | 1.095 | 0.412 |
| wgcna_switch_only | go_invisible_modules | 6 | 1.021 | 0.914 | 1.14 | 7.12e-01 | 1.036 | 0.456 |
| wgcna_multiplex | go_invisible_modules | 7 | 0.993 | 0.913 | 1.079 | 8.63e-01 | 0.993 | 0.0 |
| isograph | go_visible_modules | 8 | 1.02 | 0.925 | 1.124 | 6.90e-01 | 0.996 | 0.493 |
| wgcna_switch_only | go_visible_modules | 8 | 0.959 | 0.886 | 1.038 | 3.03e-01 | 0.959 | 0.087 |
| wgcna_multiplex | go_visible_modules | 8 | 0.993 | 0.946 | 1.042 | 7.72e-01 | 0.995 | 0.674 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, strongest for the GO-invisible modules — splicing-QTL is spared relative to expression-QTL exactly where the DTU-without-DGE value concentrates.
- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL signal is a property of the switch features (which both methods share), not of IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. IsoGraph's value is the switch-feature *representation* and its finer modules, not a unique network-inference effect.
- High I2 flags between-tissue heterogeneity; prefer RE there. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
