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
| wgcna_switch_only | pheno_sig_modules | eQTL | 9 | 30231.0 | 0.92 | 0.87 | 0.96 | 2.06e-04 | 0.85 | 1.62e-02 | 0.86 |
| wgcna_switch_only | pheno_sig_modules | sQTL | 9 | 25255.0 | 0.96 | 0.89 | 1.03 | 2.82e-01 | 0.91 | 1.45e-01 | 0.55 |
| wgcna_switch_only | go_invisible_modules | eQTL | 8 | 9898.0 | 1.02 | 0.97 | 1.08 | 3.65e-01 | 1.0 | 9.62e-01 | 0.8 |
| wgcna_switch_only | go_invisible_modules | sQTL | 8 | 8873.0 | 1.1 | 1.02 | 1.18 | 9.09e-03 | 1.04 | 6.31e-01 | 0.71 |
| wgcna_switch_only | go_visible_modules | eQTL | 10 | 20465.0 | 0.92 | 0.88 | 0.96 | 5.61e-05 | 0.84 | 1.59e-03 | 0.82 |
| wgcna_switch_only | go_visible_modules | sQTL | 10 | 16494.0 | 0.9 | 0.85 | 0.96 | 1.60e-03 | 0.9 | 1.41e-02 | 0.32 |
| wgcna_multiplex | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | pheno_sig_modules | eQTL | 10 | 57388.0 | 0.92 | 0.9 | 0.94 | 1.21e-10 | 0.82 | 3.13e-03 | 0.96 |
| wgcna_multiplex | pheno_sig_modules | sQTL | 10 | 44300.0 | 0.98 | 0.95 | 1.02 | 4.21e-01 | 0.93 | 3.07e-01 | 0.9 |
| wgcna_multiplex | go_invisible_modules | eQTL | 6 | 5428.0 | 0.97 | 0.91 | 1.02 | 2.41e-01 | 0.93 | 5.50e-01 | 0.92 |
| wgcna_multiplex | go_invisible_modules | sQTL | 6 | 4590.0 | 1.01 | 0.94 | 1.09 | 7.69e-01 | 0.98 | 8.56e-01 | 0.71 |
| wgcna_multiplex | go_visible_modules | eQTL | 10 | 54638.0 | 0.92 | 0.9 | 0.95 | 5.69e-10 | 0.82 | 2.11e-03 | 0.95 |
| wgcna_multiplex | go_visible_modules | sQTL | 10 | 42083.0 | 0.99 | 0.95 | 1.02 | 4.42e-01 | 0.93 | 2.82e-01 | 0.89 |

## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 17 | 1.068 | 1.037 | 1.1 | 1.35e-05 | 1.066 | 0.291 |
| isograph | pheno_sig_modules | 12 | 1.111 | 1.048 | 1.177 | 3.57e-04 | 1.107 | 0.229 |
| isograph | go_invisible_modules | 11 | 1.068 | 0.993 | 1.148 | 7.72e-02 | 1.068 | 0.0 |
| isograph | go_visible_modules | 11 | 1.084 | 1.0 | 1.174 | 5.02e-02 | 1.082 | 0.284 |
| wgcna_switch_only | pheno_sig_modules | 9 | 1.025 | 0.94 | 1.118 | 5.77e-01 | 1.013 | 0.228 |
| wgcna_switch_only | go_invisible_modules | 8 | 1.066 | 0.977 | 1.164 | 1.49e-01 | 1.054 | 0.355 |
| wgcna_switch_only | go_visible_modules | 10 | 0.968 | 0.897 | 1.045 | 4.08e-01 | 0.968 | 0.0 |
| wgcna_multiplex | pheno_sig_modules | 10 | 1.046 | 0.999 | 1.095 | 5.31e-02 | 1.079 | 0.748 |
| wgcna_multiplex | go_invisible_modules | 6 | 1.019 | 0.924 | 1.123 | 7.08e-01 | 1.019 | 0.0 |
| wgcna_multiplex | go_visible_modules | 10 | 1.044 | 0.997 | 1.093 | 6.50e-02 | 1.074 | 0.751 |

## Cross-method contrast on the 13 tissues all methods share

The matched WGCNA baselines (`wgcna_switch_only`, `wgcna_multiplex`) consume the SAME switch / switch+abundance features as IsoGraph, so comparing their splicing-specificity on the SAME tissues isolates whether the signal lives in the switch features or in IsoGraph's VAE + Leiden inference.

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 8 | 1.11 | 1.063 | 1.159 | 2.37e-06 | 1.11 | 0.0 |
| isograph | pheno_sig_modules | 7 | 1.112 | 1.046 | 1.182 | 6.53e-04 | 1.112 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | 7 | 1.015 | 0.927 | 1.113 | 7.41e-01 | 0.994 | 0.39 |
| wgcna_multiplex | pheno_sig_modules | 8 | 1.044 | 0.996 | 1.093 | 7.02e-02 | 1.071 | 0.798 |
| isograph | go_invisible_modules | 7 | 1.066 | 0.985 | 1.153 | 1.12e-01 | 1.061 | 0.157 |
| wgcna_switch_only | go_invisible_modules | 7 | 1.065 | 0.975 | 1.162 | 1.61e-01 | 1.049 | 0.441 |
| wgcna_multiplex | go_invisible_modules | 6 | 1.019 | 0.924 | 1.123 | 7.08e-01 | 1.019 | 0.0 |
| isograph | go_visible_modules | 7 | 1.083 | 0.998 | 1.175 | 5.51e-02 | 1.083 | 0.0 |
| wgcna_switch_only | go_visible_modules | 8 | 0.957 | 0.883 | 1.036 | 2.77e-01 | 0.957 | 0.0 |
| wgcna_multiplex | go_visible_modules | 8 | 1.041 | 0.994 | 1.091 | 8.52e-02 | 1.065 | 0.8 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, carried by the phenotype-associated modules — splicing-QTL is spared relative to expression-QTL in the modules that track the trait.
- **Do NOT read a GO-invisible localisation off this table.** On the 2026-08-29 refresh the `go_invisible_modules` and `go_visible_modules` arms are statistically indistinguishable, so neither is the site of the effect and neither is a null. The DTU-without-DGE claim rests on the GO-invisible content gate, not on these genetics. An earlier framing made `go_invisible` the headline; it came from a stale `module_enrichment` join that relabelled the partition and must not be restored.
- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL signal is a property of the switch features (which both methods share), not of IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. IsoGraph's value is the switch-feature *representation* and its finer modules, not a unique network-inference effect.
- **Primary internal control = the matched WGCNA baselines**, not `go_visible_modules`. The baselines hold the switch features fixed and vary only the inference, so a null there localises the effect to IsoGraph's inference. `go_visible_modules` is a secondary, exploratory stratification on module CONTENT and is only ever a relative contrast — never an on/off null.
- High I2 flags between-tissue heterogeneity; prefer RE there. A nominally significant ratio carrying high I2 is driven by a few tissues, not by a consistent effect, and is weaker evidence than a smaller ratio at I2 near 0 — compare module sets on consistency as well as magnitude. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
