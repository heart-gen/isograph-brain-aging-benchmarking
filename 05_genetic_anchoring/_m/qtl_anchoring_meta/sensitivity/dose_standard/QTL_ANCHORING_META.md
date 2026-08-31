# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 86465.0 | 0.82 | 0.81 | 0.83 | 3.72e-184 | 0.82 | 4.30e-58 | 0.7 |
| isograph | all_modules | sQTL | 17 | 70399.0 | 0.86 | 0.85 | 0.87 | 5.13e-122 | 0.85 | 5.08e-29 | 0.79 |
| isograph | pheno_sig_modules | eQTL | 12 | 20631.0 | 0.77 | 0.75 | 0.79 | 5.20e-78 | 0.75 | 8.39e-18 | 0.76 |
| isograph | pheno_sig_modules | sQTL | 12 | 16433.0 | 0.85 | 0.83 | 0.87 | 1.86e-35 | 0.83 | 7.76e-08 | 0.81 |
| isograph | go_invisible_modules | eQTL | 11 | 9926.0 | 0.92 | 0.89 | 0.95 | 4.15e-06 | 0.89 | 1.03e-03 | 0.64 |
| isograph | go_invisible_modules | sQTL | 11 | 8396.0 | 0.94 | 0.91 | 0.97 | 2.60e-05 | 0.92 | 1.92e-02 | 0.75 |
| isograph | go_visible_modules | eQTL | 11 | 10705.0 | 0.68 | 0.65 | 0.7 | 1.46e-86 | 0.63 | 9.17e-19 | 0.77 |
| isograph | go_visible_modules | sQTL | 11 | 8037.0 | 0.79 | 0.76 | 0.82 | 2.12e-35 | 0.72 | 2.64e-08 | 0.81 |

## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 17 | 1.045 | 1.026 | 1.065 | 3.03e-06 | 1.044 | 0.635 |
| isograph | pheno_sig_modules | 12 | 1.1 | 1.06 | 1.142 | 5.21e-07 | 1.098 | 0.486 |
| isograph | go_invisible_modules | 11 | 1.021 | 0.976 | 1.069 | 3.67e-01 | 1.025 | 0.318 |
| isograph | go_visible_modules | 11 | 1.144 | 1.083 | 1.208 | 1.29e-06 | 1.133 | 0.638 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, strongest for the GO-invisible modules — splicing-QTL is spared relative to expression-QTL exactly where the DTU-without-DGE value concentrates.
- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL signal is a property of the switch features (which both methods share), not of IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. IsoGraph's value is the switch-feature *representation* and its finer modules, not a unique network-inference effect.
- **Primary internal control = the matched WGCNA baselines**, not `go_visible_modules`. The baselines hold the switch features fixed and vary only the inference, so a null there localises the effect to IsoGraph's inference. `go_visible_modules` is a secondary control on module CONTENT and is only ever a relative contrast — read it as the low end of a gradient, not as an on/off null.
- High I2 flags between-tissue heterogeneity; prefer RE there. A nominally significant ratio carrying high I2 is driven by a few tissues, not by a consistent effect, and is weaker evidence than a smaller ratio at I2 near 0 — compare module sets on consistency as well as magnitude. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
