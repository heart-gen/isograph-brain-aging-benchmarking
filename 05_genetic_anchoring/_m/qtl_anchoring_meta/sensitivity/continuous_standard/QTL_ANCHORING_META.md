# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph, wgcna_switch_only, wgcna_multiplex.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 120520.0 | 0.9 | 0.89 | 0.91 | 2.09e-166 | 0.9 | 1.72e-27 | 0.84 |
| isograph | all_modules | sQTL | 17 | 104145.0 | 0.95 | 0.94 | 0.96 | 8.24e-35 | 0.95 | 7.05e-20 | 0.45 |
| isograph | pheno_sig_modules | eQTL | 10 | 21586.0 | 0.91 | 0.89 | 0.92 | 1.69e-37 | 0.9 | 4.22e-07 | 0.83 |
| isograph | pheno_sig_modules | sQTL | 10 | 18252.0 | 0.95 | 0.94 | 0.97 | 5.59e-10 | 0.94 | 1.19e-02 | 0.85 |
| isograph | go_invisible_modules | eQTL | 7 | 5455.0 | 0.98 | 0.95 | 1.01 | 1.28e-01 | 0.97 | 2.17e-01 | 0.57 |
| isograph | go_invisible_modules | sQTL | 7 | 4971.0 | 1.01 | 0.98 | 1.04 | 6.33e-01 | 0.99 | 7.85e-01 | 0.7 |
| isograph | go_visible_modules | eQTL | 9 | 16131.0 | 0.89 | 0.87 | 0.9 | 1.27e-42 | 0.88 | 6.53e-14 | 0.63 |
| isograph | go_visible_modules | sQTL | 9 | 13281.0 | 0.93 | 0.91 | 0.95 | 7.17e-14 | 0.94 | 2.54e-03 | 0.73 |
| wgcna_switch_only | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | pheno_sig_modules | eQTL | 6 | 22366.0 | 1.01 | 0.99 | 1.03 | 4.17e-01 | 1.01 | 5.24e-01 | 0.39 |
| wgcna_switch_only | pheno_sig_modules | sQTL | 6 | 20534.0 | 0.98 | 0.96 | 1.01 | 1.68e-01 | 0.98 | 3.82e-01 | 0.81 |
| wgcna_switch_only | go_invisible_modules | eQTL | 6 | 20788.0 | 1.01 | 0.98 | 1.03 | 5.71e-01 | 1.01 | 6.51e-01 | 0.37 |
| wgcna_switch_only | go_invisible_modules | sQTL | 6 | 19102.0 | 0.98 | 0.96 | 1.0 | 9.31e-02 | 0.97 | 3.11e-01 | 0.78 |
| wgcna_switch_only | go_visible_modules | eQTL | 2 | 1578.0 | 1.02 | 0.96 | 1.07 | 5.49e-01 | 1.02 | 5.49e-01 | 0.0 |
| wgcna_switch_only | go_visible_modules | sQTL | 2 | 1432.0 | 1.02 | 0.96 | 1.08 | 5.27e-01 | 1.02 | 5.20e-01 | 0.14 |
| wgcna_multiplex | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | pheno_sig_modules | eQTL | 8 | 35511.0 | 0.94 | 0.93 | 0.96 | 2.96e-13 | 0.9 | 6.83e-03 | 0.96 |
| wgcna_multiplex | pheno_sig_modules | sQTL | 8 | 28196.0 | 0.98 | 0.97 | 1.0 | 7.84e-02 | 0.98 | 2.37e-01 | 0.76 |
| wgcna_multiplex | go_invisible_modules | eQTL | 4 | 5804.0 | 0.92 | 0.89 | 0.95 | 6.41e-09 | 0.89 | 1.81e-03 | 0.75 |
| wgcna_multiplex | go_invisible_modules | sQTL | 4 | 5162.0 | 0.97 | 0.94 | 1.0 | 5.78e-02 | 0.97 | 3.29e-01 | 0.61 |
| wgcna_multiplex | go_visible_modules | eQTL | 8 | 32960.0 | 0.96 | 0.94 | 0.97 | 1.01e-08 | 0.91 | 1.90e-02 | 0.96 |
| wgcna_multiplex | go_visible_modules | sQTL | 8 | 26053.0 | 0.99 | 0.98 | 1.01 | 4.93e-01 | 0.98 | 4.34e-01 | 0.86 |

## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 17 | 1.054 | 1.042 | 1.065 | 1.93e-20 | 1.054 | 0.687 |
| isograph | pheno_sig_modules | 10 | 1.045 | 1.023 | 1.069 | 7.85e-05 | 1.045 | 0.0 |
| isograph | go_invisible_modules | 7 | 1.028 | 0.989 | 1.069 | 1.62e-01 | 1.028 | 0.0 |
| isograph | go_visible_modules | 9 | 1.048 | 1.022 | 1.075 | 2.57e-04 | 1.052 | 0.252 |
| wgcna_switch_only | pheno_sig_modules | 6 | 0.974 | 0.942 | 1.007 | 1.19e-01 | 0.969 | 0.51 |
| wgcna_switch_only | go_invisible_modules | 6 | 0.973 | 0.94 | 1.006 | 1.08e-01 | 0.968 | 0.501 |
| wgcna_switch_only | go_visible_modules | 2 | 1.001 | 0.926 | 1.083 | 9.72e-01 | 1.001 | 0.0 |
| wgcna_multiplex | pheno_sig_modules | 8 | 1.04 | 1.016 | 1.064 | 8.36e-04 | 1.06 | 0.826 |
| wgcna_multiplex | go_invisible_modules | 4 | 1.055 | 1.012 | 1.099 | 1.08e-02 | 1.056 | 0.062 |
| wgcna_multiplex | go_visible_modules | 8 | 1.036 | 1.013 | 1.06 | 2.20e-03 | 1.057 | 0.818 |

## Cross-method contrast on the 13 tissues all methods share

The matched WGCNA baselines (`wgcna_switch_only`, `wgcna_multiplex`) consume the SAME switch / switch+abundance features as IsoGraph, so comparing their splicing-specificity on the SAME tissues isolates whether the signal lives in the switch features or in IsoGraph's VAE + Leiden inference.

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 5 | 1.04 | 1.019 | 1.061 | 1.37e-04 | 1.04 | 0.382 |
| isograph | pheno_sig_modules | 5 | 1.039 | 1.013 | 1.065 | 2.83e-03 | 1.039 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | 5 | 0.964 | 0.929 | 1.001 | 5.47e-02 | 0.959 | 0.546 |
| wgcna_multiplex | pheno_sig_modules | 5 | 1.03 | 1.005 | 1.056 | 1.68e-02 | 1.032 | 0.864 |
| isograph | go_invisible_modules | 4 | 1.026 | 0.972 | 1.083 | 3.57e-01 | 1.026 | 0.0 |
| wgcna_switch_only | go_invisible_modules | 5 | 0.962 | 0.927 | 0.999 | 4.35e-02 | 0.958 | 0.516 |
| wgcna_multiplex | go_invisible_modules | 3 | 1.051 | 1.009 | 1.095 | 1.76e-02 | 1.051 | 0.0 |
| isograph | go_visible_modules | 5 | 1.039 | 1.011 | 1.067 | 5.66e-03 | 1.039 | 0.0 |
| wgcna_switch_only | go_visible_modules | 1 | 1.032 | 0.874 | 1.218 | 7.10e-01 | 1.032 | 0.0 |
| wgcna_multiplex | go_visible_modules | 5 | 1.028 | 1.003 | 1.053 | 2.63e-02 | 1.032 | 0.857 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, carried by the phenotype-associated modules — splicing-QTL is spared relative to expression-QTL in the modules that track the trait.
- **Do NOT read a GO-invisible localisation off this table.** On the 2026-08-29 refresh the `go_invisible_modules` and `go_visible_modules` arms are statistically indistinguishable, so neither is the site of the effect and neither is a null. The DTU-without-DGE claim rests on the GO-invisible content gate, not on these genetics. An earlier framing made `go_invisible` the headline; it came from a stale `module_enrichment` join that relabelled the partition and must not be restored.
- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL signal is a property of the switch features (which both methods share), not of IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. IsoGraph's value is the switch-feature *representation* and its finer modules, not a unique network-inference effect.
- **Primary internal control = the matched WGCNA baselines**, not `go_visible_modules`. The baselines hold the switch features fixed and vary only the inference, so a null there localises the effect to IsoGraph's inference. `go_visible_modules` is a secondary, exploratory stratification on module CONTENT and is only ever a relative contrast — never an on/off null.
- High I2 flags between-tissue heterogeneity; prefer RE there. A nominally significant ratio carrying high I2 is driven by a few tissues, not by a consistent effect, and is weaker evidence than a smaller ratio at I2 near 0 — compare module sets on consistency as well as magnitude. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
