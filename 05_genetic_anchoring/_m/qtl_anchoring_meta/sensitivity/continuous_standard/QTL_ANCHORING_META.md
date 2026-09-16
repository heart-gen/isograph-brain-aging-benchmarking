# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph, wgcna_switch_only, wgcna_multiplex.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 87530.0 | 0.88 | 0.87 | 0.89 | 1.56e-212 | 0.88 | 6.13e-39 | 0.82 |
| isograph | all_modules | sQTL | 17 | 75492.0 | 0.93 | 0.92 | 0.94 | 7.61e-57 | 0.93 | 2.06e-44 | 0.22 |
| isograph | pheno_sig_modules | eQTL | 8 | 19371.0 | 0.89 | 0.87 | 0.9 | 1.71e-51 | 0.89 | 5.09e-12 | 0.74 |
| isograph | pheno_sig_modules | sQTL | 8 | 16601.0 | 0.92 | 0.9 | 0.93 | 3.13e-25 | 0.92 | 9.51e-25 | 0.02 |
| isograph | go_invisible_modules | eQTL | 8 | 5228.0 | 0.94 | 0.91 | 0.96 | 5.35e-06 | 0.94 | 2.16e-03 | 0.36 |
| isograph | go_invisible_modules | sQTL | 8 | 4736.0 | 0.96 | 0.93 | 0.99 | 3.40e-03 | 0.96 | 3.40e-03 | 0.0 |
| isograph | go_visible_modules | eQTL | 7 | 14143.0 | 0.88 | 0.86 | 0.89 | 2.67e-46 | 0.88 | 4.92e-14 | 0.68 |
| isograph | go_visible_modules | sQTL | 7 | 11865.0 | 0.91 | 0.89 | 0.93 | 3.53e-23 | 0.91 | 3.53e-23 | 0.0 |
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
| isograph | all_modules | 17 | 1.057 | 1.045 | 1.069 | 1.76e-20 | 1.056 | 0.62 |
| isograph | pheno_sig_modules | 8 | 1.031 | 1.008 | 1.055 | 9.01e-03 | 1.031 | 0.0 |
| isograph | go_invisible_modules | 8 | 1.022 | 0.981 | 1.063 | 2.98e-01 | 1.022 | 0.0 |
| isograph | go_visible_modules | 7 | 1.033 | 1.006 | 1.06 | 1.61e-02 | 1.033 | 0.0 |
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
| isograph | all_modules | 5 | 1.036 | 1.014 | 1.059 | 1.15e-03 | 1.036 | 0.0 |
| isograph | pheno_sig_modules | 5 | 1.027 | 1.0 | 1.053 | 4.71e-02 | 1.027 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | 5 | 0.964 | 0.929 | 1.001 | 5.47e-02 | 0.959 | 0.546 |
| wgcna_multiplex | pheno_sig_modules | 5 | 1.03 | 1.005 | 1.056 | 1.68e-02 | 1.032 | 0.864 |
| isograph | go_invisible_modules | 5 | 1.026 | 0.979 | 1.077 | 2.85e-01 | 1.026 | 0.0 |
| wgcna_switch_only | go_invisible_modules | 5 | 0.962 | 0.927 | 0.999 | 4.35e-02 | 0.958 | 0.516 |
| wgcna_multiplex | go_invisible_modules | 3 | 1.051 | 1.009 | 1.095 | 1.76e-02 | 1.051 | 0.0 |
| isograph | go_visible_modules | 5 | 1.025 | 0.996 | 1.055 | 9.74e-02 | 1.025 | 0.0 |
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
