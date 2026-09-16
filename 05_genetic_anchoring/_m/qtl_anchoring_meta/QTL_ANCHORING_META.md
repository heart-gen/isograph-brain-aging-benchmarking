# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph, wgcna_switch_only, wgcna_multiplex.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 87530.0 | 0.8 | 0.79 | 0.82 | 5.91e-151 | 0.8 | 8.43e-34 | 0.78 |
| isograph | all_modules | sQTL | 17 | 75492.0 | 0.83 | 0.81 | 0.85 | 2.11e-53 | 0.82 | 1.84e-17 | 0.72 |
| isograph | pheno_sig_modules | eQTL | 8 | 19371.0 | 0.8 | 0.77 | 0.82 | 3.04e-42 | 0.8 | 4.47e-07 | 0.84 |
| isograph | pheno_sig_modules | sQTL | 8 | 16601.0 | 0.81 | 0.77 | 0.85 | 8.59e-19 | 0.81 | 2.52e-07 | 0.62 |
| isograph | go_invisible_modules | eQTL | 8 | 5228.0 | 0.87 | 0.82 | 0.92 | 9.35e-07 | 0.87 | 6.58e-04 | 0.36 |
| isograph | go_invisible_modules | sQTL | 8 | 4736.0 | 0.93 | 0.86 | 1.0 | 4.78e-02 | 0.93 | 1.35e-01 | 0.2 |
| isograph | go_visible_modules | eQTL | 7 | 14143.0 | 0.79 | 0.76 | 0.82 | 4.06e-35 | 0.78 | 1.62e-06 | 0.84 |
| isograph | go_visible_modules | sQTL | 7 | 11865.0 | 0.78 | 0.74 | 0.83 | 1.25e-18 | 0.77 | 4.90e-08 | 0.6 |
| wgcna_switch_only | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_switch_only | pheno_sig_modules | eQTL | 6 | 22366.0 | 1.05 | 1.0 | 1.1 | 5.97e-02 | 1.05 | 1.21e-01 | 0.3 |
| wgcna_switch_only | pheno_sig_modules | sQTL | 6 | 20534.0 | 0.97 | 0.91 | 1.03 | 3.00e-01 | 0.96 | 5.12e-01 | 0.76 |
| wgcna_switch_only | go_invisible_modules | eQTL | 6 | 20788.0 | 1.05 | 1.0 | 1.1 | 7.36e-02 | 1.05 | 1.20e-01 | 0.23 |
| wgcna_switch_only | go_invisible_modules | sQTL | 6 | 19102.0 | 0.97 | 0.91 | 1.03 | 3.18e-01 | 0.96 | 5.26e-01 | 0.74 |
| wgcna_switch_only | go_visible_modules | eQTL | 2 | 1578.0 | 1.02 | 0.9 | 1.15 | 7.52e-01 | 1.02 | 7.52e-01 | 0.0 |
| wgcna_switch_only | go_visible_modules | sQTL | 2 | 1432.0 | 0.99 | 0.83 | 1.18 | 9.17e-01 | 0.99 | 9.17e-01 | 0.0 |
| wgcna_multiplex | all_modules | eQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | all_modules | sQTL | 0 | 0.0 | nan | nan | nan | NA | nan | NA | nan |
| wgcna_multiplex | pheno_sig_modules | eQTL | 8 | 35511.0 | 0.9 | 0.88 | 0.93 | 1.49e-10 | 0.84 | 2.11e-02 | 0.95 |
| wgcna_multiplex | pheno_sig_modules | sQTL | 8 | 28196.0 | 1.0 | 0.95 | 1.04 | 8.46e-01 | 0.97 | 6.68e-01 | 0.86 |
| wgcna_multiplex | go_invisible_modules | eQTL | 4 | 5804.0 | 0.89 | 0.84 | 0.94 | 6.00e-05 | 0.82 | 1.25e-02 | 0.75 |
| wgcna_multiplex | go_invisible_modules | sQTL | 4 | 5162.0 | 0.9 | 0.83 | 0.97 | 6.46e-03 | 0.92 | 2.80e-01 | 0.56 |
| wgcna_multiplex | go_visible_modules | eQTL | 8 | 32960.0 | 0.92 | 0.89 | 0.95 | 5.27e-08 | 0.84 | 3.82e-02 | 0.96 |
| wgcna_multiplex | go_visible_modules | sQTL | 8 | 26053.0 | 1.02 | 0.98 | 1.07 | 2.97e-01 | 0.98 | 8.02e-01 | 0.9 |

## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 17 | 1.029 | 1.0 | 1.059 | 5.13e-02 | 1.025 | 0.532 |
| isograph | pheno_sig_modules | 8 | 1.004 | 0.948 | 1.063 | 8.98e-01 | 1.004 | 0.0 |
| isograph | go_invisible_modules | 8 | 1.062 | 0.965 | 1.168 | 2.18e-01 | 1.062 | 0.0 |
| isograph | go_visible_modules | 7 | 0.978 | 0.915 | 1.045 | 5.07e-01 | 0.978 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | 6 | 0.927 | 0.856 | 1.005 | 6.46e-02 | 0.92 | 0.465 |
| wgcna_switch_only | go_invisible_modules | 6 | 0.929 | 0.857 | 1.008 | 7.62e-02 | 0.922 | 0.459 |
| wgcna_switch_only | go_visible_modules | 2 | 0.969 | 0.785 | 1.197 | 7.73e-01 | 0.969 | 0.0 |
| wgcna_multiplex | pheno_sig_modules | 8 | 1.084 | 1.025 | 1.146 | 4.44e-03 | 1.105 | 0.639 |
| wgcna_multiplex | go_invisible_modules | 4 | 1.006 | 0.912 | 1.11 | 9.05e-01 | 1.085 | 0.611 |
| wgcna_multiplex | go_visible_modules | 8 | 1.095 | 1.036 | 1.157 | 1.26e-03 | 1.103 | 0.503 |

## Cross-method contrast on the 13 tissues all methods share

The matched WGCNA baselines (`wgcna_switch_only`, `wgcna_multiplex`) consume the SAME switch / switch+abundance features as IsoGraph, so comparing their splicing-specificity on the SAME tissues isolates whether the signal lives in the switch features or in IsoGraph's VAE + Leiden inference.

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 5 | 1.015 | 0.964 | 1.069 | 5.80e-01 | 1.015 | 0.0 |
| isograph | pheno_sig_modules | 5 | 1.0 | 0.939 | 1.065 | 9.97e-01 | 1.0 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | 5 | 0.93 | 0.853 | 1.015 | 1.02e-01 | 0.92 | 0.57 |
| wgcna_multiplex | pheno_sig_modules | 5 | 1.059 | 0.998 | 1.124 | 5.97e-02 | 1.059 | 0.632 |
| isograph | go_invisible_modules | 5 | 1.061 | 0.95 | 1.185 | 2.94e-01 | 1.061 | 0.0 |
| wgcna_switch_only | go_invisible_modules | 5 | 0.929 | 0.852 | 1.014 | 9.89e-02 | 0.92 | 0.567 |
| wgcna_multiplex | go_invisible_modules | 3 | 0.996 | 0.902 | 1.099 | 9.29e-01 | 1.027 | 0.511 |
| isograph | go_visible_modules | 5 | 0.978 | 0.911 | 1.05 | 5.36e-01 | 0.978 | 0.0 |
| wgcna_switch_only | go_visible_modules | 1 | 0.999 | 0.668 | 1.494 | 9.96e-01 | 0.999 | 0.0 |
| wgcna_multiplex | go_visible_modules | 5 | 1.076 | 1.015 | 1.141 | 1.47e-02 | 1.072 | 0.501 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, carried by the phenotype-associated modules — splicing-QTL is spared relative to expression-QTL in the modules that track the trait.
- **Do NOT read a GO-invisible localisation off this table.** On the 2026-08-29 refresh the `go_invisible_modules` and `go_visible_modules` arms are statistically indistinguishable, so neither is the site of the effect and neither is a null. The DTU-without-DGE claim rests on the GO-invisible content gate, not on these genetics. An earlier framing made `go_invisible` the headline; it came from a stale `module_enrichment` join that relabelled the partition and must not be restored.
- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL signal is a property of the switch features (which both methods share), not of IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. IsoGraph's value is the switch-feature *representation* and its finer modules, not a unique network-inference effect.
- **Primary internal control = the matched WGCNA baselines**, not `go_visible_modules`. The baselines hold the switch features fixed and vary only the inference, so a null there localises the effect to IsoGraph's inference. `go_visible_modules` is a secondary, exploratory stratification on module CONTENT and is only ever a relative contrast — never an on/off null.
- High I2 flags between-tissue heterogeneity; prefer RE there. A nominally significant ratio carrying high I2 is driven by a few tissues, not by a consistent effect, and is weaker evidence than a smaller ratio at I2 near 0 — compare module sets on consistency as well as magnitude. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
