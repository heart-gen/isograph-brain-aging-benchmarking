# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph, wgcna_switch_only, wgcna_multiplex.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 120520.0 | 0.84 | 0.83 | 0.85 | 3.32e-113 | 0.84 | 7.19e-24 | 0.8 |
| isograph | all_modules | sQTL | 17 | 104145.0 | 0.87 | 0.85 | 0.89 | 7.28e-35 | 0.86 | 3.66e-12 | 0.71 |
| isograph | pheno_sig_modules | eQTL | 10 | 21586.0 | 0.85 | 0.82 | 0.87 | 1.63e-26 | 0.83 | 7.59e-06 | 0.83 |
| isograph | pheno_sig_modules | sQTL | 10 | 18252.0 | 0.92 | 0.88 | 0.96 | 8.85e-05 | 0.88 | 4.27e-02 | 0.84 |
| isograph | go_invisible_modules | eQTL | 7 | 5455.0 | 0.95 | 0.9 | 1.01 | 8.07e-02 | 0.95 | 2.35e-01 | 0.55 |
| isograph | go_invisible_modules | sQTL | 7 | 4971.0 | 1.04 | 0.97 | 1.11 | 2.96e-01 | 1.01 | 8.31e-01 | 0.61 |
| isograph | go_visible_modules | eQTL | 9 | 16131.0 | 0.82 | 0.79 | 0.85 | 1.50e-28 | 0.8 | 8.41e-07 | 0.79 |
| isograph | go_visible_modules | sQTL | 9 | 13281.0 | 0.87 | 0.83 | 0.92 | 1.18e-07 | 0.88 | 5.29e-02 | 0.79 |
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
| isograph | all_modules | 17 | 1.033 | 1.006 | 1.062 | 1.79e-02 | 1.03 | 0.604 |
| isograph | pheno_sig_modules | 10 | 1.065 | 1.01 | 1.123 | 1.98e-02 | 1.066 | 0.104 |
| isograph | go_invisible_modules | 7 | 1.08 | 0.987 | 1.182 | 9.27e-02 | 1.08 | 0.0 |
| isograph | go_visible_modules | 9 | 1.049 | 0.987 | 1.116 | 1.26e-01 | 1.063 | 0.338 |
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
| isograph | all_modules | 5 | 1.03 | 0.981 | 1.082 | 2.28e-01 | 1.028 | 0.587 |
| isograph | pheno_sig_modules | 5 | 1.04 | 0.979 | 1.105 | 2.00e-01 | 1.04 | 0.0 |
| wgcna_switch_only | pheno_sig_modules | 5 | 0.93 | 0.853 | 1.015 | 1.02e-01 | 0.92 | 0.57 |
| wgcna_multiplex | pheno_sig_modules | 5 | 1.059 | 0.998 | 1.124 | 5.97e-02 | 1.059 | 0.632 |
| isograph | go_invisible_modules | 4 | 1.072 | 0.948 | 1.211 | 2.67e-01 | 1.072 | 0.0 |
| wgcna_switch_only | go_invisible_modules | 5 | 0.929 | 0.852 | 1.014 | 9.89e-02 | 0.92 | 0.567 |
| wgcna_multiplex | go_invisible_modules | 3 | 0.996 | 0.902 | 1.099 | 9.29e-01 | 1.027 | 0.511 |
| isograph | go_visible_modules | 5 | 1.025 | 0.96 | 1.095 | 4.52e-01 | 1.025 | 0.0 |
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
