# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 120520.0 | 0.85 | 0.84 | 0.86 | 9.02e-153 | 0.84 | 1.85e-24 | 0.86 |
| isograph | all_modules | sQTL | 17 | 104145.0 | 0.83 | 0.82 | 0.84 | 3.00e-203 | 0.82 | 7.13e-27 | 0.89 |
| isograph | pheno_sig_modules | eQTL | 10 | 21586.0 | 0.88 | 0.86 | 0.9 | 1.28e-26 | 0.84 | 8.02e-05 | 0.9 |
| isograph | pheno_sig_modules | sQTL | 10 | 18252.0 | 0.89 | 0.87 | 0.91 | 1.65e-24 | 0.87 | 3.05e-03 | 0.93 |
| isograph | go_invisible_modules | eQTL | 7 | 5455.0 | 0.99 | 0.95 | 1.03 | 5.80e-01 | 0.98 | 5.00e-01 | 0.53 |
| isograph | go_invisible_modules | sQTL | 7 | 4971.0 | 0.97 | 0.94 | 1.01 | 1.36e-01 | 0.95 | 3.29e-01 | 0.86 |
| isograph | go_visible_modules | eQTL | 9 | 16131.0 | 0.84 | 0.82 | 0.87 | 2.79e-33 | 0.81 | 1.16e-06 | 0.84 |
| isograph | go_visible_modules | sQTL | 9 | 13281.0 | 0.87 | 0.85 | 0.89 | 2.19e-26 | 0.91 | 1.13e-01 | 0.93 |

## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 17 | 0.979 | 0.962 | 0.995 | 1.23e-02 | 0.977 | 0.839 |
| isograph | pheno_sig_modules | 10 | 1.01 | 0.978 | 1.043 | 5.60e-01 | 1.021 | 0.536 |
| isograph | go_invisible_modules | 7 | 0.981 | 0.93 | 1.035 | 4.83e-01 | 0.981 | 0.0 |
| isograph | go_visible_modules | 9 | 1.026 | 0.988 | 1.066 | 1.81e-01 | 1.094 | 0.801 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, carried by the phenotype-associated modules — splicing-QTL is spared relative to expression-QTL in the modules that track the trait.
- **Do NOT read a GO-invisible localisation off this table.** On the 2026-08-29 refresh the `go_invisible_modules` and `go_visible_modules` arms are statistically indistinguishable, so neither is the site of the effect and neither is a null. The DTU-without-DGE claim rests on the GO-invisible content gate, not on these genetics. An earlier framing made `go_invisible` the headline; it came from a stale `module_enrichment` join that relabelled the partition and must not be restored.
- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL signal is a property of the switch features (which both methods share), not of IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. IsoGraph's value is the switch-feature *representation* and its finer modules, not a unique network-inference effect.
- **Primary internal control = the matched WGCNA baselines**, not `go_visible_modules`. The baselines hold the switch features fixed and vary only the inference, so a null there localises the effect to IsoGraph's inference. `go_visible_modules` is a secondary, exploratory stratification on module CONTENT and is only ever a relative contrast — never an on/off null.
- High I2 flags between-tissue heterogeneity; prefer RE there. A nominally significant ratio carrying high I2 is driven by a few tissues, not by a consistent effect, and is weaker evidence than a smaller ratio at I2 near 0 — compare module sets on consistency as well as magnitude. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
