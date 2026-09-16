# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (graph method, module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity. Graph methods: isograph.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the qtl_anchoring arrays for each method complete).

## Pooled odds ratios (module genes vs background, per QTL type)

| graph_method | module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | eQTL | 17 | 87530.0 | 0.82 | 0.81 | 0.83 | 1.31e-192 | 0.81 | 1.18e-25 | 0.89 |
| isograph | all_modules | sQTL | 17 | 75492.0 | 0.82 | 0.81 | 0.83 | 1.82e-221 | 0.8 | 2.73e-25 | 0.91 |
| isograph | pheno_sig_modules | eQTL | 8 | 19371.0 | 0.82 | 0.8 | 0.85 | 1.71e-45 | 0.82 | 2.18e-05 | 0.9 |
| isograph | pheno_sig_modules | sQTL | 8 | 16601.0 | 0.81 | 0.79 | 0.83 | 3.67e-63 | 0.8 | 2.90e-14 | 0.79 |
| isograph | go_invisible_modules | eQTL | 8 | 5228.0 | 0.9 | 0.86 | 0.94 | 5.02e-06 | 0.9 | 2.81e-03 | 0.43 |
| isograph | go_invisible_modules | sQTL | 8 | 4736.0 | 0.91 | 0.88 | 0.95 | 1.82e-06 | 0.91 | 3.15e-02 | 0.78 |
| isograph | go_visible_modules | eQTL | 7 | 14143.0 | 0.81 | 0.78 | 0.83 | 4.14e-40 | 0.78 | 4.42e-06 | 0.89 |
| isograph | go_visible_modules | sQTL | 7 | 11865.0 | 0.79 | 0.77 | 0.82 | 8.48e-57 | 0.79 | 2.51e-11 | 0.8 |

## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)

| graph_method | module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| isograph | all_modules | 17 | 0.995 | 0.978 | 1.014 | 6.14e-01 | 0.992 | 0.78 |
| isograph | pheno_sig_modules | 8 | 0.981 | 0.946 | 1.017 | 2.99e-01 | 0.986 | 0.397 |
| isograph | go_invisible_modules | 8 | 1.0 | 0.942 | 1.061 | 9.92e-01 | 1.004 | 0.361 |
| isograph | go_visible_modules | 7 | 0.977 | 0.936 | 1.019 | 2.83e-01 | 0.989 | 0.489 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1**, carried by the phenotype-associated modules — splicing-QTL is spared relative to expression-QTL in the modules that track the trait.
- **Do NOT read a GO-invisible localisation off this table.** On the 2026-08-29 refresh the `go_invisible_modules` and `go_visible_modules` arms are statistically indistinguishable, so neither is the site of the effect and neither is a null. The DTU-without-DGE claim rests on the GO-invisible content gate, not on these genetics. An earlier framing made `go_invisible` the headline; it came from a stale `module_enrichment` join that relabelled the partition and must not be restored.
- **If the matched WGCNA baselines show the SAME specificity**, the splicing-QTL signal is a property of the switch features (which both methods share), not of IsoGraph's inference — the genetic-anchoring analog of the three-baseline result. IsoGraph's value is the switch-feature *representation* and its finer modules, not a unique network-inference effect.
- **Primary internal control = the matched WGCNA baselines**, not `go_visible_modules`. The baselines hold the switch features fixed and vary only the inference, so a null there localises the effect to IsoGraph's inference. `go_visible_modules` is a secondary, exploratory stratification on module CONTENT and is only ever a relative contrast — never an on/off null.
- High I2 flags between-tissue heterogeneity; prefer RE there. A nominally significant ratio carrying high I2 is driven by a few tissues, not by a consistent effect, and is weaker evidence than a smaller ratio at I2 near 0 — compare module sets on consistency as well as magnitude. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
