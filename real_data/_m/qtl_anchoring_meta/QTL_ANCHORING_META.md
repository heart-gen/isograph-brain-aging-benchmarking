# Cross-tissue meta-analysis — co-switch module xQTL anchoring

Inverse-variance meta-analysis of the matched module-membership log-OR across 17 brain region/cohort analyses, per (module set, xQTL kind). FE = fixed effect, RE = DerSimonian-Laird random effects; I2 = heterogeneity.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring_meta` (after the 13.qtl_anchoring.sh array completes).

## Pooled odds ratios (module genes vs background, per QTL type)

| module_set | xqtl_kind | k | n_fg_total | or_fe | or_fe_low | or_fe_high | p_fe | or_re | p_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_modules | eQTL | 17 | 80925.0 | 0.82 | 0.8 | 0.83 | 5.16e-122 | 0.82 | 2.74e-74 | 0.39 |
| all_modules | sQTL | 17 | 65578.0 | 0.88 | 0.86 | 0.9 | 1.63e-25 | 0.88 | 1.66e-12 | 0.55 |
| pheno_sig_modules | eQTL | 12 | 16598.0 | 0.76 | 0.73 | 0.79 | 1.21e-55 | 0.75 | 1.49e-18 | 0.64 |
| pheno_sig_modules | sQTL | 12 | 12867.0 | 0.86 | 0.82 | 0.91 | 2.18e-08 | 0.86 | 6.22e-03 | 0.69 |
| go_invisible_modules | eQTL | 10 | 8646.0 | 0.88 | 0.84 | 0.92 | 1.12e-08 | 0.89 | 9.12e-04 | 0.54 |
| go_invisible_modules | sQTL | 10 | 7197.0 | 0.99 | 0.93 | 1.05 | 7.92e-01 | 0.99 | 7.98e-01 | 0.41 |
| go_visible_modules | eQTL | 11 | 7952.0 | 0.67 | 0.64 | 0.71 | 1.22e-56 | 0.6 | 1.52e-12 | 0.84 |
| go_visible_modules | sQTL | 11 | 5670.0 | 0.72 | 0.66 | 0.78 | 1.06e-15 | 0.68 | 1.73e-07 | 0.51 |

## Splicing-specificity contrast (pooled sQTL OR / eQTL OR, paired within analysis)

| module_set | k | ratio_fe | ratio_fe_low | ratio_fe_high | p_fe | ratio_re | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all_modules | 17 | 1.071 | 1.039 | 1.104 | 7.78e-06 | 1.071 | 0.349 |
| pheno_sig_modules | 12 | 1.127 | 1.059 | 1.198 | 1.50e-04 | 1.131 | 0.222 |
| go_invisible_modules | 10 | 1.125 | 1.042 | 1.215 | 2.61e-03 | 1.12 | 0.151 |
| go_visible_modules | 11 | 1.037 | 0.942 | 1.14 | 4.62e-01 | 1.035 | 0.486 |

## Reading

- Co-switch module genes are cis-QTL **depleted** for both QTL types (OR < 1) — expected: coordinated/network genes are more constrained and carry fewer common-variant cis-QTL. This shared baseline is NOT the result.
- **The result is the splicing-specificity contrast: sQTL OR / eQTL OR > 1** in every module set (paired within analysis, removing the shared constraint baseline) => splicing-QTL is spared relative to expression-QTL in co-switch genes. The **GO-invisible disease modules** retain sQTL at background rate (OR ~ 1.0, ns) while still losing eQTL — the sharpest splicing-specific signal.
- High I2 flags between-tissue heterogeneity; prefer RE there. The contrast se is conservative (treats sQTL/eQTL estimates as independent though they share the foreground genes).
- Scope unchanged: cis-sQTL anchors member-gene splicing to genetics, not the co-switching coordination itself.
