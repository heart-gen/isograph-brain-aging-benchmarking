# Switch vs abundance QTL effect sizes, without selection asymmetry

The meta stage compares |slope| on genes significant on **both** axes. That subset is selected asymmetrically: the precise abundance phenotype clears the bar at modest effects, the noisier switch phenotype only at large ones, so the switch slopes there are inflated by a winner's curse that the abundance slopes do not share. This table uses every gene tested on both axes and measures each axis at its own lead, at the other axis's lead, and both at a common variant.

Both phenotypes are rank-INT transformed before mapping, so slopes are in SD units of the transformed phenotype and |z| is comparable across axes. `frac_switch_larger` is the fraction of genes where the switch axis is larger; Wilcoxon is paired.

| arm | region | contrast | subset | n | median S | median A | frac S > A | Wilcoxon P |
|---|---|---|---|---|---|---|---|---|
| all_samples | caudate | own_lead_z | all_tested | 12,702 | 3.746 | 4.794 | 0.237 | 0 |
| all_samples | caudate | own_lead_z | both_significant | 1,504 | 6.852 | 8.577 | 0.397 | 7.82e-25 |
| all_samples | caudate | own_lead_slope | all_tested | 12,702 | 0.449 | 0.284 | 0.715 | 0 |
| all_samples | caudate | own_lead_slope | both_significant | 1,504 | 0.509 | 0.371 | 0.665 | 1.07e-43 |
| all_samples | caudate | cross_lead_z | all_tested | 12,702 | 0.947 | 0.972 | 0.468 | 2.62e-27 |
| all_samples | caudate | cross_lead_z | both_significant | 1,504 | 4.230 | 5.021 | 0.372 | 2.77e-33 |
| all_samples | caudate | cross_lead_slope | all_tested | 12,702 | 0.093 | 0.063 | 0.594 | 8.56e-130 |
| all_samples | caudate | cross_lead_slope | both_significant | 1,504 | 0.294 | 0.213 | 0.583 | 1.14e-13 |
| all_samples | caudate | common_ablead_z | all_tested | 12,702 | 0.947 | 4.794 | 0.019 | 0 |
| all_samples | caudate | common_ablead_z | both_significant | 1,504 | 4.230 | 8.577 | 0.124 | 4.92e-170 |
| all_samples | caudate | common_swlead_z | all_tested | 12,702 | 3.746 | 0.972 | 0.928 | 0 |
| all_samples | caudate | common_swlead_z | both_significant | 1,504 | 6.852 | 5.021 | 0.686 | 2.56e-42 |
| all_samples | dlpfc | own_lead_z | all_tested | 12,537 | 3.742 | 4.671 | 0.240 | 0 |
| all_samples | dlpfc | own_lead_z | both_significant | 1,320 | 6.737 | 8.196 | 0.381 | 1.88e-25 |
| all_samples | dlpfc | own_lead_slope | all_tested | 12,537 | 0.479 | 0.321 | 0.696 | 0 |
| all_samples | dlpfc | own_lead_slope | both_significant | 1,320 | 0.553 | 0.434 | 0.636 | 7.87e-26 |
| all_samples | dlpfc | cross_lead_z | all_tested | 12,537 | 0.943 | 0.940 | 0.481 | 2.52e-17 |
| all_samples | dlpfc | cross_lead_z | both_significant | 1,320 | 4.171 | 4.823 | 0.383 | 1.37e-25 |
| all_samples | dlpfc | cross_lead_slope | all_tested | 12,537 | 0.100 | 0.070 | 0.584 | 5.99e-102 |
| all_samples | dlpfc | cross_lead_slope | both_significant | 1,320 | 0.322 | 0.232 | 0.583 | 1.82e-10 |
| all_samples | dlpfc | common_ablead_z | all_tested | 12,537 | 0.943 | 4.671 | 0.018 | 0 |
| all_samples | dlpfc | common_ablead_z | both_significant | 1,320 | 4.171 | 8.196 | 0.132 | 1.61e-152 |
| all_samples | dlpfc | common_swlead_z | all_tested | 12,537 | 3.742 | 0.940 | 0.935 | 0 |
| all_samples | dlpfc | common_swlead_z | both_significant | 1,320 | 6.737 | 4.823 | 0.711 | 4.1e-41 |
| all_samples | hippocampus | own_lead_z | all_tested | 12,552 | 3.703 | 4.228 | 0.284 | 0 |
| all_samples | hippocampus | own_lead_z | both_significant | 853 | 6.589 | 7.545 | 0.427 | 1.14e-08 |
| all_samples | hippocampus | own_lead_slope | all_tested | 12,552 | 0.472 | 0.226 | 0.825 | 0 |
| all_samples | hippocampus | own_lead_slope | both_significant | 853 | 0.544 | 0.315 | 0.782 | 9.56e-69 |
| all_samples | hippocampus | cross_lead_z | all_tested | 12,552 | 0.879 | 0.867 | 0.487 | 4.29e-08 |
| all_samples | hippocampus | cross_lead_z | both_significant | 853 | 4.472 | 5.080 | 0.373 | 6.32e-17 |
| all_samples | hippocampus | cross_lead_slope | all_tested | 12,552 | 0.095 | 0.048 | 0.675 | 0 |
| all_samples | hippocampus | cross_lead_slope | both_significant | 853 | 0.342 | 0.191 | 0.703 | 4.78e-42 |
| all_samples | hippocampus | common_ablead_z | all_tested | 12,552 | 0.879 | 4.228 | 0.014 | 0 |
| all_samples | hippocampus | common_ablead_z | both_significant | 853 | 4.472 | 7.545 | 0.137 | 2.23e-93 |
| all_samples | hippocampus | common_swlead_z | all_tested | 12,552 | 3.703 | 0.867 | 0.953 | 0 |
| all_samples | hippocampus | common_swlead_z | both_significant | 853 | 6.589 | 5.080 | 0.693 | 1.05e-28 |

**How to read it.** `both_significant` reproduces the meta-stage estimand and carries its selection. The unselected statistics are `own_lead_*` and `cross_lead_*` on `all_tested`; `common_ablead_z` favours abundance and `common_swlead_z` favours switch, so a real difference must survive both. Per-decile rows (`power_decile` 0–9, deciles of variants tested) are in `effect_size_symmetric.parquet`.

FDR for the `both_significant` subset: BH q < 0.05 on the permutation p-values.

## Why |slope| and |z| disagree: lead-variant allele frequency

A slope's standard error scales with 1/√(2pq), so an axis whose lead variants are rarer carries more noise in |slope|. Read |slope| as the effect estimate and |z| as the conservative statistic; this table gives the allele-frequency context for both.

| arm | region | median lead MAF S | median lead MAF A | lead MAF < 0.05 S | lead MAF < 0.05 A | median SE ratio S/A | r(|slope S|, 1/√(2pq)) |
|---|---|---|---|---|---|---|---|
| all_samples | caudate | 0.086 | 0.160 | 0.39 | 0.25 | 2.10 | 0.77 |
| all_samples | dlpfc | 0.087 | 0.162 | 0.38 | 0.25 | 2.00 | 0.81 |
| all_samples | hippocampus | 0.080 | 0.138 | 0.39 | 0.28 | 2.52 | 0.82 |
