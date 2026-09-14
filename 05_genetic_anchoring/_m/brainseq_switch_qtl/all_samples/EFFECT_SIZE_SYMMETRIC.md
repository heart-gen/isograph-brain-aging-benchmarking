# Switch vs abundance QTL effect sizes, without selection asymmetry

The meta stage compares |slope| on genes significant on **both** axes. That subset is selected asymmetrically: the precise abundance phenotype clears the bar at modest effects, the noisier switch phenotype only at large ones, so the switch slopes there are inflated by a winner's curse that the abundance slopes do not share. This table uses every gene tested on both axes and measures each axis at its own lead, at the other axis's lead, and both at a common variant.

Both phenotypes are rank-INT transformed before mapping, so slopes are in SD units of the transformed phenotype and |z| is comparable across axes. `frac_switch_larger` is the fraction of genes where the switch axis is larger; Wilcoxon is paired.

| arm | region | contrast | subset | n | median S | median A | frac S > A | Wilcoxon P |
|---|---|---|---|---|---|---|---|---|
| all_samples | caudate | own_lead_z | all_tested | 11,135 | 3.745 | 4.785 | 0.232 | 0 |
| all_samples | caudate | own_lead_z | both_significant | 1,292 | 6.468 | 8.304 | 0.372 | 2.26e-30 |
| all_samples | caudate | own_lead_slope | all_tested | 11,135 | 0.433 | 0.258 | 0.742 | 0 |
| all_samples | caudate | own_lead_slope | both_significant | 1,292 | 0.502 | 0.365 | 0.687 | 5e-46 |
| all_samples | caudate | cross_lead_z | all_tested | 11,135 | 0.926 | 0.932 | 0.469 | 4.19e-28 |
| all_samples | caudate | cross_lead_z | both_significant | 1,292 | 3.785 | 4.696 | 0.361 | 1.66e-35 |
| all_samples | caudate | cross_lead_slope | all_tested | 11,135 | 0.088 | 0.057 | 0.612 | 1.72e-164 |
| all_samples | caudate | cross_lead_slope | both_significant | 1,292 | 0.266 | 0.191 | 0.614 | 3.8e-18 |
| all_samples | caudate | common_ablead_z | all_tested | 11,135 | 0.926 | 4.785 | 0.015 | 0 |
| all_samples | caudate | common_ablead_z | both_significant | 1,292 | 3.785 | 8.304 | 0.103 | 4.6e-157 |
| all_samples | caudate | common_swlead_z | all_tested | 11,135 | 3.745 | 0.932 | 0.929 | 0 |
| all_samples | caudate | common_swlead_z | both_significant | 1,292 | 6.468 | 4.696 | 0.693 | 1.44e-38 |
| all_samples | dlpfc | own_lead_z | all_tested | 11,080 | 3.742 | 4.657 | 0.236 | 0 |
| all_samples | dlpfc | own_lead_z | both_significant | 1,133 | 6.510 | 8.115 | 0.369 | 4.45e-29 |
| all_samples | dlpfc | own_lead_slope | all_tested | 11,080 | 0.479 | 0.288 | 0.740 | 0 |
| all_samples | dlpfc | own_lead_slope | both_significant | 1,133 | 0.556 | 0.415 | 0.676 | 3.72e-36 |
| all_samples | dlpfc | cross_lead_z | all_tested | 11,080 | 0.936 | 0.923 | 0.474 | 3.32e-19 |
| all_samples | dlpfc | cross_lead_z | both_significant | 1,133 | 4.131 | 4.684 | 0.350 | 9.45e-34 |
| all_samples | dlpfc | cross_lead_slope | all_tested | 11,080 | 0.099 | 0.062 | 0.610 | 1.12e-158 |
| all_samples | dlpfc | cross_lead_slope | both_significant | 1,133 | 0.315 | 0.238 | 0.597 | 9.42e-11 |
| all_samples | dlpfc | common_ablead_z | all_tested | 11,080 | 0.936 | 4.657 | 0.017 | 0 |
| all_samples | dlpfc | common_ablead_z | both_significant | 1,133 | 4.131 | 8.115 | 0.117 | 3.13e-142 |
| all_samples | dlpfc | common_swlead_z | all_tested | 11,080 | 3.742 | 0.923 | 0.932 | 0 |
| all_samples | dlpfc | common_swlead_z | both_significant | 1,133 | 6.510 | 4.684 | 0.692 | 1.32e-28 |
| all_samples | hippocampus | own_lead_z | all_tested | 10,607 | 3.713 | 4.199 | 0.293 | 0 |
| all_samples | hippocampus | own_lead_z | both_significant | 774 | 6.496 | 7.641 | 0.403 | 1.5e-13 |
| all_samples | hippocampus | own_lead_slope | all_tested | 10,607 | 0.466 | 0.198 | 0.860 | 0 |
| all_samples | hippocampus | own_lead_slope | both_significant | 774 | 0.550 | 0.307 | 0.818 | 3.01e-76 |
| all_samples | hippocampus | cross_lead_z | all_tested | 10,607 | 0.884 | 0.871 | 0.484 | 1.82e-09 |
| all_samples | hippocampus | cross_lead_z | both_significant | 774 | 4.031 | 4.784 | 0.348 | 3.92e-22 |
| all_samples | hippocampus | cross_lead_slope | all_tested | 10,607 | 0.094 | 0.043 | 0.702 | 0 |
| all_samples | hippocampus | cross_lead_slope | both_significant | 774 | 0.321 | 0.182 | 0.702 | 4.41e-39 |
| all_samples | hippocampus | common_ablead_z | all_tested | 10,607 | 0.884 | 4.199 | 0.012 | 0 |
| all_samples | hippocampus | common_ablead_z | both_significant | 774 | 4.031 | 7.641 | 0.103 | 1.21e-95 |
| all_samples | hippocampus | common_swlead_z | all_tested | 10,607 | 3.713 | 0.871 | 0.951 | 0 |
| all_samples | hippocampus | common_swlead_z | both_significant | 774 | 6.496 | 4.784 | 0.700 | 2.75e-24 |

**How to read it.** `both_significant` reproduces the meta-stage estimand and carries its selection. The unselected statistics are `own_lead_*` and `cross_lead_*` on `all_tested`; `common_ablead_z` favours abundance and `common_swlead_z` favours switch, so a real difference must survive both. Per-decile rows (`power_decile` 0–9, deciles of variants tested) are in `effect_size_symmetric.parquet`.

FDR for the `both_significant` subset: BH q < 0.05 on the permutation p-values.

## Why |slope| and |z| disagree: lead-variant allele frequency

A slope's standard error scales with 1/√(2pq), so an axis whose lead variants are rarer carries more noise in |slope|. Read |slope| as the effect estimate and |z| as the conservative statistic; this table gives the allele-frequency context for both.

| arm | region | median lead MAF S | median lead MAF A | lead MAF < 0.05 S | lead MAF < 0.05 A | median SE ratio S/A | r(|slope S|, 1/√(2pq)) |
|---|---|---|---|---|---|---|---|
| all_samples | caudate | 0.085 | 0.158 | 0.39 | 0.26 | 2.23 | 0.77 |
| all_samples | dlpfc | 0.088 | 0.162 | 0.38 | 0.25 | 2.16 | 0.82 |
| all_samples | hippocampus | 0.083 | 0.137 | 0.39 | 0.29 | 2.79 | 0.83 |
