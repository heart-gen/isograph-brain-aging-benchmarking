# Switch vs abundance QTL effect sizes, without selection asymmetry

The meta stage compares |slope| on genes significant on **both** axes. That subset is selected asymmetrically: the precise abundance phenotype clears the bar at modest effects, the noisier switch phenotype only at large ones, so the switch slopes there are inflated by a winner's curse that the abundance slopes do not share. This table uses every gene tested on both axes and measures each axis at its own lead, at the other axis's lead, and both at a common variant.

Both phenotypes are rank-INT transformed before mapping, so slopes are in SD units of the transformed phenotype and |z| is comparable across axes. `frac_switch_larger` is the fraction of genes where the switch axis is larger; Wilcoxon is paired.

| arm | region | contrast | subset | n | median S | median A | frac S > A | Wilcoxon P |
|---|---|---|---|---|---|---|---|---|
| ea_only | caudate | own_lead_z | all_tested | 11,133 | 3.526 | 4.028 | 0.297 | 0 |
| ea_only | caudate | own_lead_z | both_significant | 640 | 6.513 | 7.712 | 0.434 | 2.42e-07 |
| ea_only | caudate | own_lead_slope | all_tested | 11,133 | 0.536 | 0.288 | 0.779 | 0 |
| ea_only | caudate | own_lead_slope | both_significant | 640 | 0.631 | 0.411 | 0.733 | 7.18e-35 |
| ea_only | caudate | cross_lead_z | all_tested | 11,133 | 0.851 | 0.836 | 0.497 | 5.02e-05 |
| ea_only | caudate | cross_lead_z | both_significant | 640 | 4.583 | 4.956 | 0.405 | 1.06e-09 |
| ea_only | caudate | cross_lead_slope | all_tested | 11,133 | 0.112 | 0.065 | 0.635 | 1.31e-234 |
| ea_only | caudate | cross_lead_slope | both_significant | 640 | 0.404 | 0.258 | 0.673 | 4.07e-22 |
| ea_only | caudate | common_ablead_z | all_tested | 11,133 | 0.851 | 4.028 | 0.015 | 0 |
| ea_only | caudate | common_ablead_z | both_significant | 640 | 4.583 | 7.712 | 0.156 | 3.39e-65 |
| ea_only | caudate | common_swlead_z | all_tested | 11,133 | 3.526 | 0.836 | 0.953 | 0 |
| ea_only | caudate | common_swlead_z | both_significant | 640 | 6.513 | 4.956 | 0.703 | 4.39e-20 |
| ea_only | dlpfc | own_lead_z | all_tested | 11,099 | 3.528 | 3.934 | 0.317 | 0 |
| ea_only | dlpfc | own_lead_z | both_significant | 381 | 6.795 | 8.376 | 0.428 | 2.87e-06 |
| ea_only | dlpfc | own_lead_slope | all_tested | 11,099 | 0.644 | 0.350 | 0.781 | 0 |
| ea_only | dlpfc | own_lead_slope | both_significant | 381 | 0.747 | 0.542 | 0.709 | 1.36e-21 |
| ea_only | dlpfc | cross_lead_z | all_tested | 11,099 | 0.822 | 0.815 | 0.495 | 0.00242 |
| ea_only | dlpfc | cross_lead_z | both_significant | 381 | 4.964 | 5.517 | 0.420 | 1.9e-06 |
| ea_only | dlpfc | cross_lead_slope | all_tested | 11,099 | 0.131 | 0.081 | 0.625 | 2.37e-218 |
| ea_only | dlpfc | cross_lead_slope | both_significant | 381 | 0.535 | 0.375 | 0.654 | 1.98e-12 |
| ea_only | dlpfc | common_ablead_z | all_tested | 11,099 | 0.822 | 3.934 | 0.012 | 0 |
| ea_only | dlpfc | common_ablead_z | both_significant | 381 | 4.964 | 8.376 | 0.160 | 1.45e-39 |
| ea_only | dlpfc | common_swlead_z | all_tested | 11,099 | 3.528 | 0.815 | 0.964 | 0 |
| ea_only | dlpfc | common_swlead_z | both_significant | 381 | 6.795 | 5.517 | 0.661 | 2.62e-10 |
| ea_only | hippocampus | own_lead_z | all_tested | 10,605 | 3.512 | 3.763 | 0.357 | 2.65e-300 |
| ea_only | hippocampus | own_lead_z | both_significant | 297 | 6.557 | 7.800 | 0.411 | 6.08e-05 |
| ea_only | hippocampus | own_lead_slope | all_tested | 10,605 | 0.628 | 0.254 | 0.871 | 0 |
| ea_only | hippocampus | own_lead_slope | both_significant | 297 | 0.685 | 0.423 | 0.795 | 7.06e-31 |
| ea_only | hippocampus | cross_lead_z | all_tested | 10,605 | 0.797 | 0.763 | 0.500 | 0.783 |
| ea_only | hippocampus | cross_lead_z | both_significant | 297 | 5.132 | 5.497 | 0.377 | 2.41e-07 |
| ea_only | hippocampus | cross_lead_slope | all_tested | 10,605 | 0.129 | 0.056 | 0.706 | 0 |
| ea_only | hippocampus | cross_lead_slope | both_significant | 297 | 0.523 | 0.278 | 0.731 | 1.11e-22 |
| ea_only | hippocampus | common_ablead_z | all_tested | 10,605 | 0.797 | 3.763 | 0.010 | 0 |
| ea_only | hippocampus | common_ablead_z | both_significant | 297 | 5.132 | 7.800 | 0.172 | 2.09e-31 |
| ea_only | hippocampus | common_swlead_z | all_tested | 10,605 | 3.512 | 0.763 | 0.971 | 0 |
| ea_only | hippocampus | common_swlead_z | both_significant | 297 | 6.557 | 5.497 | 0.630 | 2.04e-06 |

**How to read it.** `both_significant` reproduces the meta-stage estimand and carries its selection. The unselected statistics are `own_lead_*` and `cross_lead_*` on `all_tested`; `common_ablead_z` favours abundance and `common_swlead_z` favours switch, so a real difference must survive both. Per-decile rows (`power_decile` 0–9, deciles of variants tested) are in `effect_size_symmetric.parquet`.

FDR for the `both_significant` subset: BH q < 0.05 on the permutation p-values.

## Why |slope| and |z| disagree: lead-variant allele frequency

A slope's standard error scales with 1/√(2pq), so an axis whose lead variants are rarer carries more noise in |slope|. Read |slope| as the effect estimate and |z| as the conservative statistic; this table gives the allele-frequency context for both.

| arm | region | median lead MAF S | median lead MAF A | lead MAF < 0.05 S | lead MAF < 0.05 A | median SE ratio S/A | r(|slope S|, 1/√(2pq)) |
|---|---|---|---|---|---|---|---|
| ea_only | caudate | 0.088 | 0.162 | 0.39 | 0.26 | 2.27 | 0.84 |
| ea_only | dlpfc | 0.084 | 0.160 | 0.39 | 0.27 | 2.20 | 0.88 |
| ea_only | hippocampus | 0.080 | 0.133 | 0.40 | 0.32 | 2.79 | 0.89 |
