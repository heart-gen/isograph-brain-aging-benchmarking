# Switch vs abundance QTL effect sizes, without selection asymmetry

The meta stage compares |slope| on genes significant on **both** axes. That subset is selected asymmetrically: the precise abundance phenotype clears the bar at modest effects, the noisier switch phenotype only at large ones, so the switch slopes there are inflated by a winner's curse that the abundance slopes do not share. This table uses every gene tested on both axes and measures each axis at its own lead, at the other axis's lead, and both at a common variant.

Both phenotypes are rank-INT transformed before mapping, so slopes are in SD units of the transformed phenotype and |z| is comparable across axes. `frac_switch_larger` is the fraction of genes where the switch axis is larger; Wilcoxon is paired.

| arm | region | contrast | subset | n | median S | median A | frac S > A | Wilcoxon P |
|---|---|---|---|---|---|---|---|---|
| ea_only | caudate | own_lead_z | all_tested | 12,688 | 3.531 | 4.025 | 0.299 | 0 |
| ea_only | caudate | own_lead_z | both_significant | 728 | 6.614 | 8.131 | 0.397 | 1.7e-12 |
| ea_only | caudate | own_lead_slope | all_tested | 12,688 | 0.543 | 0.319 | 0.753 | 0 |
| ea_only | caudate | own_lead_slope | both_significant | 728 | 0.619 | 0.441 | 0.677 | 1.45e-25 |
| ea_only | caudate | cross_lead_z | all_tested | 12,688 | 0.864 | 0.875 | 0.481 | 1.08e-12 |
| ea_only | caudate | cross_lead_z | both_significant | 728 | 4.909 | 5.475 | 0.374 | 2.26e-15 |
| ea_only | caudate | cross_lead_slope | all_tested | 12,688 | 0.117 | 0.076 | 0.601 | 2.36e-173 |
| ea_only | caudate | cross_lead_slope | both_significant | 728 | 0.437 | 0.287 | 0.609 | 1.12e-13 |
| ea_only | caudate | common_ablead_z | all_tested | 12,688 | 0.864 | 4.025 | 0.016 | 0 |
| ea_only | caudate | common_ablead_z | both_significant | 728 | 4.909 | 8.131 | 0.157 | 5.18e-75 |
| ea_only | caudate | common_swlead_z | all_tested | 12,688 | 3.531 | 0.875 | 0.949 | 0 |
| ea_only | caudate | common_swlead_z | both_significant | 728 | 6.614 | 5.475 | 0.670 | 2.24e-16 |
| ea_only | dlpfc | own_lead_z | all_tested | 12,533 | 3.526 | 3.935 | 0.317 | 0 |
| ea_only | dlpfc | own_lead_z | both_significant | 502 | 6.677 | 7.944 | 0.434 | 9.42e-07 |
| ea_only | dlpfc | own_lead_slope | all_tested | 12,533 | 0.634 | 0.384 | 0.741 | 0 |
| ea_only | dlpfc | own_lead_slope | both_significant | 502 | 0.721 | 0.560 | 0.673 | 1.35e-21 |
| ea_only | dlpfc | cross_lead_z | all_tested | 12,533 | 0.819 | 0.815 | 0.494 | 0.00174 |
| ea_only | dlpfc | cross_lead_z | both_significant | 502 | 4.782 | 5.202 | 0.420 | 4.82e-07 |
| ea_only | dlpfc | cross_lead_slope | all_tested | 12,533 | 0.132 | 0.088 | 0.605 | 4.96e-168 |
| ea_only | dlpfc | cross_lead_slope | both_significant | 502 | 0.478 | 0.347 | 0.610 | 1.28e-10 |
| ea_only | dlpfc | common_ablead_z | all_tested | 12,533 | 0.819 | 3.935 | 0.012 | 0 |
| ea_only | dlpfc | common_ablead_z | both_significant | 502 | 4.782 | 7.944 | 0.155 | 2.55e-52 |
| ea_only | dlpfc | common_swlead_z | all_tested | 12,533 | 3.526 | 0.815 | 0.964 | 0 |
| ea_only | dlpfc | common_swlead_z | both_significant | 502 | 6.677 | 5.202 | 0.683 | 7.13e-18 |
| ea_only | hippocampus | own_lead_z | all_tested | 12,488 | 3.502 | 3.769 | 0.351 | 0 |
| ea_only | hippocampus | own_lead_z | both_significant | 325 | 6.952 | 7.386 | 0.505 | 0.274 |
| ea_only | hippocampus | own_lead_slope | all_tested | 12,488 | 0.631 | 0.290 | 0.839 | 0 |
| ea_only | hippocampus | own_lead_slope | both_significant | 325 | 0.704 | 0.395 | 0.843 | 3.11e-36 |
| ea_only | hippocampus | cross_lead_z | all_tested | 12,488 | 0.786 | 0.783 | 0.498 | 0.298 |
| ea_only | hippocampus | cross_lead_z | both_significant | 325 | 5.163 | 5.420 | 0.443 | 0.00641 |
| ea_only | hippocampus | cross_lead_slope | all_tested | 12,488 | 0.131 | 0.063 | 0.678 | 0 |
| ea_only | hippocampus | cross_lead_slope | both_significant | 325 | 0.526 | 0.269 | 0.766 | 2.01e-29 |
| ea_only | hippocampus | common_ablead_z | all_tested | 12,488 | 0.786 | 3.769 | 0.011 | 0 |
| ea_only | hippocampus | common_ablead_z | both_significant | 325 | 5.163 | 7.386 | 0.212 | 3.87e-27 |
| ea_only | hippocampus | common_swlead_z | all_tested | 12,488 | 3.502 | 0.783 | 0.975 | 0 |
| ea_only | hippocampus | common_swlead_z | both_significant | 325 | 6.952 | 5.420 | 0.714 | 4.87e-16 |

**How to read it.** `both_significant` reproduces the meta-stage estimand and carries its selection. The unselected statistics are `own_lead_*` and `cross_lead_*` on `all_tested`; `common_ablead_z` favours abundance and `common_swlead_z` favours switch, so a real difference must survive both. Per-decile rows (`power_decile` 0–9, deciles of variants tested) are in `effect_size_symmetric.parquet`.

FDR for the `both_significant` subset: BH q < 0.05 on the permutation p-values.

## Why |slope| and |z| disagree: lead-variant allele frequency

A slope's standard error scales with 1/√(2pq), so an axis whose lead variants are rarer carries more noise in |slope|. Read |slope| as the effect estimate and |z| as the conservative statistic; this table gives the allele-frequency context for both.

| arm | region | median lead MAF S | median lead MAF A | lead MAF < 0.05 S | lead MAF < 0.05 A | median SE ratio S/A | r(|slope S|, 1/√(2pq)) |
|---|---|---|---|---|---|---|---|
| ea_only | caudate | 0.090 | 0.164 | 0.38 | 0.26 | 2.11 | 0.84 |
| ea_only | dlpfc | 0.092 | 0.163 | 0.38 | 0.27 | 1.99 | 0.87 |
| ea_only | hippocampus | 0.081 | 0.132 | 0.40 | 0.31 | 2.50 | 0.87 |
