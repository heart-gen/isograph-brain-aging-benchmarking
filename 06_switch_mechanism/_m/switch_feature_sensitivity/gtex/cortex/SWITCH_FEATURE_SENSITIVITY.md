# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/cortex** (RSEM, n=270, 25 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12449 | 0 | 0.9955 | 0.02426 | 0.9956 | 7 | 6 | 6 | 0 |
| pseudocount | pseudocount=0.25 | False | 12449 | 0 | 0.999 | 0.01205 | 0.9991 | 6 | 6 | 6 | 0 |
| pseudocount | pseudocount=0.5 | True | 12449 | 0 | 1 | 0 | 1 | 6 | 6 | 6 | 0 |
| pseudocount | pseudocount=1 | False | 12449 | 0 | 0.9986 | 0.01526 | 0.995 | 4 | 6 | 4 | 0 |
| pseudocount | pseudocount=2 | False | 12449 | 0 | 0.9936 | 0.03542 | 0.979 | 4 | 6 | 4 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12449 | 0 | 1 | 0 | 1 | 6 | 6 | 6 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13508 | 0 | 1 | 0.1533 | 0.9612 | 4 | 6 | 4 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10264 | 2185 | 1 | 0.2034 | 0.7616 | 2 | 6 | 2 | 0 |
| expression | switching share>=0.1,frac>=0.05 | False | 13033 | 0 | 1 | 0.08418 | 0.9805 | 4 | 6 | 4 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 11179 | 1270 | 1 | 0.1539 | 0.9254 | 3 | 6 | 3 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 11008 | 2686 | 0.8499 | 0.314 | 0.6997 | 1 | 6 | 1 | 1 |
| expression | no filter | False | 16489 | 0 | 0.7554 | 0.3481 | 0.6569 | 1 | 6 | 1 | 2 |
| minor_isoform | min_usage=0 | True | 12449 | 0 | 1 | 0 | 1 | 6 | 6 | 6 | 0 |
| minor_isoform | min_usage=0.01 | False | 12449 | 0 | 1 | 0 | 1 | 6 | 6 | 6 | 0 |
| minor_isoform | min_usage=0.05 | False | 12311 | 138 | 1 | 0.03793 | 0.9964 | 6 | 6 | 6 | 0 |
| minor_isoform | min_usage=0.1 | False | 11102 | 1347 | 1 | 0.1874 | 0.8089 | 2 | 6 | 2 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4809 | 0.07957 | 0.1185 | 0.4891 | gtex | cortex |
| [3, 5) | 5507 | 0.0737 | 0.09642 | 0.4507 | gtex | cortex |
| [5, 10) | 2105 | 0.06858 | 0.06651 | 0.3625 | gtex | cortex |
| [10, 20) | 28 | 0.05533 | 0.03571 | 0.3214 | gtex | cortex |
