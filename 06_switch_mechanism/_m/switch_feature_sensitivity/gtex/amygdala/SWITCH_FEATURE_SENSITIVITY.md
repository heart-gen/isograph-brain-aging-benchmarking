# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/amygdala** (RSEM, n=181, 27 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 5.77316e-14 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16359 | 0 | 0.9878 | 0.1013 | 0.9604 | 10 | 11 | 9 | 0 |
| pseudocount | pseudocount=0.25 | False | 16359 | 0 | 0.9972 | 0.05349 | 0.9685 | 8 | 11 | 8 | 0 |
| pseudocount | pseudocount=0.5 | True | 16359 | 0 | 1 | 0 | 1 | 11 | 11 | 11 | 0 |
| pseudocount | pseudocount=1 | False | 16359 | 0 | 0.996 | 0.0618 | 0.9276 | 4 | 11 | 4 | 0 |
| pseudocount | pseudocount=2 | False | 16359 | 0 | 0.9809 | 0.1272 | 0.8558 | 4 | 11 | 4 | 0 |
| expression | no filter | True | 16359 | 0 | 1 | 0 | 1 | 11 | 11 | 11 | 0 |
| expression | count>5,frac>=0.5 | False | 12786 | 3573 | 0.8077 | 0.3247 | 0.7601 | 6 | 11 | 6 | 2 |
| expression | count>10,frac>=0.5 | False | 11458 | 4901 | 0.7007 | 0.357 | 0.6958 | 1 | 11 | 1 | 2 |
| expression | count>10,frac>=0.7 | False | 10031 | 6328 | 0.4015 | 0.4343 | 0.4633 | 9 | 11 | 7 | 2 |
| expression | count>20,frac>=0.7 | False | 8579 | 7780 | 0.3229 | 0.458 | 0.3395 | 8 | 11 | 6 | 2 |
| expression | count>10,frac>=0.9 | False | 7808 | 8551 | 0.2278 | 0.4928 | 0.2066 | 15 | 11 | 7 | 3 |
| minor_isoform | min_usage=0 | True | 16359 | 0 | 1 | 0 | 1 | 11 | 11 | 11 | 0 |
| minor_isoform | min_usage=0.01 | False | 14942 | 1417 | 0.9763 | 0.2057 | 0.9165 | 3 | 11 | 3 | 0 |
| minor_isoform | min_usage=0.05 | False | 12814 | 3545 | 0.7594 | 0.3641 | 0.6079 | 1 | 11 | 1 | 2 |
| minor_isoform | min_usage=0.1 | False | 10829 | 5530 | 0.4665 | 0.4292 | 0.3754 | 11 | 11 | 6 | 2 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1122 | 0.073 | 0.04813 | 0.1569 | gtex | amygdala |
| [3, 5) | 2491 | 0.06649 | 0.03934 | 0.157 | gtex | amygdala |
| [5, 10) | 5421 | 0.06604 | 0.03486 | 0.1627 | gtex | amygdala |
| [10, 20) | 4698 | 0.07309 | 0.04257 | 0.2476 | gtex | amygdala |
| [20, 10000) | 2627 | 0.07902 | 0.0708 | 0.3856 | gtex | amygdala |
