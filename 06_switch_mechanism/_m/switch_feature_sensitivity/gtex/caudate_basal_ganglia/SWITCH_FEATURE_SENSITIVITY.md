# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/caudate_basal_ganglia** (RSEM, n=300, 42 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 8.13412e-14 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16551 | 0 | 0.9879 | 0.09885 | 0.9802 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.25 | False | 16551 | 0 | 0.9972 | 0.05353 | 0.9929 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.5 | True | 16551 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=1 | False | 16551 | 0 | 0.996 | 0.06127 | 0.9341 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=2 | False | 16551 | 0 | 0.9809 | 0.1288 | 0.8252 | 0 | 0 | 0 | 0 |
| expression | no filter | True | 16551 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| expression | count>5,frac>=0.5 | False | 13305 | 3246 | 0.8124 | 0.3099 | 0.7683 | 1 | 0 | 0 | 0 |
| expression | count>10,frac>=0.5 | False | 12025 | 4526 | 0.7069 | 0.3486 | 0.6873 | 0 | 0 | 0 | 0 |
| expression | count>10,frac>=0.7 | False | 10607 | 5944 | 0.3984 | 0.4345 | 0.5394 | 0 | 0 | 0 | 0 |
| expression | count>20,frac>=0.7 | False | 9130 | 7421 | 0.3301 | 0.4533 | 0.4897 | 0 | 0 | 0 | 0 |
| expression | count>10,frac>=0.9 | False | 8211 | 8340 | 0.2304 | 0.4912 | 0.4038 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0 | True | 16551 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.01 | False | 15223 | 1328 | 0.9727 | 0.2107 | 0.9259 | 1 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.05 | False | 13081 | 3470 | 0.7079 | 0.3678 | 0.6688 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.1 | False | 11032 | 5519 | 0.4466 | 0.4295 | 0.5512 | 0 | 0 | 0 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1149 | 0.04991 | 0.007833 | 0.2489 | gtex | caudate_basal_ganglia |
| [3, 5) | 2537 | 0.04516 | 0.001971 | 0.2755 | gtex | caudate_basal_ganglia |
| [5, 10) | 5481 | 0.0466 | 0.004014 | 0.3069 | gtex | caudate_basal_ganglia |
| [10, 20) | 4729 | 0.04791 | 0.00148 | 0.4045 | gtex | caudate_basal_ganglia |
| [20, 10000) | 2655 | 0.04811 | 0.001883 | 0.5394 | gtex | caudate_basal_ganglia |
