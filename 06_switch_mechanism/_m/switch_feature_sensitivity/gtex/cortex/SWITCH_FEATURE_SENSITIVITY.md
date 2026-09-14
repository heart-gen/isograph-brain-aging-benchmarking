# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/cortex** (RSEM, n=270, 37 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 3.37508e-14 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16489 | 0 | 0.989 | 0.1006 | 0.9964 | 24 | 24 | 24 | 0 |
| pseudocount | pseudocount=0.25 | False | 16489 | 0 | 0.9975 | 0.05252 | 0.9982 | 24 | 24 | 24 | 0 |
| pseudocount | pseudocount=0.5 | True | 16489 | 0 | 1 | 0 | 1 | 24 | 24 | 24 | 0 |
| pseudocount | pseudocount=1 | False | 16489 | 0 | 0.9965 | 0.06016 | 0.9873 | 24 | 24 | 23 | 0 |
| pseudocount | pseudocount=2 | False | 16489 | 0 | 0.9828 | 0.1263 | 0.975 | 22 | 24 | 21 | 0 |
| expression | no filter | True | 16489 | 0 | 1 | 0 | 1 | 24 | 24 | 24 | 0 |
| expression | count>5,frac>=0.5 | False | 13474 | 3015 | 0.8355 | 0.3136 | 0.9575 | 22 | 24 | 20 | 1 |
| expression | count>10,frac>=0.5 | False | 12217 | 4272 | 0.7107 | 0.3479 | 0.9423 | 22 | 24 | 19 | 1 |
| expression | count>10,frac>=0.7 | False | 11008 | 5481 | 0.3484 | 0.4334 | 0.8051 | 9 | 24 | 8 | 4 |
| expression | count>20,frac>=0.7 | False | 9542 | 6947 | 0.2803 | 0.4526 | 0.7462 | 6 | 24 | 6 | 4 |
| expression | count>10,frac>=0.9 | False | 9078 | 7411 | 0.1975 | 0.4891 | 0.3406 | 6 | 24 | 6 | 8 |
| minor_isoform | min_usage=0 | True | 16489 | 0 | 1 | 0 | 1 | 24 | 24 | 24 | 0 |
| minor_isoform | min_usage=0.01 | False | 15198 | 1291 | 0.9785 | 0.2057 | 0.962 | 19 | 24 | 19 | 1 |
| minor_isoform | min_usage=0.05 | False | 13083 | 3406 | 0.6508 | 0.368 | 0.7844 | 11 | 24 | 10 | 4 |
| minor_isoform | min_usage=0.1 | False | 11023 | 5466 | 0.3613 | 0.4339 | 0.504 | 8 | 24 | 8 | 7 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1168 | 0.07135 | 0.06421 | 0.2269 | gtex | cortex |
| [3, 5) | 2523 | 0.06252 | 0.04043 | 0.1371 | gtex | cortex |
| [5, 10) | 5426 | 0.06081 | 0.03907 | 0.1209 | gtex | cortex |
| [10, 20) | 4710 | 0.0672 | 0.06985 | 0.1907 | gtex | cortex |
| [20, 10000) | 2662 | 0.08008 | 0.1168 | 0.2791 | gtex | cortex |
