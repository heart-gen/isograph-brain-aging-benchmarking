# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/putamen_basal_ganglia** (RSEM, n=254, 41 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 1.19904e-14 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16166 | 0 | 0.9879 | 0.09675 | 0.9854 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.25 | False | 16166 | 0 | 0.9972 | 0.0498 | 0.9952 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.5 | True | 16166 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=1 | False | 16166 | 0 | 0.996 | 0.06606 | 0.9236 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=2 | False | 16166 | 0 | 0.9805 | 0.1259 | 0.8212 | 0 | 0 | 0 | 0 |
| expression | no filter | True | 16166 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| expression | count>5,frac>=0.5 | False | 12842 | 3324 | 0.7963 | 0.3257 | 0.8164 | 0 | 0 | 0 | 0 |
| expression | count>10,frac>=0.5 | False | 11547 | 4619 | 0.7051 | 0.3546 | 0.8221 | 0 | 0 | 0 | 0 |
| expression | count>10,frac>=0.7 | False | 10105 | 6061 | 0.4103 | 0.4344 | 0.6802 | 0 | 0 | 0 | 0 |
| expression | count>20,frac>=0.7 | False | 8677 | 7489 | 0.3367 | 0.4573 | 0.6422 | 0 | 0 | 0 | 0 |
| expression | count>10,frac>=0.9 | False | 7728 | 8438 | 0.2353 | 0.4939 | 0.5017 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0 | True | 16166 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.01 | False | 14804 | 1362 | 0.9708 | 0.2149 | 0.9263 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.05 | False | 12680 | 3486 | 0.7314 | 0.3698 | 0.7074 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.1 | False | 10679 | 5487 | 0.4679 | 0.4395 | 0.5579 | 0 | 0 | 0 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1102 | 0.05485 | 0.009982 | 0.2913 | gtex | putamen_basal_ganglia |
| [3, 5) | 2459 | 0.05314 | 0.005287 | 0.2847 | gtex | putamen_basal_ganglia |
| [5, 10) | 5347 | 0.05416 | 0.007481 | 0.3226 | gtex | putamen_basal_ganglia |
| [10, 20) | 4647 | 0.05635 | 0.007747 | 0.4153 | gtex | putamen_basal_ganglia |
| [20, 10000) | 2611 | 0.06383 | 0.008809 | 0.55 | gtex | putamen_basal_ganglia |
