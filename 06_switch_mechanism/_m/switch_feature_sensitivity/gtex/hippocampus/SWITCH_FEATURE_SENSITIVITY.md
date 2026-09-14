# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/hippocampus** (RSEM, n=255, 41 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 9.76996e-15 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16508 | 0 | 0.9877 | 0.09723 | 0.9847 | 16 | 17 | 16 | 0 |
| pseudocount | pseudocount=0.25 | False | 16508 | 0 | 0.9972 | 0.05107 | 0.9915 | 17 | 17 | 17 | 0 |
| pseudocount | pseudocount=0.5 | True | 16508 | 0 | 1 | 0 | 1 | 17 | 17 | 17 | 0 |
| pseudocount | pseudocount=1 | False | 16508 | 0 | 0.9959 | 0.06415 | 0.9431 | 18 | 17 | 16 | 0 |
| pseudocount | pseudocount=2 | False | 16508 | 0 | 0.9804 | 0.1273 | 0.8319 | 16 | 17 | 14 | 0 |
| expression | no filter | True | 16508 | 0 | 1 | 0 | 1 | 17 | 17 | 17 | 0 |
| expression | count>5,frac>=0.5 | False | 12979 | 3529 | 0.7995 | 0.3179 | 0.7741 | 15 | 17 | 12 | 0 |
| expression | count>10,frac>=0.5 | False | 11671 | 4837 | 0.7004 | 0.3518 | 0.7168 | 14 | 17 | 11 | 2 |
| expression | count>10,frac>=0.7 | False | 10180 | 6328 | 0.4047 | 0.4332 | 0.613 | 21 | 17 | 11 | 3 |
| expression | count>20,frac>=0.7 | False | 8702 | 7806 | 0.3398 | 0.4543 | 0.5703 | 18 | 17 | 11 | 3 |
| expression | count>10,frac>=0.9 | False | 7742 | 8766 | 0.2361 | 0.4952 | 0.534 | 16 | 17 | 10 | 4 |
| minor_isoform | min_usage=0 | True | 16508 | 0 | 1 | 0 | 1 | 17 | 17 | 17 | 0 |
| minor_isoform | min_usage=0.01 | False | 15166 | 1342 | 0.9741 | 0.2107 | 0.9178 | 14 | 17 | 13 | 0 |
| minor_isoform | min_usage=0.05 | False | 13060 | 3448 | 0.7479 | 0.3586 | 0.7242 | 16 | 17 | 11 | 1 |
| minor_isoform | min_usage=0.1 | False | 11022 | 5486 | 0.4781 | 0.4358 | 0.5542 | 19 | 17 | 11 | 4 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1139 | 0.07029 | 0.06146 | 0.2599 | gtex | hippocampus |
| [3, 5) | 2525 | 0.06167 | 0.04277 | 0.2848 | gtex | hippocampus |
| [5, 10) | 5465 | 0.06663 | 0.0441 | 0.3233 | gtex | hippocampus |
| [10, 20) | 4716 | 0.07256 | 0.06891 | 0.4186 | gtex | hippocampus |
| [20, 10000) | 2663 | 0.08717 | 0.1138 | 0.5794 | gtex | hippocampus |
