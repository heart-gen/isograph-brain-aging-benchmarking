# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/hypothalamus** (RSEM, n=257, 27 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 1.75415e-14 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16907 | 0 | 0.9874 | 0.1035 | 0.9881 | 7 | 8 | 7 | 0 |
| pseudocount | pseudocount=0.25 | False | 16907 | 0 | 0.9971 | 0.053 | 0.9972 | 6 | 8 | 6 | 0 |
| pseudocount | pseudocount=0.5 | True | 16907 | 0 | 1 | 0 | 1 | 8 | 8 | 8 | 0 |
| pseudocount | pseudocount=1 | False | 16907 | 0 | 0.9959 | 0.06163 | 0.9147 | 6 | 8 | 6 | 0 |
| pseudocount | pseudocount=2 | False | 16907 | 0 | 0.9796 | 0.125 | 0.8097 | 5 | 8 | 4 | 0 |
| expression | no filter | True | 16907 | 0 | 1 | 0 | 1 | 8 | 8 | 8 | 0 |
| expression | count>5,frac>=0.5 | False | 13461 | 3446 | 0.8145 | 0.3145 | 0.7657 | 3 | 8 | 3 | 0 |
| expression | count>10,frac>=0.5 | False | 12151 | 4756 | 0.7102 | 0.3455 | 0.723 | 4 | 8 | 3 | 1 |
| expression | count>10,frac>=0.7 | False | 10771 | 6136 | 0.4155 | 0.4181 | 0.4422 | 2 | 8 | 1 | 1 |
| expression | count>20,frac>=0.7 | False | 9313 | 7594 | 0.3545 | 0.4384 | 0.3716 | 1 | 8 | 1 | 3 |
| expression | count>10,frac>=0.9 | False | 8358 | 8549 | 0.2484 | 0.4812 | 0.1185 | 1 | 8 | 1 | 5 |
| minor_isoform | min_usage=0 | True | 16907 | 0 | 1 | 0 | 1 | 8 | 8 | 8 | 0 |
| minor_isoform | min_usage=0.01 | False | 15520 | 1387 | 0.9714 | 0.2062 | 0.9465 | 5 | 8 | 5 | 0 |
| minor_isoform | min_usage=0.05 | False | 13389 | 3518 | 0.7257 | 0.3569 | 0.5753 | 1 | 8 | 1 | 1 |
| minor_isoform | min_usage=0.1 | False | 11322 | 5585 | 0.4628 | 0.4261 | 0.2199 | 1 | 8 | 1 | 5 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1210 | 0.0639 | 0.02727 | 0.1926 | gtex | hypothalamus |
| [3, 5) | 2625 | 0.05969 | 0.02171 | 0.1855 | gtex | hypothalamus |
| [5, 10) | 5591 | 0.06461 | 0.01914 | 0.2266 | gtex | hypothalamus |
| [10, 20) | 4793 | 0.06859 | 0.03401 | 0.3265 | gtex | hypothalamus |
| [20, 10000) | 2688 | 0.07938 | 0.04836 | 0.4888 | gtex | hypothalamus |
